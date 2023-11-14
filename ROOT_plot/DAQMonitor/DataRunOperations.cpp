/**
 * @file DataRunOperations.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <string>
#include <fstream>
#include <thread>
#include <sstream>

#include "macros/DAQState.cpp"
#include "macros/ErrorLogger.cpp"

#include "DAQMonitor/LockableStream.cpp"
#include "DAQMonitor/EthernetCapture/DataCaptureOperations.cpp"
#include "DAQMonitor/PacketDecoding/PacketDecodingOperations.cpp"
#include "DAQMonitor/Views/RunView.cpp"

#include "macros/UIFramework/UIException.cpp"
#include "macros/UIFramework/UISignals.cpp"

#include "src/Geometry.cpp"
#include "src/ProgramControl/Terminator.cpp"
#include "src/ProgramControl/Threads.cpp"
#include "src/DataModel/DAQData.cpp"

using namespace std;
using namespace Muon;
using namespace State;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

namespace Muon {
namespace DataRun {

    void startRun      ();
    void stopRun       ();
    void openRunViewer ();
    void closeRunViewer();

}
namespace DataRunIMPL {

    void initializeDataStream(LockableStream &dataStream);

    class RunWindows {

    public:

        RunWindows    (      RunWindows &other) = delete;
        void operator=(const RunWindows &other) = delete;

        void openRunViewer       ();
        void closeRunViewer      ();
        void closeRunViewerWindow();

        static RunWindows &getInstance();

    private:

        RunWindows();

        RunView *viewer;

    };

}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

DataRunIMPL::RunWindows::RunWindows() : viewer(nullptr) {}

void DataRunIMPL::RunWindows::openRunViewer() {

    // TODO: Open it with the currently open plotting windows selected.

    if(viewer != nullptr) return;

    viewer = new RunView(gClient->GetRoot());

    viewer->SetWindowName("Run Viewer");
    viewer->MapSubwindows();
    viewer->Resize(viewer->GetDefaultSize());
    viewer->MapWindow();

    viewer->Connect("CloseWindow()", "DataRunIMPL::RunWindows", this, "closeRunViewer()");

}

void DataRunIMPL::RunWindows::closeRunViewer() {

    if(viewer == nullptr) return;

    viewer = nullptr;

}

void DataRunIMPL::RunWindows::closeRunViewerWindow() {

    if(viewer == nullptr) return;

    RunView *temp = viewer;

    closeRunViewer();

    temp->CloseWindow();

}

DataRunIMPL::RunWindows &DataRunIMPL::RunWindows::getInstance() {

    static RunWindows instance;

    return instance;

}

void DataRunIMPL::initializeDataStream(LockableStream &dataStream) {

    dataStream.stream = nullptr;

    State::DAQState state = State::DAQState::getState();

    if(state.persistentState.dataSource == DAT_FILE_SOURCE) {

        cout << "Reading data from file: " << state.persistentState.inputFilename << endl;

        fstream *fileStream = new fstream(state.persistentState.inputFilename);
        if(!fileStream->is_open()) {
                            
            delete fileStream;
            fileStream = nullptr;

            ErrorLogger::getInstance().logError(
                string("Couldn't open file ") + state.persistentState.inputFilename
            );

            cout << "Aborted run!" << endl;

            // TODO: We need to link this exception up to an alert. But we want
            //       that alert to be created in the UI thread.
            throw UIException(
                string("File \"") 
                + state.persistentState.inputFilename 
                + "\" could not be opened."
            );

        }

        dataStream.stream = fileStream;

    } else if(state.persistentState.dataSource == NETWORK_DEVICE_SOURCE) {

        dataStream.stream = new stringstream();

    }

}

void DataRun::openRunViewer() {

    DataRunIMPL::RunWindows::getInstance().openRunViewer();

}

void DataRun::closeRunViewer() {

    DataRunIMPL::RunWindows::getInstance().closeRunViewerWindow();

}

void DataRun::stopRun() {

    using namespace DataRunIMPL;

    if(!State::DAQState::getState().tempState.runStarted) {

        throw UIException(
            "Please start a run."
        );

    }

    Terminator::getInstance().terminate("RUN_FLAG");

    // Emit a signal that the run stopped
    UISignalBus::getInstance().onRunStop();

}

// TODO: startRun should be coupled with monitor functionality, not the
//       entry view. Move the logic.
// TODO: Add controls for stopping a run without terminating GUI. We'll
//       want another kind of termination signal. Perhaps we can expand
//       the Terminator to have multiple kinds of signals our threads
//       can condition on.
void DataRun::startRun() {

    using namespace DataRunIMPL;

    // TODO: Since ErrorLogger is application-wide, we don't want to do
    //       this. But we also need to count per-run errors.
    ErrorLogger::getInstance().clear();

    DAQState state = DAQState::getState();

    if(state.tempState.runStarted) {

        throw UIException(
            "Please finish the current run before starting a new run."
        );

    }

    if(state.persistentState.dataSource == DAT_FILE_SOURCE) {

        if(state.persistentState.inputFilename == "") {

            throw UIException("Please select a DAT file.");

        }

        cout << endl << "File source selected" << endl;
        cout << "Filename: " << state.persistentState.inputFilename << endl;

    } else if (state.persistentState.dataSource == NETWORK_DEVICE_SOURCE) {

        if(state.persistentState.inputDevicename == "") {

            throw UIException("Please select a network device.");

        }

        cout << endl << "Network source selected" << endl;
        cout << "Ethernet device: " << state.persistentState.inputDevicename << endl;

    }

    state.tempState.runStarted = true;
    state.commit(); // NOTE: This shouldn't fail, but better if it's robust

    // Clear the DAQData of any data from a previous run
    DAQData &data = DAQData::getInstance();

    data.lock  ();
    data.clear ();
    data.unlock();

    Geometry::getInstance().SetRunN(getRunNumber());

    // TODO: Hook up error handling on a per-thread basis. Threads should
    //       report to a threadsafe error handler that does the error handling
    // TODO: It's not obvious why we make a parent thread. It's because we want
    //       a few things to happen after the capture and decode threads are
    //       joined. Make that more obvious.
    ProgramFlow::threadLock.lock();
    ProgramFlow::threads.emplace_back(thread([&data]() {

        LockableStream dataStream;
        initializeDataStream(dataStream);

        // TODO: Put the thread termination conditions here
        thread dataCaptureThread([&dataStream, &data]() {

            if(DAQState::getState().persistentState.dataSource == NETWORK_DEVICE_SOURCE) {

                DataCapture::runDataCapture(dataStream, data);

            }

        });

        thread decodeThread([&dataStream, &data](){

            Decode::runDecoding(dataStream, data);

        });

        dataCaptureThread.join();
        decodeThread     .join();

        // TODO: Again, I would rather avoid caring about the type of stream.
        // TODO: We don't really need to lock here
        dataStream.lock();
        fstream *temp = dynamic_cast<fstream*>(dataStream.stream);
        if(temp) {
            temp->close();
        }
        dataStream.unlock();

        if(dataStream.stream) delete dataStream.stream;
        dataStream.stream = nullptr;

        DAQState state = DAQState::getState();
        state.tempState.runStarted = false;
        state.commit();

        // Once the run is shut down, we can clear the
        // run flag if it exists
        Terminator::getInstance().clearFlag("RUN_FLAG");

        cout << "Run finished!" << endl;

    }));

    ProgramFlow::threadLock.unlock();

    UISignalBus::getInstance().onRunStart();

}