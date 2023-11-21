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

#include <sys/stat.h>

#include "macros/DAQState.cpp"
#include "macros/ErrorLogger.cpp"

#include "DAQMonitor/LockableStream.cpp"
#include "DAQMonitor/EthernetCapture/DataCaptureOperations.cpp"
#include "DAQMonitor/PacketDecoding/PacketDecodingOperations.cpp"

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

}
namespace DataRunIMPL {

    void   initializeDataStream(LockableStream &dataStream);
    string getCurrentTimestamp (const string   &format    );

}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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

    string runLabel = "";

    if(state.tempState.runStarted) {

        throw UIException(
            "Please finish the current run before starting a new run."
        );

    }

    if(state.persistentState.dataSource == DAT_FILE_SOURCE) {

        if(state.persistentState.inputFilename == "") {

            throw UIException("Please select a DAT file.");

        }

        string filename = state.persistentState.inputFilename;

        cout << endl << "File source selected" << endl;
        cout << "Filename: " << filename << endl;

        // TODO: Extract timestamp

        size_t extensionPos = filename.find_last_of(".");

        if(extensionPos == string::npos) {

            runLabel = filename;

        } else {

            runLabel = filename.substr(0, extensionPos);

        }

    } else if (state.persistentState.dataSource == NETWORK_DEVICE_SOURCE) {

        if(state.persistentState.inputDevicename == "") {

            throw UIException("Please select a network device.");

        }

        cout << endl << "Network source selected" << endl;
        cout << "Ethernet device: " << state.persistentState.inputDevicename << endl;

        runLabel = string("run_") + getCurrentTimestamp("%Y%m%d_%H%M%S");

    }

    state.tempState.runStarted = true;
    state.tempState.runLabel = runLabel;
    state.commit(); // NOTE: This shouldn't fail, but better if it's robust

    // Clear the DAQData of any data from a previous run
    DAQData &data = DAQData::getInstance();

    data.lock  ();
    data.clear ();
    data.unlock();

    // TODO: Hook up error handling on a per-thread basis. Threads should
    //       report to a threadsafe error handler that does the error handling
    // TODO: It's not obvious why we make a parent thread. It's because we want
    //       a few things to happen after the capture and decode threads are
    //       joined. Make that more obvious.
    ProgramFlow::threadLock.lock();
    ProgramFlow::threads.emplace_back(thread([&data, runLabel]() {

        // TODO: Add the run number to this
        cout << endl << "Starting run: " << runLabel << endl; 

        // NOTE: Following the legacy code, runN is in YYYYMMDD format and does
        //       not include hours/minutes/seconds
        // NOTE: This assumes the DAT filenames are formatted as "run_YYYYMMDD_HHMMSS.dat"
        int runN = (
            (TObjString*)(TString(
                runLabel.substr(3, runLabel.size()).data()
            ).Tokenize("_")->At(0))
        )->String().Atoi();
        Geometry::getInstance().SetRunN(runN);

        LockableStream dataStream;
        initializeDataStream(dataStream);

        // TODO: Put the thread termination conditions here
        thread dataCaptureThread([&dataStream, &data, runLabel]() {

            if(DAQState::getState().persistentState.dataSource == NETWORK_DEVICE_SOURCE) {

                DataCapture::runDataCapture(dataStream, data, runLabel);

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

string Muon::DataRunIMPL::getCurrentTimestamp(const string &format) {

    char formatBuffer[40];
    time_t sys_time;
    struct tm *timeinfo;
    sys_time = time(0);
    timeinfo = localtime(&sys_time);
    memset(formatBuffer, 0, sizeof(formatBuffer));
    strftime(formatBuffer, 40, format.data(), timeinfo);

    return string(formatBuffer);

}