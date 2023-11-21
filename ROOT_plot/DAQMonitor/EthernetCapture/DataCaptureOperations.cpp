/**
 * @file DataCaptureOperations.cpp
 *
 * @brief Top-level logic flow for the data capture module.
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <string>
#include <fstream>

#include "macros/ErrorLogger.cpp"

#include "DAQMonitor/EthernetCapture/src/PCapSessionHandler.cpp"
#include "DAQMonitor/EthernetCapture/src/NetworkDeviceException.cpp"

#include "src/ProgramControl/Terminator.cpp"
#include "src/DataModel/DAQData.cpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

namespace Muon {
namespace DataCapture {

    void runDataCapture(
        LockableStream &dataStream, 
        DAQData &data, 
        string runLabel
    );

}
namespace DataCaptureIMPL {

    bool   directoryExists             (const string       &path          );
    bool   createDirectory             (const string       &path          );
    void   createIfMissing             (const string       &directoryName );

    void   initializePCapSessionHandler(PCapSessionHandler &sessionHandler);

}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void Muon::DataCapture::runDataCapture(
    LockableStream &dataStream, 
    DAQData &data, 
    string runLabel
) {

    using namespace DataCaptureIMPL;

    // Resets all the static data stored in session handlers
    PCapSessionHandler::reset();

    PCapSessionHandler sessionHandler;
    initializePCapSessionHandler(sessionHandler);

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

    createIfMissing("./data");

    string outputFile("data/");
    outputFile += runLabel;
    outputFile += ".dat";

    ofstream fileWriter(outputFile);
    if(!fileWriter.is_open()) {

        ErrorLogger::getInstance().logError(
            string("Failed to open output file: ") + outputFile
        );
        cout << "Aborted run!" << endl;

        throw logic_error("Data capture could not open output .DAT file");

    }

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

    cout << "Saving packet data to: " << outputFile << endl;

    // TODO: Main shouldn't need to care about packet counting
    // TODO: If the data stream is a stringstream, it needs to be reset
    //       each loop, or else it will act like a memory leak
    // TODO: Place the termination condition with the thread
    int packets   = 0;
    int thousands = 0;
    while(!Terminator::getInstance().isTerminated("RUN_FLAG")) {

        packets += sessionHandler.bufferPackets(          ); // Retrieves and buffers packets from device
        sessionHandler.writePackets            (dataStream); // Writes buffered packets to dataStream
        sessionHandler.writePackets            (fileWriter); // Writes buffered packets to the .dat file
        sessionHandler.clearBuffer             (          ); // Clears the packet buffer

        // The while loop ensures that if we e.g. buffer 2000 packets at once,
        // we'll enumerate both thousands in cout.
        while(packets / 1000 > thousands) {

            ++thousands;
            cout << "Recorded " << thousands * 1000 << " packets" << endl;

        }

        data.lock();
        data.packetCount = packets;
        data.unlock();

    }

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////

    fileWriter.close();

    cout << "Suspended data capture." << endl;
    cout << data.packetCount << " packets recorded." << endl;

}

bool Muon::DataCaptureIMPL::directoryExists(const string &path) {

    struct stat sb;

    if(stat(path.data(), &sb) == 0) {

        return true;

    }

    return false;

}

bool Muon::DataCaptureIMPL::createDirectory(const string &path) {

    if(mkdir(path.data(), 0777) == 0) return true;

    return false;

}

void Muon::DataCaptureIMPL::createIfMissing(const string &directoryName) {

    if(!directoryExists(directoryName)) {

        createDirectory(directoryName);

        cout << "Created output directory: " << directoryName << endl;

    }

}

void Muon::DataCaptureIMPL::initializePCapSessionHandler(
    PCapSessionHandler &sessionHandler
) {

    State::DAQState state = State::DAQState::getState();

    try {

        sessionHandler.initializeSession(
            state.persistentState.inputDevicename
        );

        sessionHandler.setCheckPackets(true);

    } catch(NetworkDeviceException &e) {

        ErrorLogger::getInstance().logError(e.what());
        cout << "Aborted run!" << endl;

        throw UIException(
            string("Unable to listen on network device: \"" 
                + state.persistentState.inputDevicename 
                + "\". Does the device exist?")
        );

    }

}