/**
 * @file StartMonitor.cpp
 *
 * @brief Entry point for the online monitor.
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 * 
 * NOTE: This assumes triggerless mode.
 */

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// DEPENDENCIES /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <stdio.h>
#include <thread>
#include <sstream>
#include <fstream>

#include <sys/stat.h>

#include "analysis/MonitorHooks.cpp"

#include "macros/ErrorLogger.cpp"

#include "DAQMonitor/PacketDecoding/src/Decoder.cpp"
//#include "DAQMonitor/EthernetCapture/src/DeviceSelector.cpp"
//#include "DAQMonitor/EthernetCapture/src/PCapSessionHandler.cpp"

#include "src/Geometry.cpp"
#include "src/ProgramControl/Terminator.cpp"
#include "src/ProgramControl/SigHandlers.cpp"
#include "src/DataModel/DAQData.cpp"

#include "monitorConfig.cpp"

using namespace std ;
using namespace Muon;

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// INTERFACE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * Macro defining the entry command for the monitor.
 */
void   StartMonitor       (const string &filename = "");

bool   directoryExists    (const string &path         );
bool   createDirectory    (const string &path         );
string getCurrentTimestamp(const string &format       );
void   createIfMissing    (const string &directoryName);

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// IMPLEMENTATION ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void StartMonitor(const string &filename = "") {

	ErrorLogger::getInstance().setOutputStream(cerr);

	///////////////////////////////////////////////////////////////////////////
	///////////////////////// DATA STREAM SETUP ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	LockableStream dataStream;
	dataStream.stream = nullptr;

	if(filename == "") {

		dataStream.stream = new stringstream();

	} else {

		cout << "Reading data from file: " << filename << endl;

		fstream *fileStream = new fstream(filename);
		if(!fileStream->is_open()) {
			
			delete fileStream;
			fileStream = nullptr;

			ErrorLogger::getInstance().logError(
				string("Couldn't open file ") + filename
			);

			cout << "Aborted run!" << endl;
			exit(1);

		}

		dataStream.stream = fileStream;

	}

	///////////////////////////////////////////////////////////////////////////
	//////////////////////////// DATA CAPTURE /////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	//DeviceManager      devices       ;
	//PCapSessionHandler sessionHandler;

	if(filename == "") {
		/*
		try {

			devices.initialize();

			PCapDevice networkDevice = runDeviceSelector(devices);

			sessionHandler.initializeSession(networkDevice);

			 Sets the handler to notify the user when a packet is lost
			sessionHandler.setCheckPackets(true);

		} catch(NetworkDeviceException &e) {

			ErrorLogger::getInstance().logError(e.what());
			cout << "Aborted run!" << endl;
			return 1; // NOTE: We don't call the terminator here because we
			          //       don't need to close down any threads and we
			          //       want the program to exit immediately.

		}
		*/
	}

	// THIS MUST BE CALLED BEFORE STARTING ANY THREADS.
	// It intercepts SIGINT/SIGTERM/SIGQUIT to cleanly terminate threads.
	// It must also be called after we're done getting user input or we'll
	// get weird behavior where code will run expecting data that does not
	// exist, since it does not interrupt anything
	setTerminationHandlers(flagForTermination);

	DAQData &data = DAQData::getInstance();

	MonitorHooks::beforeStartRun(data);
	/*
	thread dataCaptureThread([&sessionHandler, &dataStream]() {

		// End the thread if the session handler isn't ready, which means
		// we're reading from a file.
		if(!sessionHandler.isReady()) return;

		string runTimestamp = getCurrentTimestamp("%Y%m%d_%H%M%S");

		cout << endl << "Starting run: " << runTimestamp << endl; // TODO: Add the run number to this

		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////

		createIfMissing("./data");

		string outputFile("data/run_");
		outputFile += runTimestamp;
		outputFile += ".dat";

		ofstream fileWriter(outputFile);
		if(!fileWriter.is_open()) {

			ErrorLogger::getInstance().logError(
				string("Failed to open output file: ") + outputFile
			);
			cout << "Aborted run!" << endl;
			Terminator::getInstance().terminate();
			return;

		}
	
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////

		cout << "Saving packet data to: " << outputFile << endl;

		int packets   = 0;
		int thousands = 0;
		while(!Terminator::getInstance().isTerminated()) {

			packets += sessionHandler.bufferPackets(          ); // Retrieves and buffers packets from device
			sessionHandler.writePackets            (dataStream); // Writes buffered packets to dataStream
			sessionHandler.writePackets            (fileWriter); // Writes buffered packets to the .dat file
			sessionHandler.clearBuffer             (          ); // Clears the packet buffer

			while(packets / 1000 > thousands) {

				++thousands;
				cout << "Recorded " << thousands * 1000 << " packets" << endl;

			}

		}

		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////

		fileWriter.close();

		cout << endl << "Data capture finished!" << endl;
		cout << packets << " packets recorded." << endl;

	});
	*/
	///////////////////////////////////////////////////////////////////////////
	/////////////////////////// DATA PROCESSING ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	thread decodeThread([&dataStream, &data]() {

		Geometry::getInstance().SetRunN(getRunNumber());

		Decoder decoder(dataStream, DAQData::getInstance());

		while(!Terminator::getInstance().isTerminated()) {

			MonitorHooks::beforeUpdateData(data);

			decoder.refresh();

			// Hacky fix for clearing the stringstream every loop. 
			// TODO: Avoid the dynamic cast
			dataStream.lock();
			stringstream *temp = dynamic_cast<stringstream*>(dataStream.stream);
			if(temp) {
				string unread = temp->eof() ?
					"" : temp->str().substr(temp->tellg());
				temp->str(unread);
			}
			dataStream.unlock();

			MonitorHooks::updatedData(data);

			this_thread::sleep_for(chrono::milliseconds((int)(1000 / DATA_REFRESH_RATE)));

		}

		cout << "Suspended data decoding." << endl;

	});

	MonitorHooks::startedRun(data);

	decodeThread     .join();
	//dataCaptureThread.join();

	MonitorHooks::finishedRun(data);

	cout << endl << "Processed " << DAQData::getInstance().processedEvents.size() << " events." << endl;

	dataStream.lock();
	fstream *temp = dynamic_cast<fstream*>(dataStream.stream);
	if(temp) {
		temp->close();
	}
	dataStream.unlock();

	if(dataStream.stream) delete dataStream.stream;
	dataStream.stream = nullptr;

	cout << "Shut down complete!" << endl;

	gApplication->Terminate(0);

}

bool directoryExists(const string &path) {

	struct stat sb;

	if(stat(path.data(), &sb) == 0) {

		return true;

	}

	return false;

}

bool createDirectory(const string &path) {

	if(mkdir(path.data(), 0777) == 0) return true;

	return false;

}

string getCurrentTimestamp(const string &format) {

	char formatBuffer[40];
	time_t sys_time;
	struct tm *timeinfo;
	sys_time = time(0);
	timeinfo = localtime(&sys_time);
	memset(formatBuffer, 0, sizeof(formatBuffer));
	strftime(formatBuffer, 40, format.data(), timeinfo);

	return string(formatBuffer);

}

void createIfMissing(const string &directoryName) {

	if(!directoryExists(directoryName)) {

		createDirectory(directoryName);

		cout << "Created output directory: " << directoryName << endl;

	}

}