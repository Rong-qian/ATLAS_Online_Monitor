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

#include "macros/Monitor.cpp"
#include "macros/ErrorLogger.cpp"

#include "src/Geometry.cpp"
#include "src/EthernetCapture/DeviceSelector.cpp"
#include "src/EthernetCapture/PCapSessionHandler.cpp"
#include "src/ProgramControl/Terminator.cpp"
#include "src/ProgramControl/SigHandlers.cpp"
#include "src/DataModel/DAQData.cpp"

#include "monitorConfig.cpp"

using namespace std;
using namespace Muon;

///////////////////////////////////////////////////////////////////////////////
////////////////////////////// INTERFACE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// TODO: Rethink the semantics of Monitor. We've made our labels as if it's a
//       display element when really it does data decoding. It should be the
//       UI element that does the 'refreshing', for example.

// TODO: Examine this for ROOT tips:
//       https://mightynotes.wordpress.com/2020/02/15/cern-root-tips-and-tricks/

// TODO: Call getpid() for the run? Do we need to give it a pid?

// TODO: Config GUI for options like whether to show lost packets

// TODO: Split out all cout calls to a console logger object. Easy to make it
//       threadsafe or switch the stream we're logging to

// TODO: Make the error logger threadsafe

// TODO: Set Geometry::runN as in DecodeOffline.cpp

/**
 * Macro defining the entry command for the monitor.
 */
void StartMonitor(const string &filename = "");

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

///////////////////////////////////////////////////////////////////////////////
/////////////////////////// IMPLEMENTATION ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void StartMonitor(const string &filename = "") {

	ErrorLogger::getInstance().setOutputStream(cerr);

	///////////////////////////////////////////////////////////////////////////
	///////////////////////// DATA STREAM SETUP ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	MonitorHooks::beforeStartRun();

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

	DeviceManager      devices       ;
	PCapSessionHandler sessionHandler;

	if(filename == "") {

		try {

			devices.initialize();

			PCapDevice networkDevice = runDeviceSelector(devices);

			sessionHandler.initializeSession(networkDevice);

			 // Sets the handler to notify the user when a packet is lost
			sessionHandler.setCheckPackets(true);

		} catch(NetworkDeviceException &e) {

			ErrorLogger::getInstance().logError(e.what());
			cout << "Aborted run!" << endl;
			return 1; // NOTE: We don't call the terminator here because we
			          //       don't need to close down any threads and we
			          //       want the program to exit immediately.

		}

	}

	// THIS MUST BE CALLED BEFORE STARTING ANY THREADS.
	// It intercepts SIGINT/SIGTERM/SIGQUIT to cleanly terminate threads.
	// It must also be called after we're done getting user input or we'll
	// get weird behavior where code will run expecting data that does not
	// exist, since it does not interrupt anything
	setTerminationHandlers(flagForTermination);

	// TODO: Make session handler non-blocking
	//         -- thread it and implement a thread-safe data stream object
	// TODO: Make data capture run on the main thread to avoid lots of problems
	// TODO: Make sure we make sessionHandler stop waiting for a packet if we
	//       ctrl+c
	// TODO: This lambda is a bit cluttered.
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

		// TODO: Main shouldn't need to care about packet counting
		// TODO: If the data stream is a stringstream, it needs to be reset
		//       each loop, or else it will act like a memory leak
		int i = 0;
		while(!Terminator::getInstance().isTerminated()) {

			sessionHandler.bufferPackets(          ); // Retrieves and buffers packets from device
			sessionHandler.writePackets (dataStream); // Writes buffered packets to dataStream
			sessionHandler.writePackets (fileWriter); // Writes buffered packets to the .dat file
			sessionHandler.clearBuffer  (          ); // Clears the packet buffer

			++i;

			if(i % 1000 == 0) {
				cout << "Recorded " << i << " packets" << endl; // TODO: mutex
			}

		}

		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////

		fileWriter.close();

		cout << endl << "Data capture finished!" << endl;
		cout << i << " packets recorded." << endl;

	});

	///////////////////////////////////////////////////////////////////////////
	/////////////////////////// DATA PROCESSING ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	thread decodeThread([&dataStream]() {

		Geometry::getInstance().SetRunN(getRunNumber());

		Monitor monitor(dataStream);

		while(!Terminator::getInstance().isTerminated()) {

			MonitorHooks::beforeUpdateData();

			// TODO: Performance analysis. I'd like this loop to run faster
			//         -- I think binning and drawing is our weak point. Let's
			//            bin every event before drawing

			// FIXME: In file reading mode, this will read the whole file before
			//        terminating on ctrl+c

			monitor.refresh();

			// TODO: This is hacky; fix it. The idea here is to clear processed
			//       data from the dataStream so we don't produce a de facto
			//       memory leak. But we'd rather the code not have to care what
			//       kind of stream dataStream is. Perhaps we make our own kind
			//       of iostream that clears after read?
			//       https://stackoverflow.com/questions/63034484/how-to-create-stream-which-handles-both-input-and-output-in-c
			//       https://stackoverflow.com/questions/12410961/c-connect-output-stream-to-input-stream
			//       https://stackoverflow.com/questions/26346320/how-to-redirect-input-stream-to-output-stream-in-one-line
			// TODO: We might be able to hook our file and data output streams
			//       together so we only have to write to one of them:
			//       https://stackoverflow.com/questions/1760726/how-can-i-compose-output-streams-so-output-goes-multiple-places-at-once
			// TODO: Might it make sense to make the data stream unbuffered? See:
			//       https://stackoverflow.com/questions/52581080/usage-of-output-stream-buffer-in-context-to-stdcout-and-stdendl
			// TODO: Anyway, we can revisit how we want to handle the data streams
			//       later. For now, this will suffice.
			dataStream.lock();
			stringstream *temp = dynamic_cast<stringstream*>(dataStream.stream);
			if(temp) {
				string unread = temp->eof() ?
					"" : temp->str().substr(temp->tellg());
				temp->str(unread);
			}
			dataStream.unlock();

			MonitorHooks::updatedData();

			this_thread::sleep_for(chrono::milliseconds((int)(1000 / DATA_REFRESH_RATE)));

		}

		cout << "Suspended data decoding." << endl; // TODO: mutex

	});

	MonitorHooks::startedRun();

	decodeThread     .join();
	dataCaptureThread.join();

	MonitorHooks::finishedRun();

	cout << endl << "Processed " << DAQData::getInstance().processedEvents.size() << " events." << endl;

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

	cout << "Shut down complete!" << endl;

	gApplication->Terminate(0);

}