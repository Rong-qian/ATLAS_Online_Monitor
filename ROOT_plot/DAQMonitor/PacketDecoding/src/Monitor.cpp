/**
 * @file Monitor.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

///////////////////////////////////////////////////////////////////////////////
//////////////////////////// DEPENDENCIES /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include <istream>

#include "macros/ErrorLogger.cpp"

#include "DAQMonitor/PacketDecoding/src/SignalDecoding.cpp"
#include "DAQMonitor/PacketDecoding/src/EventDecoding.cpp"
#include "DAQMonitor/LockableStream.cpp"

#include "src/Signal.cpp"
#include "src/Event.cpp"
#include "src/DataModel/DAQData.cpp"

///////////////////////////////////////////////////////////////////////////////
//////////////////////// INTERFACE ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// TODO: Preserve ability to read from a file source

// TODO: Let the monitor just take in a data stream, and choose between an
//       ifstream and another kind of stream at a higher level.

/**
 * Provides state and functionality for the online monitor.
 */
class Monitor {

public:

	Monitor(LockableStream &in, DAQData &dataOut);

	bool isStale();
	void refresh();

private:

	SignalReader reader;

	DAQData &data;

	vector<Signal> signalBuffer;
	vector<Event > eventBuffer ;

};

///////////////////////////////////////////////////////////////////////////////
//////////////////////// HELPER FUNCTIONS /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/**
 * Returns true if signals represents a complete event (i.e. if signals ends
 * with a trailer).
 */
bool isEvent(const vector<Signal> &signals);

///////////////////////////////////////////////////////////////////////////////
/////////////////////// IMPLEMENTATION ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Monitor::Monitor(LockableStream &in, DAQData &dataOut) 
	: reader(in), data(dataOut) {}

bool Monitor::isStale() {

	return reader.isReady();

}

void Monitor::refresh() {

	while(isStale()) {

		Signal sig = reader.extractSignal();

		if(validateSignalErrors(sig)) {

			// TODO: We should keep some metadata if it's a TDC header or
			//       trailer. See DecodeOffline.cpp

			signalBuffer.push_back(sig);

		} else {

			cerr << "Dropped signal" << endl;

		}

		validateSignalWarnings(sig);

		if(isEvent(signalBuffer)) {

			Event e = assembleEvent(signalBuffer);

			if(validateEventErrors(e)) {

				eventBuffer.push_back(e);

			} else {

				cerr << "Dropped event" << endl;

			}

			signalBuffer.clear();

			validateEventWarnings(e);

		}

	}

	data.lock();
	data.newEvents.clear();
	data.unlock();

	for(Event &e : eventBuffer) {

		// Condition ignores empty events
		if(e.Trailer().HitCount() != 0) {

			processEvent(e);

			data.lock();
			data.processedEvents.push_back(e);
			data.newEvents.push_back(e);
			data.plots.binEvent(e);
			data.unlock();

		}

	}

	// TODO: Make it clear that totalEventCount does not include events that
	//       were dropped
	data.lock();
	data.totalEventCount += eventBuffer.size();
	data.unlock();

	eventBuffer.clear();

	// NOTE: processedEvents will maintain a record of each event received
	//       by the monitor

}

bool isEvent(const vector<Signal> &signals) {

	if(signals.empty()) return false;

	return signals.back().isEventTrailer();

}