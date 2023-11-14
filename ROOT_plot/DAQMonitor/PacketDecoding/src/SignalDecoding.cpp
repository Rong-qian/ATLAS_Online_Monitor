/**
 * @file SignalDecoding.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <stdio.h>
#include <string>
#include <istream>

#include "macros/ErrorLogger.cpp"

#include "DAQMonitor/LockableStream.cpp"

#include "src/Signal.cpp"
#include "src/Geometry.cpp"

using namespace Muon;
using namespace std;

// TODO: This does three different things. Break it up.

const string SIGNAL_ERROR = "signalErr" ;
const string SIGNAL_WARN  = "signalWarn";

// Byte swap from big-endian to little-endian or vice versa
uint64_t byteSwap(uint64_t data, uint8_t dataSize);

/**
 * Reads signals on command from a filestream.
 */
class SignalReader {

public:

	SignalReader(LockableStream &in);

	bool   isReady      ();
	Signal extractSignal();

private:

	LockableStream &dataStream;

};

/**
 * Provides validation of 'sig', logging an error to the error logger
 * and returning false should sig produce an error.
 * 
 * Signals that fail error validation should be dropped.
 * 
 * @param sig The signal to validate.
 */
bool validateSignalErrors(const Signal &sig);

/**
 * Provides validation of 'sig', logging any warnings to the error logger.
 * 
 * Signals that fail warning validation should not be dropped.
 * 
 * @param sig The signal to validate.
 */
void validateSignalWarnings(const Signal &sig);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

SignalReader::SignalReader(LockableStream &in) :
	dataStream(in) {}

bool SignalReader::isReady() {

	// TODO: What if no dataStream.end?

	// Get initial and end positions for the stream
	dataStream.lock();
	streampos initialPos = dataStream.stream->tellg();
	dataStream.stream->seekg(0, dataStream.stream->end);
	streampos endPos = dataStream.stream->tellg();

	// Reset the stream
	dataStream.stream->clear();
	dataStream.stream->seekg(initialPos, dataStream.stream->beg);
	dataStream.unlock();

	// Return true if there are at least WORD_SIZE characters left on the
	// stream
	return (endPos - initialPos) >= Signal::WORD_SIZE;

}

// TODO: It's faster if we read all available signals at once and parse them
//       into signals later
// TODO: Remember the power of memcpy and static casting for parsing things
//         -- build the right struct and you can memcpy right into it very fast
Signal SignalReader::extractSignal() {

	char readBuffer[Signal::WORD_SIZE];

	dataStream.lock();
	dataStream.stream->read(readBuffer, Signal::WORD_SIZE);
	dataStream.unlock();

	uint64_t word = 0;
	memcpy(&word, readBuffer, Signal::WORD_SIZE);

	// TODO: Consider moving the byteswap to PCapSessionhandler?
	word = byteSwap(word, Signal::WORD_SIZE);

	return Signal(word);

}

// TODO: Pull out validators
// TODO: Test validators
// TODO: We don't have event numbers, which we need to have for error reporting
// TODO: Combine this with signal buffer validation. Either this can be done
//       with access to the buffer or buffer validation can be done here.
//       If it happens in buffer validation, we can still keep track of signal
//       count.
bool validateSignalErrors(const Signal &sig) {

	if(sig.isEventHeader ()) return true;
	if(sig.isEventTrailer()) return true;

	ErrorLogger &logger = ErrorLogger::getInstance();
	Geometry    &geo    = Geometry::getInstance   ();

	if(!geo.IsActiveTDC(sig.TDC())) {

		string msg = "ERROR -- Unexpected data TDCID = ";
		msg += to_string(sig.TDC());
		msg += ", Channel = ";
		msg += to_string(sig.Channel());

		logger.logError(msg, SIGNAL_ERROR);

		return false;

	}

	if(sig.isTDCError()) {

		string msg = "";

		TDCErrorData errorData = sig.getTDCError();

		msg = "ERROR -- Received TDC Error Signal:";

		if(errorData.LSBFlag1 > 0) {

			msg += "\n         \tTDCID = ";
			msg += to_string(errorData.TDC);
			msg += ", Channel = ";
			msg += to_string(errorData.LSBChannel1);
			msg += " overflowed!";

		} if(errorData.LSBFlag2 > 0) {

			msg += "\n         \tTDCID = ";
			msg += to_string(errorData.TDC);
			msg += ", Channel = ";
			msg += to_string(errorData.LSBChannel2);
			msg += " overflowed!";

		} 

		logger.logError(msg, SIGNAL_ERROR);

		return false;

	}

	if(sig.isTDCHeader()) {

		// TODO: Validation

		return true;

	}

	if(sig.isTDCTrailer()) {

		// TODO: Validation

		return true;

	}

	if(!geo.IsActiveTDCChannel(sig.TDC(), sig.Channel())) {

		string msg = "ERROR -- Unexpected data TDCID = ";
		msg += to_string(sig.TDC());
		msg += ", Channel = ";
		msg += to_string(sig.Channel());

		logger.logError(msg, SIGNAL_ERROR);

		return false;

	}

	return true;

}

void validateSignalWarnings(const Signal &sig) {



}

uint64_t byteSwap(uint64_t data, uint8_t dataSize) {

	uint64_t result = 0;
	for(uint8_t byte = 0; byte < dataSize; ++byte) {

		uint64_t temp = data >> 8 * (dataSize - byte - 1);
		temp = temp & 0xff;

		result += temp << 8 * byte;

	}

	return result;

}