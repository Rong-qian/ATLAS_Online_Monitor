/**
 * @file EventDecoding.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <vector>
#include <string>
#include <algorithm>

#include "macros/ErrorLogger.cpp"

#include "src/Signal.cpp"
#include "src/Event.cpp"
#include "src/HitFinder.cpp"
#include "src/TimeCorrection.cpp"

const string EVENT_ERROR = "eventErr" ;
const string EVENT_WARN  = "eventWarn";

Event assembleEvent        (      vector<Signal>  signals);
bool  validateEventErrors  (const Event          &e      );
void  validateEventWarnings(const Event          &e      );
void  processEvent         (      Event          &e      );

// TODO: This does three different things. Break it up.

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// TODO: We also want some metadata, as described in DecodeOffline
Event assembleEvent(vector<Signal> signals) {

	return Event(
		signals.front(), 
		signals.back(),
		vector<Signal>(
			signals.begin() + 1, 
			signals.end  () - 1
		)
	);

}

// TODO: Separate errors that we drop the event for and errors that we 
//       report without dropping the event
//         -- validateEventErrors and validateEventWarnings
// TODO: Use exceptions?
// Validation that, when failed, drops the event
bool validateEventErrors(const Event &e) {

	ErrorLogger &logger = ErrorLogger::getInstance();

	if(!e.Header().isEventHeader()) {

		logger.logError(
			"ERROR -- Found event with no header",
			EVENT_ERROR
		);

		return false;

	}

	if(!e.Trailer().isEventTrailer()) {

		logger.logError(
			"ERROR -- Found event with no trailer",
			EVENT_ERROR
		);

		return false;

	}

	if(e.Header().HeaderEID() != e.Trailer().TrailerEID()) {

		logger.logError(
			"ERROR -- Found event with mismatched event IDs",
			EVENT_ERROR
		);

	}

	for(const Signal &sig : e.Signals()) {

		if(sig.isEventHeader()) {

			logger.logError(
				"ERROR -- Found event with multiple headers",
				EVENT_ERROR
			);

			return false;

		}

		if(sig.isEventTrailer()) {

			logger.logError(
				"ERROR -- Found event with multiple trailers",
				EVENT_ERROR
			);

			return false;

		}

	}

	return true;

}

// Validation that, when failed, keeps the event but warns the user.
void validateEventWarnings(const Event &e) {

	ErrorLogger &logger = ErrorLogger::getInstance();

	if(e.Trailer().HeaderCountErr()) {

		logger.logError(
			string("WARNING -- Header count error flag. Got ")
			+ to_string(e.Trailer().TDCHdrCount())
			+ " header(s)!",
			EVENT_WARN
		);

	}

	if(e.Trailer().TrailerCountErr()) {

		logger.logError(
			string("WARNING -- Trailer count error flag. Got ")
			+ to_string(e.Trailer().TDCTlrCount())
			+ " trailer(s)!",
			EVENT_WARN
		);

	}

	/*
	if(e.Trailer().HitCount() != e.Signals().size() + [error_word_count]) {

		logger.logError(
			string("WARNING -- Hit count in trailer = ")
			+ to_string(e.Trailer().HitCount())
			+ ", real hit count = "
			+ to_string(e.Signals().size())
			+ ", error word count = TODO: IMPLEMENT",
			EVENT_WARN
		);

	}
	*/

	int TDCHeaderCount  = 0;
	int TDCTrailerCount = 0;

	for(const Signal &sig : e.Signals()) {

		if(sig.isTDCHeader ()) ++TDCHeaderCount ;
		if(sig.isTDCTrailer()) ++TDCTrailerCount;

	}

	if(TDCHeaderCount != e.Trailer().TDCHdrCount()) {

		logger.logError(
			string("WARNING -- ")
			+ to_string(TDCHeaderCount)
			+ " found in data, event trailer indicates "
			+ to_string(e.Trailer().TDCHdrCount())
			+ "!",
			EVENT_WARN
		);

	}

	if(TDCTrailerCount != e.Trailer().TDCTlrCount()) {

		logger.logError(
			string("WARNING -- ")
			+ to_string(TDCTrailerCount)
			+ " found in data, event trailer indicates "
			+ to_string(e.Trailer().TDCTlrCount())
			+ "!",
			EVENT_WARN
		);

	}

	// TODO: Validate trailer hit count against e.signals().size() + error_word_count
	//         -- make sure to exclude TDC headers and trailers

}

void processEvent(Event &e) {

	DoHitFinding(&e, *TimeCorrection::getInstance(), 0, 0);
	// TODO: No hit clustering?

	e.SetPassCheck(true);
	e.CheckClusterTime();

}