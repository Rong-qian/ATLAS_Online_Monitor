/**
 * @file DAQState.cpp
 *
 * @brief Global state for the DAQ interface.
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 * 
 * To add new temporary state, simply add a member to the TemporaryState
 * struct. To add new persistent state, add a member to the PersistantState 
 * struct, then update DAQState::readPersistentState() and 
 * DAQState::writePersistentState with readers and writers for the new member.
 */

#pragma once

#include <string>
#include <mutex>
#include <iostream>
#include <algorithm>
#include <map>

using namespace std;

namespace Muon {
namespace State {

	string trim(const string &str);

	// DAQState::dataSource values
	const int DAT_FILE_SOURCE       = 1;
	const int NETWORK_DEVICE_SOURCE = 2;

	struct PersistentState {

		int    dataSource         = NETWORK_DEVICE_SOURCE;

		bool   warnSlowFrames     = true                 ;
		int    slowFrameTolerance = 5                    ;
		double slowFrameDecayRate = 0.2                  ;

		string inputFilename      = ""                   ;
		string inputDevicename    = ""                   ;

	};

	struct TemporaryState {

		bool   runStarted = false;
		string runLabel   = ""   ;

	};

	class DAQState {

	public:

		PersistentState persistentState;
		TemporaryState  tempState      ;

		void update(); // Sets this state to match the current committed state
		bool commit(bool force = false);

		bool load(const string &filename);
		bool save(const string &filename);

		bool readPersistentState (istream &in );
		void writePersistentState(ostream &out);

		static DAQState getState();

	private:

		unsigned int commitNumber = 0;

		static DAQState masterState;

		static mutex stateLock;

		DAQState();

		bool isOutdated();

	};

	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	bool DAQState::load(const string &filename) {

		ifstream stateSource(filename);
		if(!stateSource.is_open()) {

			return false;

		}

		if(!readPersistentState(stateSource)) {

			stateSource.close();
			return false;

		}

		stateSource.close();
		return true;

	}

	bool DAQState::save(const string &filename) {

		ofstream stateDest(filename);
		if(!stateDest.is_open()) {

			return false;

		}

		writePersistentState(stateDest);

		stateDest.close();

		return true;

	}


	void DAQState::writePersistentState(ostream &out) {

		out << "Data Source: "           << persistentState.dataSource         << endl;
		out << "Input Device Name: "     << persistentState.inputDevicename    << endl;
		out << "Input File Name: "       << persistentState.inputFilename      << endl;
		out << "Warn Slow Frames: "      << persistentState.warnSlowFrames     << endl;
		out << "Slow Frame Tolerance: "  << persistentState.slowFrameTolerance << endl;
		out << "Slow Frame Decay Rate: " << persistentState.slowFrameDecayRate << endl;

	}

	bool DAQState::readPersistentState(istream &in) {

		map<string, string> tokens;

		string key;
		while(getline(in, key, ':')) {

			string value;
			getline(in, value);

			key   = trim(key  );
			value = trim(value);

			tokens[key] = value;

		}

		if(tokens["Data Source"] == "") {

			persistentState.dataSource = 0;

		} else {

			persistentState.dataSource = stoi(tokens["Data Source"]);

		}

		if(tokens["Warn Slow Frames"] == "") {

			persistentState.warnSlowFrames = true;

		} else {

			persistentState.warnSlowFrames = stoi(tokens["Warn Slow Frames"]);

		}

		if(tokens["Slow Frame Tolerance"] == "") {

			persistentState.slowFrameTolerance = 5;

		} else {

			persistentState.slowFrameTolerance = stoi(tokens["Slow Frame Tolerance"]);

		}

		if(tokens["Slow Frame Decay Rate"] == "") {

			persistentState.slowFrameDecayRate = 0.2;

		} else {

			persistentState.slowFrameDecayRate = stod(tokens["Slow Frame Decay Rate"]);

		}

		persistentState.inputDevicename = tokens["Input Device Name"];
		persistentState.inputFilename   = tokens["Input File Name"  ];

		return true;

	}

	DAQState DAQState::masterState;

	mutex DAQState::stateLock;

	DAQState::DAQState() {}

	bool DAQState::commit(bool force = false) {

		stateLock.lock();

		if(!force) {

			if(isOutdated()) {

				// TODO: Better to throw an exception here?

				stateLock.unlock();
				return false;

			}

		}

		++commitNumber;
		masterState = *this;

		stateLock.unlock();

		return true;

	}

	void DAQState::update() {

		stateLock.lock();
		*this = masterState;
		stateLock.unlock();

	}

	DAQState DAQState::getState() {

		stateLock.lock();
		DAQState state = masterState;
		stateLock.unlock();

		return state;

	}

	bool DAQState::isOutdated() {

		return commitNumber < masterState.commitNumber;

	}

	string trim(const string &str) {

		string result = str;

		// Trim left whitespace
		result.erase(result.begin(), find_if(result.begin(), result.end(), [](unsigned char c) {

			return !isspace(c);

		}));

		// Trim right whitespace
		result.erase(find_if(result.rbegin(), result.rend(), [](unsigned char c) {

			return !isspace(c);

		}).base(), result.end());

		return result;

	}

}
}