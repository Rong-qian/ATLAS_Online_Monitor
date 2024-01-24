/**
 * @file DAQData.cpp
 *
 * @brief Singleton storage for data retrieved from the DAQ. 
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 * 
 * CAUTION: DAQData is not innately threadsafe, and must be locked
 * before modification to ensure thread safety. Locking and unlocking
 * functions are provided.
 */

#pragma once

#include <vector>

#include "src/DataModel/Plots.cpp"
#include "src/Event.cpp"

using namespace std;
using namespace Muon;

struct DAQData {

    // A vector of fully processed and populated events, with hit finding
    // already performed.
    vector<Event> processedEvents;

    // A vector of the fully processed and populated events decoded in the
    // last decode iteration, with hit finding already performed.
    // Excludes empty events.
    vector<Event> newEvents;

    // A struct of histograms containing aggregate event data. See
    // src/DataModel/Plots.cpp.
    Plots plots;

    int totalEventCount = 0;
    int packetCount     = 0;

    void lock  () const; // Locks an internal mutex
    void unlock() const; // Unlocks an internal mutex

    /**
     * Clears all data stored in the DAQData. Clear is called at the
     * beginning of each run to clear any latent data from the last run.
     */
    void clear();

    // DAQData should not be copied
    DAQData       (      DAQData &other) = delete;
    void operator=(const DAQData &other) = delete;

    // Access the singleton DAQData instance by reference here.
    static DAQData &getInstance();

private:

    DAQData();
    mutable mutex dataLock;

};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

DAQData::DAQData() {}

DAQData &DAQData::getInstance() {

    static DAQData data;

    return data;

}

void DAQData::lock  () const { dataLock.lock  (); }
void DAQData::unlock() const { dataLock.unlock(); }

void DAQData::clear () {

    processedEvents.clear();
    plots.clear();
    totalEventCount = 0;
    packetCount     = 0;

}