/**
 * @file StartDAQ.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 * 
 * NOTE: StartDAQ blocks SIGINT. Termination can be forced by
 *       typing CTRL+C twice.
 */

// TODO: Consider an operator pattern where the operations file maintains a
//       global pointer to the associated view and implements an initializer

#pragma once

#include <thread>
#include <mutex>
#include <chrono>
#include <cmath>

#include "monitorConfig.cpp"

#include "macros/DAQState.cpp"
#include "macros/ErrorLogger.cpp"
#include "macros/UIFramework/UIException.cpp"
#include "macros/UIFramework/UISignals.cpp"
#include "macros/UIWindows/AlertBox/AlertOperations.cpp"

#include "MainMenu/Views/EntryView.cpp"

#include "src/ProgramControl/Terminator.cpp"
#include "src/ProgramControl/SigHandlers.cpp"
#include "src/ProgramControl/Threads.cpp"

using namespace std;
using namespace Muon;


// COPIED TODOS FROM START_MONITOR:
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

// TODO: Look at this for using QT in root:
//       https://root.cern/doc/v606/QtFileDialog_8C.html

const string STATE_STORAGE = "settings.txt";

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void startDAQ();

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void StartDAQ() {

    // TODO: We want a console logger class to replace cout too.
    ErrorLogger::getInstance().setOutputStream(cerr);

    // Read in DAQState from file
    State::DAQState state = State::DAQState::getState();
    state.load(STATE_STORAGE);
    state.commit();

    // THIS MUST BE CALLED BEFORE STARTING ANY THREADS.
    // It intercepts SIGINT/SIGTERM/SIGQUIT to cleanly terminate threads.
    setTerminationHandlers(flagForTermination);

    ///////////////////////////////////////////////////////////////////////////
    //////////////////////////// SET UP UI ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // TODO: Move this to EntryOperations
    ProgramFlow::threadLock.lock();
    ProgramFlow::threads.emplace_back(

        thread([](){

            // Create the main window
            TGMainFrame *mainFrame = new TGMainFrame(gClient->GetRoot());

            // Create the menu displayed on the main window
            EntryView *menu = new EntryView(mainFrame, 500, 350);
            mainFrame->AddFrame(menu, new TGLayoutHints(kLHintsCenterX));

            // Set up the main window now that it has all its components
            mainFrame->SetWindowName("DAQ Manager");
            mainFrame->MapSubwindows();
            mainFrame->Resize(mainFrame->GetDefaultSize());
            mainFrame->MapWindow();

            // Wire up the main window's close button to the terminator
            mainFrame->Connect(
                "CloseWindow()", 
                "EntryOperations", 
                nullptr, 
                "exitAll()"
            );

            // Run the UI update loop
            // TODO: Check Resource Red to remember how to properly track
            //       dropped frames
            // TODO: Now that the logic is a little more complicated, pull it
            //       out
            // TODO: Include a setting regarding whether to warn for slow
            //       frames
            int frameNum      = 0;
            double slowFrames = 0;
            while(!Terminator::getInstance().isTerminated()) {

                // Record the loop start time
                auto UIUpdateStartTime = chrono::high_resolution_clock::now();

                try {

                    gSystem->ProcessEvents();

                    if(frameNum % UI_UPDATE_FRAMES == 0) {

                        // Reset it every time to prevent overflows
                        frameNum = 0;

                        UISignalBus::getInstance().onUpdate();

                    }

                    ++frameNum;

                } catch (UIException &e) {

                    // TODO: More sophisticated error handling. We should
                    //       be logging these errors too.
                    // TODO: Multiple types of error titles
                    Error::popupAlert(
                        gClient->GetRoot(), "Error", e.what()
                    );

                }

                // Record the loop end time
                auto UIUpdateEndTime = chrono::high_resolution_clock::now();

                // Compute how long it took the loop to run
                chrono::duration<double> loopTime = UIUpdateEndTime - UIUpdateStartTime;
                chrono::milliseconds updateTime = chrono::duration_cast<chrono::milliseconds>(loopTime);

                DAQState state = DAQState::getState();

                bool warnSlowFrames = state.persistentState.warnSlowFrames;

                // Compute how long the thread should sleep to maintain a 
                // consistent frame rate
                int sleepTime = (int)(1000. / GUI_REFRESH_RATE) - (int)(updateTime.count());

                // Deal with the case where the update loop took too long
                if(sleepTime < 0) {

                    if(warnSlowFrames) ++slowFrames;
                    sleepTime = 0;

                }

                if(warnSlowFrames) {

                    // If we are slow more than 5 times before slowFrames decays
                    // we should alert the user.
                    if(slowFrames > state.persistentState.slowFrameTolerance) {

                        cerr << "UI WARNING: GUI refresh rate too fast for data updates. "
                             << "Consider reducing GUI update framerate." << endl;

                    }

                    // Linear decay rate for slowFrames
                    slowFrames -= state.persistentState.slowFrameDecayRate;
                    slowFrames = max(slowFrames, 0.);

                }

                this_thread::sleep_for(
                    chrono::milliseconds(sleepTime)
                );

            }

            // TODO: Make sure all windows are closed without double deletion.
            //       Otherwise we will find that the GUI will stay open until
            //       every thread terminates.
            //         -- Really though we want the data capture thread to
            //            terminate immediately
            //         -- A window manager might be useful here. ExitAll
            //            can delete all the windows, and then we can just
            //            rewire the normal CloseWindow signal

        })

    );
    ProgramFlow::threadLock.unlock();

    ///////////////////////////////////////////////////////////////////////////
    ////////////////////////////// CLEANUP ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // Wait for the threads to be done before terminating
    ProgramFlow::joinAllThreads();
    
    // Save DAQState to file
    //   -- TODO: Consider saving DAQState on startRun too
    state.update();
    state.save(STATE_STORAGE);

    cout << "Shutdown" << endl;

    // This makes sure the GUI closes when the main thread terminates
    gApplication->Terminate(0);

}