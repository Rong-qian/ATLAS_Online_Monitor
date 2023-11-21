/**
 * @file RunView.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include "DAQMonitor/Views/Components/RunStats.cpp"
#include "DAQMonitor/Views/Components/RunMenuBar.cpp"

using namespace std;

class RunView : public TGMainFrame {

public:

	RunView(const TGWindow *p);

private:

	// VIEW

	RunMenuBar         *menuBar       ;
	RunStats           *runStats      ;

	// CONNECTIONS

	void makeConnections();

};

void RunView::makeConnections() {

	// TODO: Before we reimplement this, we need to handle the case where 
	//       RunView stops the run while something else is stopping it or
	//       after something else has already stopped it.
	// Connect("CloseWindow()", "DataRun", nullptr, "stopRun()");
	Connect("CloseWindow()", "RunStats", runStats, "teardown()");

}

RunView::RunView(const TGWindow *p) 
	: TGMainFrame(p, 1, 1, kFixedSize) {

	menuBar = new RunMenuBar(this);
	AddFrame(menuBar, new TGLayoutHints(kLHintsTop | kLHintsLeft));

	runStats = new RunStats(this);
	AddFrame(
		runStats, new TGLayoutHints(
			kLHintsExpandX | kLHintsExpandY
		)
	);

	makeConnections();

	Resize(300, 250);

}