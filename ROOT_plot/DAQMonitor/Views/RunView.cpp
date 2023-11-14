/**
 * @file RunView.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include "DAQMonitor/DataRunOperations.cpp"

#include "DAQMonitor/Views/Components/CanvasSelector.cpp"
#include "DAQMonitor/Views/Components/RunStats.cpp"

using namespace std;

class RunView : public TGMainFrame {

public:

	RunView(const TGWindow *p);

private:

	// VIEW

	RunStats       *runStats      ;
	CanvasSelector *canvasSelector;

	// CONNECTIONS

	void makeConnections();

};

void RunView::makeConnections() {

	// TODO: Before we reimplement this, we need to handle the case where 
	//       RunView stops the run while something else is stopping it or
	//       after something else has already stopped it.
	// Connect("CloseWindow()", "DataRun", nullptr, "stopRun()");

}

RunView::RunView(const TGWindow *p) : TGMainFrame(p) {

	runStats = new RunStats(this);
	AddFrame(runStats, new TGLayoutHints(kLHintsCenterX, 5, 5, 5, 5));

	canvasSelector = new CanvasSelector(this, "Select Plots");
	AddFrame(canvasSelector, new TGLayoutHints(kLHintsCenterX));

	makeConnections();

}