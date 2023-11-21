/**
 * @file HistogramPlotter.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <vector>
#include <string>

#include "macros/UIFramework/UISignals.cpp"

#include "src/DataModel/DAQData.cpp"

using namespace std;

class HistogramPlotter : public TGMainFrame {

public:

	HistogramPlotter(
		const TGWindow *p, 
		vector<TH1*> histograms, 
		const string &title = "", 
		int w = 1, 
		int h = 1,
		int rows = 1
	);

	~HistogramPlotter();

	void update();

private:

	// DATA

	vector<TH1*> histograms;

	// VIEW

	TRootEmbeddedCanvas *canvas;

	// CONNECTIONS

	void makeConnections ();
	void breakConnections();

};

HistogramPlotter::HistogramPlotter(
	const TGWindow *p, 
	vector<TH1*> histograms, 
	const string &title = "", 
	int w = 1, 
	int h = 1,
	int rows = 1
) : TGMainFrame(p), histograms(histograms) {

	canvas = new TRootEmbeddedCanvas(title.data(), this, w, h);
	AddFrame(canvas, new TGLayoutHints(kLHintsCenterX));

	int count = histograms.size();
	int cols = count / rows;
    if(count % rows != 0) ++cols;

    canvas->GetCanvas()->Divide(cols, rows);

	makeConnections();

	update();

}

HistogramPlotter::~HistogramPlotter() {

	breakConnections();

}

void HistogramPlotter::update() {

	for(int i = 0; i < histograms.size(); ++i) {

		canvas->GetCanvas()->cd(i + 1);

		// TODO: Try to break this dependence on DAQData. It's messy as heck,
		//       and makes it weird to use the HistogramPlotter for anything 
		//       that's not in DAQData

		DAQData &data = DAQData::getInstance();

		data.lock();
		histograms[i]->Draw();
		data.unlock();

	}

	canvas->GetCanvas()->Update();

}

void HistogramPlotter::makeConnections() {

	UISignalBus::getInstance().Connect(
		"onUpdate()", 
		"HistogramPlotter", 
		this, 
		"update()"
	);

}

void HistogramPlotter::breakConnections() {

	UISignalBus::getInstance().Disconnect("onUpdate()", this, "update()");

}