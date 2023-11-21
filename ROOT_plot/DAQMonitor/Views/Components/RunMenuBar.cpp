/**
 * @file RunMenuBar.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <vector>
#include <string>

#include "macros/UIFramework/UIException.cpp"

#include "DAQMonitor/Plotting/Views/HistogramPlotter.cpp"

#include "DAQMonitor/PlotSaving/PlotSavingOperations.cpp"

#include "DAQMonitor/ErrorHandling/Views/ErrorView.cpp"

#include "src/DataModel/DAQData.cpp"

using namespace std;

// TODO: Pull TDC and ADC plot submenus to other components
// TODO: Break direct dependence on DAQData. Interface with the data model
//       through an operation.

const int SAVE_PLOTS_ID      = 1;
const int SHOW_ERROR_LOG_ID  = 2;
const int SHOW_NOISE_RATE_ID = 3;
const int CLOSE_ALL_ID       = 4;
const int PLOT_ENTRY_BASE    = 300000;

const int OVERVIEW_TDC       = -1;

struct PlotEntry {

public:

	PlotEntry() {

		// TODO: Do we care that this doesn't get reset when the run view
		//       is closed?
		static int instanceNum = PLOT_ENTRY_BASE;

		entryID = instanceNum;

		++instanceNum;

	}

	HistogramPlotter *plot = nullptr;
	int tdc = OVERVIEW_TDC;

	int getEntryID() const { return entryID; }

private:

	int entryID;

};

class RunMenuBar : public TGMenuBar {

public:

	RunMenuBar(const TGWindow *p);

	void handleFileMenuActivate(int id);
	void handleViewMenuActivate(int id);

	// METHODS

	void openADCOverview  ();
	void closeADCOverview ();

	void openTDCOverview  ();
	void closeTDCOverview ();

	void openErrorViewer  ();
	void closeErrorViewer ();

	void openNoiseDisplay ();
	void closeNoiseDisplay();

	void openADCChannelPlot (PlotEntry &entry);
	void closeADCChannelPlot(PlotEntry &entry);

	void openTDCChannelPlot (PlotEntry &entry);
	void closeTDCChannelPlot(PlotEntry &entry);

private:

	void closeAllWindows();

	// DATA

	DAQData &data;

	// COMPONENTS

	TGPopupMenu *fileMenu;

	TGPopupMenu *viewMenu;

		TGPopupMenu *adcPlotSelector;
		TGPopupMenu *tdcPlotSelector;

	// SUBVIEWS

	PlotEntry adc_overview;
	PlotEntry tdc_overview;

	vector<PlotEntry> adc_chnl_plots;
	vector<PlotEntry> tdc_chnl_plots;

	ErrorView *errorViewer;

	// CONNECTIONS

	void makeConnections();

};

RunMenuBar::RunMenuBar(const TGWindow *p) 
	: TGMenuBar(p, 100, 20, kHorizontalFrame), data(DAQData::getInstance()) {

	adc_chnl_plots.resize(data.plots.p_adc_time.size()          );
	tdc_chnl_plots.resize(data.plots.p_tdc_time_corrected.size());

	for(int i = 0; i < adc_chnl_plots.size(); ++i) {

		adc_chnl_plots[i].tdc = i;

	}

	for(int i = 0; i < tdc_chnl_plots.size(); ++i) {

		tdc_chnl_plots[i].tdc = i;

	}

	errorViewer  = nullptr;

	fileMenu = new TGPopupMenu(gClient->GetRoot());
	AddPopup("&File", fileMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft));

		fileMenu->AddEntry("Save Plots", SAVE_PLOTS_ID);


	viewMenu = new TGPopupMenu(gClient->GetRoot());
	AddPopup("&View", viewMenu, new TGLayoutHints(kLHintsTop | kLHintsLeft));

		viewMenu->AddSeparator();

		adcPlotSelector = new TGPopupMenu(gClient->GetRoot());
		viewMenu->AddPopup("ADC Plots", adcPlotSelector);

			adcPlotSelector->AddEntry("ADC Overview", adc_overview.getEntryID());

			adcPlotSelector->AddSeparator();

			for(const PlotEntry &entry : adc_chnl_plots) {

				adcPlotSelector->AddEntry(
					string("TDC ") + entry.tdc, entry.getEntryID()
				);

			}

		tdcPlotSelector = new TGPopupMenu(gClient->GetRoot());
		viewMenu->AddPopup("TDC Plots", tdcPlotSelector);

			tdcPlotSelector->AddEntry("TDC Overview", tdc_overview.getEntryID());

			tdcPlotSelector->AddSeparator();

			for(const PlotEntry &entry : tdc_chnl_plots) {

				tdcPlotSelector->AddEntry(
					string("TDC ") + entry.tdc, entry.getEntryID()
				);

			}

		viewMenu->AddSeparator();

		viewMenu->AddEntry("Noise Rate", SHOW_NOISE_RATE_ID);
		viewMenu->AddEntry("Error Log", SHOW_ERROR_LOG_ID);

		viewMenu->AddSeparator();

		viewMenu->AddEntry("Close All", CLOSE_ALL_ID);

	makeConnections();

}

void RunMenuBar::handleFileMenuActivate(int id) {

	switch(id) {

		case SAVE_PLOTS_ID:
			PlotSaving::savePlots();
			break;

	}

}

void RunMenuBar::handleViewMenuActivate(int id) {

	for(PlotEntry &entry : adc_chnl_plots) {

		if(entry.getEntryID() == id) {

			if(adcPlotSelector->IsEntryChecked(entry.getEntryID())) {

				entry.plot->CloseWindow();
				closeADCChannelPlot(entry);

			} else {

				openADCChannelPlot(entry);

			}

			return;

		}

	}

	for(PlotEntry &entry : tdc_chnl_plots) {

		if(entry.getEntryID() == id) {

			if(tdcPlotSelector->IsEntryChecked(entry.getEntryID())) {

				entry.plot->CloseWindow();
				closeTDCChannelPlot(entry);

			} else {

				openTDCChannelPlot(entry);

			}

			return;

		}

	}

	if(id == adc_overview.getEntryID()) {

		if(adcPlotSelector->IsEntryChecked(adc_overview.getEntryID())) {

			adc_overview.plot->CloseWindow();
			closeADCOverview();

		} else {

			openADCOverview();

		}

		return;

	}

	if(id == tdc_overview.getEntryID()) {

		if(tdcPlotSelector->IsEntryChecked(tdc_overview.getEntryID())) {

			tdc_overview.plot->CloseWindow();
			closeTDCOverview();

		} else {

			openTDCOverview();

		}

		return;

	}

	if(id == SHOW_NOISE_RATE_ID) {

		if(viewMenu->IsEntryChecked(SHOW_NOISE_RATE_ID)) {

			closeNoiseDisplay();

		} else {

			openNoiseDisplay();

		}

		return;

	}

	if(id == SHOW_ERROR_LOG_ID) {

		if(viewMenu->IsEntryChecked(SHOW_ERROR_LOG_ID)) {

			errorViewer->CloseWindow();
			closeErrorViewer();

		} else {

			openErrorViewer();

		}

		return;

	}

	if(id == CLOSE_ALL_ID) {

		closeAllWindows();

		return;

	}

	throw UIException("Menu option not recognized.");

}

void RunMenuBar::makeConnections() {

	fileMenu->Connect("Activated(Int_t)", "RunMenuBar", this, "handleFileMenuActivate(int)");
	viewMenu->Connect("Activated(Int_t)", "RunMenuBar", this, "handleViewMenuActivate(int)");

}

void RunMenuBar::openADCOverview() {

	vector<TH1*> plotList;
	plotList.reserve(data.plots.p_tdc_adc_time.size());

	for(TH1F *plot : data.plots.p_tdc_adc_time) {

		plotList.push_back(plot);

	}

	adc_overview.plot = new HistogramPlotter(
		gClient->GetRoot(), 
		plotList, 
		"ADC Plots",
		1800,
		900,
		3
	);

	adc_overview.plot->SetWindowName("ADC Plots");
	adc_overview.plot->MapSubwindows();
	adc_overview.plot->Resize(adc_overview.plot->GetDefaultSize());
	adc_overview.plot->MapWindow();

	adc_overview.plot->Connect("CloseWindow()", "RunMenuBar", this, "closeADCOverview()");

	adcPlotSelector->CheckEntry(adc_overview.getEntryID());

}

void RunMenuBar::closeADCOverview() {

	adcPlotSelector->UnCheckEntry(adc_overview.getEntryID());

	adc_overview.plot = nullptr;

}

void RunMenuBar::openTDCOverview() {

	vector<TH1*> plotList;
	plotList.reserve(data.plots.p_tdc_tdc_time_corrected.size());

	for(TH1F *plot : data.plots.p_tdc_tdc_time_corrected) {

		plotList.push_back(plot);

	}

	tdc_overview.plot = new HistogramPlotter(
		gClient->GetRoot(), 
		plotList, 
		"TDC Plots",
		1800,
		900,
		3
	);

	tdc_overview.plot->SetWindowName("TDC Plots");
	tdc_overview.plot->MapSubwindows();
	tdc_overview.plot->Resize(tdc_overview.plot->GetDefaultSize());
	tdc_overview.plot->MapWindow();

	tdc_overview.plot->Connect("CloseWindow()", "RunMenuBar", this, "closeTDCOverview()");

	tdcPlotSelector->CheckEntry(tdc_overview.getEntryID());

}

void RunMenuBar::closeTDCOverview() {

	tdcPlotSelector->UnCheckEntry(tdc_overview.getEntryID());

	tdc_overview.plot = nullptr;

}

void RunMenuBar::openADCChannelPlot(PlotEntry &entry) {

	vector<TH1*> plotList;
	plotList.reserve(data.plots.p_adc_time[entry.tdc].size());

	for(TH1F *plot : data.plots.p_adc_time[entry.tdc]) {

		plotList.push_back(plot);

	}

	string plotTitle = "ADC Channel Plots, TDC ";
	plotTitle += to_string(entry.tdc);

	entry.plot = new HistogramPlotter(
		gClient->GetRoot(), 
		plotList, 
		plotTitle.data(),
		1800,
		900,
		3
	);

	entry.plot->SetWindowName(plotTitle.data());
	entry.plot->MapSubwindows();
	entry.plot->Resize(
		entry.plot->GetDefaultSize()
	);
	entry.plot->MapWindow();

	adcPlotSelector->CheckEntry(entry.getEntryID());

}

void RunMenuBar::closeADCChannelPlot(PlotEntry &entry) {

	adcPlotSelector->UnCheckEntry(entry.getEntryID());

	entry.plot = nullptr;

}

void RunMenuBar::openTDCChannelPlot(PlotEntry &entry) {

	vector<TH1*> plotList;
	plotList.reserve(data.plots.p_tdc_time_corrected[entry.tdc].size());

	for(TH1F *plot : data.plots.p_tdc_time_corrected[entry.tdc]) {

		plotList.push_back(plot);

	}

	string plotTitle = "TDC Channel Plots, TDC ";
	plotTitle += to_string(entry.tdc);

	entry.plot = new HistogramPlotter(
		gClient->GetRoot(), 
		plotList, 
		plotTitle.data(),
		1800,
		900,
		3
	);

	entry.plot->SetWindowName(plotTitle.data());
	entry.plot->MapSubwindows();
	entry.plot->Resize(entry.plot->GetDefaultSize());
	entry.plot->MapWindow();

	tdcPlotSelector->CheckEntry(entry.getEntryID());


}
	
void RunMenuBar::closeTDCChannelPlot(PlotEntry &entry) {

	tdcPlotSelector->UnCheckEntry(entry.getEntryID());

	entry.plot = nullptr;

}

void RunMenuBar::openErrorViewer() {

    errorViewer = new ErrorView(gClient->GetRoot());

    errorViewer->SetWindowName("Error Log");
    errorViewer->MapSubwindows();
    errorViewer->Resize(errorViewer->GetDefaultSize());
    errorViewer->MapWindow();

    errorViewer->Connect("CloseWindow()", "RunMenuBar", this, "closeErrorViewer()");

    viewMenu->CheckEntry(SHOW_ERROR_LOG_ID);

}

void RunMenuBar::closeErrorViewer() {

	viewMenu->UnCheckEntry(SHOW_ERROR_LOG_ID);

	errorViewer = nullptr;

}

void RunMenuBar::openNoiseDisplay() {

	cout << "Noise display not yet implemented" << endl;
	viewMenu->CheckEntry(SHOW_NOISE_RATE_ID);

}

void RunMenuBar::closeNoiseDisplay() {

	viewMenu->UnCheckEntry(SHOW_NOISE_RATE_ID);
	cout << "Noise display not yet implemented" << endl;

}

void RunMenuBar::closeAllWindows() {

	if(errorViewer) {
		errorViewer->CloseWindow();
		closeErrorViewer();
	}

	closeNoiseDisplay();
	
	if(adc_overview.plot) {
		adc_overview.plot->CloseWindow();
		closeADCOverview();
	}

	if(tdc_overview.plot) {
		tdc_overview.plot->CloseWindow();
		closeTDCOverview();
	}

	for(PlotEntry &entry : adc_chnl_plots) {

		if(entry.plot) {

			entry.plot->CloseWindow();
			closeADCChannelPlot(entry);

		}

	}

	for(PlotEntry &entry : tdc_chnl_plots) {

		if(entry.plot) {

			entry.plot->CloseWindow();
			closeTDCChannelPlot(entry);

		}

	}

}