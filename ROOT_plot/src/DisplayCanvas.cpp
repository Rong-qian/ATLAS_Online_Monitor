// ROOT includes
#include "TDirectory.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TString.h"

//Muon Reconstruction includes
#include "src/Geometry.cpp"
#include "src/TrackParam.cpp"
#include "src/RTParam.cpp"
#include "src/TimeCorrection.cpp"
#include "src/HitFinder.cpp"
#include "src/Event.cpp"
#include "src/EventDisplay.cpp"
// New Reco Includes
#include "src/T0Fit.h"
#include "src/Observable.cpp"
#include "src/Parameterization.cpp"
#include "src/ResolutionResult.cpp"
#include "src/Track.cpp"
#include "src/IOUtility.cpp"
#include "src/Optimizer.cpp"
#include "src/DataModel/DAQData.cpp"

//#define SAVE_TRACKS_OUT_OF_ROOT // comment this line if you don't need to save plots out of rootfile
//#define ERROR_WORD 29
//#define HEADER_WORD 31
//#define TRAILER_WORD 30
//#define SPEEDFACTOR 1
//#define NEWTDC_NUMBER 17
//#define WIDTH_RES 1

#ifndef OnlineCANVAS
#define OnlineCANVAS
namespace Muon {

	class DisplayCanvas{
		public :
			DisplayCanvas(const DAQData &data);
			void UpdateCanvas(const DAQData &data);
			int displayed_event_counter;
			mutable mutex displayLock;
			void lock  () const; // Locks an internal mutex
			void unlock() const; // Unlocks an internal mutex
		private:
			TFile *p_output_rootfile;
			TCanvas *rate_canvas, *trigger_rate_canvas, *residual_canvas, *EDCanvas, *eff_canvas;
			RTParam* rtp;
			EventDisplay *ed;
			Event ed_event;
			TrackParam tp;
			TimeCorrection tc;
			unsigned long total_triggers, total_events, total_triggers_pass, total_events_pass;
			unsigned long total_signals, total_signals_pass, total_events_fail, event_print, event_signals;
			TString h_name;
			int pad_num;
			double match_window = 1.5;
	};

	DisplayCanvas::DisplayCanvas(const DAQData &data){
		residual_canvas = new TCanvas("c1", "Residuals", 2100, 900, 400, 300);
		EDCanvas = new TCanvas("c2", "Event Display", 2700, 900, 1200, 800);
		eff_canvas = new TCanvas("c3", "Efficiency", 2300, 900, 400, 300);

		ed = new EventDisplay();
		ed->SetCanvas(EDCanvas);
  		gSystem->ProcessEvents();

		residual_canvas->cd();
		data.plots.residuals->Draw();
		residual_canvas->Modified();
		residual_canvas->Update();
		eff_canvas->cd();
		data.plots.tube_efficiency->Draw("colz");
		eff_canvas->Modified();
		eff_canvas->Update();
		cout<<"Canvases created and divided.\n"<<endl;

	}


	void DisplayCanvas::UpdateCanvas(const DAQData &data)
	{
		gStyle->SetOptStat(1110); // std, mean, entris, name printed
		// gStyle->SetTitleX(999.);//hist no title
		// gStyle->SetTitleY(999.);
		gStyle->SetStatY(0.9);
		// Set y-position (fraction of pad size)
		gStyle->SetStatX(0.9);
		// Set x-position (fraction of pad size)
		//gStyle->SetStatW(0.25);
		// Set width of stat-box (fraction of pad size)
		//->SetStatH(0.25);
		// Set height of stat-box (fraction of pad size)

		// Update plots
		residual_canvas->cd();
		data.plots.residuals->Draw();
		residual_canvas->Update();
		eff_canvas->Modified();
		eff_canvas->Update();
		for (Event ed_event:*data.plots.ed_events_ptr)
		{
			EDCanvas->cd();
			ed->DrawEvent(ed_event, Geometry::getInstance(), NULL);
			gSystem->ProcessEvents();
			ed->Clear();
			gSystem->Sleep(4000);
		}
	} // End update canvas
}
void DisplayCanvas::lock  () const { displayLock.lock  (); }
void DisplayCanvas::unlock  () const { displayLock.unlock  (); }

#endif
