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
#include "src/HitCluster.cpp"
#include "src/Event.cpp"
#include "src/EventDisplay.cpp"

// New Reco Includes
#include "src/T0Fit.h"
#include "src/TrackParam.cpp"
#include "src/Observable.cpp"
#include "src/Parameterization.cpp"
#include "src/RTParam.cpp"
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
			DisplayCanvas();
			void UpdateCanvas(const DAQData &data);

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

			std::vector<bool> isDrawn;
			std::vector<std::vector<double>> p_tdc_hit_rate;
			std::vector<double> p_tdc_hit_rate_x;
			std::vector<TH1F *> p_tdc_hit_rate_graph;
			std::vector<std::vector<Bool_t>> first_signal_flag;
			std::vector<std::vector<double>> nHits;
			std::vector<std::vector<double>> nMiss;

			TH1D *residuals;
			TH2D *tube_efficiency;

	};

	DisplayCanvas::DisplayCanvas(){
		total_triggers = 0;
		total_events_pass = 0;
		rtp = new RTParam(Geometry::getInstance());
		tp = TrackParam();
		tp.SetRT(rtp);
		tp.SetGEO(&Geometry::getInstance());
		tp.setVerbose(0);
  		tp.setMaxResidual(1000000);
		TimeCorrection tc = TimeCorrection();
		residuals = new TH1D("residuals", "Residuals;Residual[#mum];Number of hits/2#mum", 500, -500, 500);



		p_tdc_hit_rate_x.resize(Geometry::MAX_TDC_CHANNEL);
		p_tdc_hit_rate.resize(Geometry::MAX_TDC);
		p_tdc_hit_rate_graph.resize(Geometry::MAX_TDC);
		first_signal_flag.resize(Geometry::MAX_TDC);

		isDrawn.resize(Geometry::MAX_TDC);
		nHits.resize(Geometry::MAX_TUBE_LAYER);
		nMiss.resize(Geometry::MAX_TUBE_LAYER);
		for (size_t i = 0; i < Geometry::MAX_TUBE_LAYER; ++i)
		{
			nHits[i].resize(Geometry::MAX_TUBE_COLUMN);
			nMiss[i].resize(Geometry::MAX_TUBE_COLUMN);
		}

		for (size_t i = 0; i < Geometry::MAX_TDC; ++i)
		{
			p_tdc_hit_rate[i].resize(Geometry::MAX_TDC_CHANNEL);
			first_signal_flag[i].resize(Geometry::MAX_TDC_CHANNEL);
		}


		for (size_t i = 0; i < Geometry::MAX_TDC; ++i)
		{
			for (size_t j = 0; j < Geometry::MAX_TDC_CHANNEL; ++j)
			{
				p_tdc_hit_rate[i][j] = 0;
			}
		}

		for (int i = 0; i < Geometry::MAX_TDC_CHANNEL; i++)
		{
			p_tdc_hit_rate_x[i] = i;
		}

		for (int i = 0; i < Geometry::MAX_TDC_CHANNEL; i++)
		{
			p_tdc_hit_rate_x[i] = i;
		}
		for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
		{
			int local_tdc_id = tdc_id % 24;
			h_name.Form("tdc_%d_hit_rate", tdc_id);
			p_tdc_hit_rate_graph[tdc_id] = new TH1F(h_name, TString::Format("tdc_%d_tdc_hit_rate;Channel No.;Rate(Hz)", local_tdc_id), 24, -0.5, 23.5);
			if (tdc_id < 24)
				p_tdc_hit_rate_graph[tdc_id]->SetFillColor(4);
			else
				p_tdc_hit_rate_graph[tdc_id]->SetFillColor(kRed);

			p_tdc_hit_rate_graph[tdc_id]->SetBarWidth(0.4);
			p_tdc_hit_rate_graph[tdc_id]->SetBarOffset(0.1 + 0.4 * (tdc_id < 24));
			p_tdc_hit_rate_graph[tdc_id]->SetStats(0);

			p_tdc_hit_rate_graph[tdc_id]->SetMaximum(1);
			p_tdc_hit_rate_graph[tdc_id]->SetMinimum(0);
		}
		tube_efficiency = new TH2D("tube_efficiency", ";Layer;Column",
				Geometry::MAX_TUBE_COLUMN, -0.5, Geometry::MAX_TUBE_COLUMN - 0.5,
				Geometry::MAX_TUBE_LAYER, -0.5, Geometry::MAX_TUBE_LAYER - 0.5);

		rate_canvas = new TCanvas("c3", "Hit Rate Plots", 2160, 0, 1800, 800);
		rate_canvas->Divide(6, 4);
		trigger_rate_canvas = new TCanvas("c4", "Trigger Board", 1440, 750, 400, 300);
		residual_canvas = new TCanvas("c5", "Residuals", 2100, 900, 400, 300);
		EDCanvas = new TCanvas("c6", "Event Display", 2700, 900, 800, 800);
		eff_canvas = new TCanvas("C7", "Efficiency", 2300, 900, 400, 300);

		ed = new EventDisplay();
		ed->SetCanvas(EDCanvas);


		for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
		{
			if (Geometry::getInstance().IsActiveTDC(tdc_id) || tdc_id == Geometry::getInstance().TRIGGER_MEZZ)
			{
				int local_tdc_id = tdc_id % 24;
				if (tdc_id == Geometry::getInstance().TRIGGER_MEZZ)
				{
					trigger_rate_canvas->cd();
        			p_tdc_hit_rate_graph[tdc_id]->Draw("B");
				}	
				else
				{	
					pad_num = Geometry::getInstance().GetPad(local_tdc_id);
					TString opts;
					if (isDrawn[tdc_id])
						opts = "same";
					else
						opts = "";
					rate_canvas->cd(pad_num+1);
					// gPad->SetLogy();
					p_tdc_hit_rate_graph[tdc_id]->Draw(opts + " B");
					
					// if (gSystem->ProcessEvents()) break;
					isDrawn[tdc_id] = 1;
				}
			}
		}

		rate_canvas->cd();
		rate_canvas->Modified();
		rate_canvas->Update();
		//gSystem->ProcessEvents();
		residual_canvas->cd();
		residuals->Draw();
		residual_canvas->Modified();
		residual_canvas->Update();
		eff_canvas->cd();
		tube_efficiency->Draw("colz");
		eff_canvas->Modified();
		eff_canvas->Update();
		cout<<"Canvases created and divided.\n"<<endl;
		p_output_rootfile = new TFile("output.root", "RECREATE");

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
		cout<<"Processing "<< data.totalEventCount <<" events \n"<<endl;

		for (Int_t tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
		{
			if (Geometry::getInstance().IsActiveTDC(tdc_id))
			{
				//TODO: if p_tdc_time_original can be used
				for (int tdc_chnl_id = 0; tdc_chnl_id != Geometry::MAX_TDC_CHANNEL; tdc_chnl_id++)
				{	
					p_tdc_hit_rate[tdc_id][tdc_chnl_id] =
						data.plots.p_tdc_time_corrected[tdc_id][tdc_chnl_id]->GetEntries() / match_window * 1000 / data.processedEvents.size();
				}
			}
		}
		float max = 0.0;

		for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; ++tdc_id)
		{
			isDrawn[tdc_id] = 0;
			max = 0.0;
			for (int iBin = 1; iBin <= p_tdc_hit_rate_graph[tdc_id]->GetNbinsX(); ++iBin)
			{
				//cout<<tdc_id<<" TestiBin:"<<iBin<<":" <<p_tdc_hit_rate[tdc_id][iBin - 1]<<"\n" ;
				p_tdc_hit_rate_graph[tdc_id]->SetBinContent(iBin, p_tdc_hit_rate[tdc_id][iBin - 1]);
				if (p_tdc_hit_rate[tdc_id][iBin - 1] > max)
					max = p_tdc_hit_rate[tdc_id][iBin - 1];
				p_tdc_hit_rate_graph[tdc_id]->SetMaximum(max);
			}
		}

		for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
		{
			string text_content;
			int CSM = tdc_id >= 24;
			int local_tdc_id = tdc_id % 24;

			if (Geometry::getInstance().IsActiveTDC(tdc_id) || tdc_id == Geometry::getInstance().TRIGGER_MEZZ)
			{

				if (tdc_id == Geometry::getInstance().TRIGGER_MEZZ)
				{
					trigger_rate_canvas->cd();
					text_content = "Entries = " + to_string((int)total_triggers);
				}
				else
				{
					pad_num = Geometry::getInstance().GetPad(local_tdc_id);
					rate_canvas->cd(pad_num+1);					
					text_content = "Entries = " + to_string((int)data.plots.p_tdc_adc_time[tdc_id]->GetEntries());
				}
				TString h_name;

				h_name.Form("tdc_%d_hit_rate", local_tdc_id);
				int brother = 0;
				// if (CSM) brother = local_tdc_id;
				brother = local_tdc_id;
				// else brother = tdc_id + 18;

				max = (p_tdc_hit_rate_graph[tdc_id]->GetMaximum() > p_tdc_hit_rate_graph[brother]->GetMaximum()) ? p_tdc_hit_rate_graph[tdc_id]->GetMaximum() : p_tdc_hit_rate_graph[brother]->GetMaximum();

				max = (1 > max) ? 1 : max;

				p_tdc_hit_rate_graph[tdc_id]->SetMaximum(1.25 * max);

				TString opts;
				if (CSM && isDrawn[local_tdc_id])
					opts = "B same";
				else
					opts = "B";

				p_tdc_hit_rate_graph[tdc_id]->Draw(opts);
				TText *xlabel = new TText();
				xlabel->SetTextColor(kBlack);
				xlabel->SetNDC();
				xlabel->SetTextFont(42);
				xlabel->SetTextSize(0.05);
				xlabel->SetTextAngle(0);
				xlabel->DrawText(0.1 + float(CSM) * 0.4, 0.85, text_content.c_str());
				text_content = "Max  = " + to_string(TMath::MaxElement(MAX_CHNL_COUNT, &p_tdc_hit_rate[tdc_id][0])).substr(0, 6) + " kHz";
				xlabel->DrawText(0.1 + float(CSM) * 0.4, 0.8, text_content.c_str());
				text_content = "Mean = " + to_string(TMath::Mean(MAX_CHNL_COUNT, &p_tdc_hit_rate[tdc_id][0])).substr(0, 6) + " kHz";
				xlabel->DrawText(0.1 + float(CSM) * 0.4, 0.75, text_content.c_str());
				xlabel->SetTextColor(kRed);
				xlabel->DrawText(0.5, 0.9, "CSM 2");
				xlabel->SetTextColor(kBlue);
				xlabel->DrawText(0.1, 0.9, "CSM 1");
				xlabel->SetTextColor(kBlack);
				TLine *l = new TLine(-0.5, 0.5, 23.5, 0.5);
				l->Draw();
          	isDrawn[local_tdc_id] = 1;


			}
		}




		for (int i = 1; i != pad_num + 1; i++)
		{
			rate_canvas->cd(i);
			gPad->Modified();
		}

		// Update plots

		rate_canvas->cd();
		rate_canvas->Update();

		// Event display
		for (Event event : data.processedEvents){
			DoHitFinding(&event, tc, 0, 0);
			DoHitClustering(&event, Geometry::getInstance());
			ed_event = event;
			bool pass_event_check = kTRUE;
			event.SetPassCheck(pass_event_check);
			event.CheckClusterTime();
			if (pass_event_check){
				TTree *optTree = new TTree("optTree", "optTree");

				optTree->Branch("event", "Event", &event);
				optTree->Fill();
				tp.setTarget(optTree);
				tp.setRangeSingle(0);
				tp.setIgnoreNone();
				total_events_pass++;
				tp.optimize();
				for (Cluster c : event.Clusters())
				{
					for (Hit h : c.Hits())
					{
						residuals->Fill(tp.Residual(h) * 1000.0);
					}
				}
				// fill efficiency distribution
				double _hitX, _hitY;
				for (int tdc_index = 0; tdc_index < Geometry::MAX_TDC; tdc_index++) {
					for (int ch_index = 0; ch_index < Geometry::MAX_TDC_CHANNEL; ch_index++) {
						if (Geometry::getInstance().IsActiveTDCChannel(tdc_index, ch_index)){
							int iL,iC;
							Geometry::getInstance().GetHitLayerColumn(tdc_index, ch_index, &iL, &iC);
							Geometry::getInstance().GetHitXY(tdc_index, ch_index, &_hitX, &_hitY);
							// get track x position and figure out what tube(s) it may go through
							double trackDist = tp.Distance_XY(_hitX,_hitY);
							if (trackDist <= Geometry::column_distance / 2)
							{
								Bool_t tubeIsHit = kFALSE;

								for (Hit h : event.Hits())
								{

									int hit_layer;
									int hit_column;
									Geometry::getInstance().GetHitLayerColumn(h.TDC(), h.Channel(), &hit_layer, &hit_column);
									if (hit_layer == iL && hit_column == iC)
										tubeIsHit = kTRUE;
								}
								int col = iC;

								if (!tubeIsHit)
								{
									nMiss[iL][col] = nMiss[iL][col] + 1.0;
								}
								else
								{
									nHits[iL][col] = nHits[iL][col] + 1.0;

								}
							} // end if: track passes through gas volume
						}   // end if: check only active tubes
					}     // end for: column
				}       // end for: layer

				delete optTree;
				if (pass_event_check && (total_events_pass % 10000 == 1))
				{
					EDCanvas->cd();
					ed_event.AddTrack(Track(tp.slope(), tp.y_int()));
					ed->DrawEvent(ed_event, Geometry::getInstance(), NULL);
					ed->Clear();
					gSystem->Sleep(2000);
				}

				for (int iL = 0; iL < Geometry::MAX_TUBE_LAYER; ++iL)
				{
					for (int iC = 0; iC < Geometry::MAX_TUBE_COLUMN; ++iC)
					{
						if (nHits.at(iL).at(iC))
						{	

							tube_efficiency->SetBinContent(iC + 1, iL + 1, nHits[iL][iC] / (nHits[iL][iC] + nMiss[iL][iC]));
						}
					}
				}
			}
		}
		eff_canvas->Modified();
		eff_canvas->Update();
		residual_canvas->cd();
		residuals->Draw();
		residual_canvas->Update();
		trigger_rate_canvas->cd();
		trigger_rate_canvas->Update();

		cout<<"Finish update canvas"<<endl;
	} // End update canvas
}
#endif
