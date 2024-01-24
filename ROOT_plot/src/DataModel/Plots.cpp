/**
 * @file Plots.cpp
 *
 * @brief Data structure containing histograms for ATLAS cosmic ray testing.
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <vector>
#include <string>

#include "macros/GlobalIncludes.h"

#include "src/Geometry.cpp"
#include "src/Event.cpp"

//Muon Reconstruction includes
#include "src/TrackParam.cpp"
#include "src/RTParam.cpp"
#include "src/TimeCorrection.cpp"
#include "src/HitFinder.cpp"
#include "src/HitCluster.cpp"
#include "src/EventDisplay.cpp"
#include "src/CheckEvent.cpp" 

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

using namespace Muon;
using namespace std;

// NOTE: There should really only be one of those. DAQData holds it.
struct Plots {

	Plots();

	TH1F *                p_leading_time          ;
	TH1F *                p_trailing_time         ;

	vector<TH1F*>         p_hits_distribution     ;

	vector<vector<TH1F*>> p_tdc_time              ;
	vector<vector<TH1F*>> p_tdc_time_original     ;
	vector<vector<TH1F*>> p_tdc_time_corrected    ;
	vector<vector<TH1F*>> p_tdc_time_selected     ;
	vector<vector<TH1F*>> p_adc_time              ;


	vector<TH1F*>         p_tdc_tdc_time_original ;
	vector<TH1F*>         p_tdc_tdc_time_corrected;
	vector<TH1F*>         p_tdc_tdc_time_selected ;
	vector<TH1F*>         p_tdc_adc_time          ;
	vector<TH1F*>         p_tdc_channel           ;

	vector<TH2F*>         p_adc_vs_tdc            ;

	TGraph *              p_tdc_hit_rate_graph    ;

	TH2D *                hitByLC                 ;
	TH2D *                badHitByLC              ;
	TH2D *                goodHitByLC             ;
	TH2D *				  tube_efficiency;

	TH1D *				   residuals               ;

	RTParam * 			   rtp					   ;
	EventDisplay *		   ed					   ;
	Event 				   ed_event				   ;
	vector<Event>          ed_events			   ;

	TrackParam 			   tp					   ;
	TimeCorrection         tc					   ;
	int 				   ed_event_counter		   ;
	int 				   ed_event_rate		   ;

	vector<vector<double>> nHits;
	vector<vector<double>> nMiss;

	void binEvent(const Event &e);
	void clear();

};

Plots::Plots() {

	TString plot_name_buffer;

	p_leading_time  = new TH1F("leading time spectrum" , "leading time spectrum", 100, 0, 1000);
	p_trailing_time = new TH1F("trailing time spectrum", "trailing time spectrum", 100, 0, 1000);

	p_hits_distribution.reserve(Geometry::MAX_TUBE_LAYER);
	for(int layer_id = 0; layer_id != Geometry::MAX_TUBE_LAYER; ++layer_id) {

		plot_name_buffer.Form("layer_%d_hits_distribution", layer_id);
		p_hits_distribution.push_back(
				new TH1F(
					plot_name_buffer, 
					plot_name_buffer, 
					Geometry::MAX_TUBE_COLUMN, 
					-0.5, 
					Geometry::MAX_TUBE_COLUMN - 0.5
					)
				);

	}

	p_tdc_time              .resize(Geometry::MAX_TDC) ;
	p_tdc_time_original     .resize(Geometry::MAX_TDC) ;
	p_tdc_time_corrected    .resize(Geometry::MAX_TDC) ;
	p_tdc_time_selected     .resize(Geometry::MAX_TDC) ;
	p_adc_time              .resize(Geometry::MAX_TDC) ;

	p_tdc_tdc_time_original .reserve(Geometry::MAX_TDC);
	p_tdc_tdc_time_corrected.reserve(Geometry::MAX_TDC);
	p_tdc_tdc_time_selected .reserve(Geometry::MAX_TDC);
	p_tdc_adc_time          .reserve(Geometry::MAX_TDC);
	p_tdc_channel           .reserve(Geometry::MAX_TDC);

	p_adc_vs_tdc            .reserve(Geometry::MAX_TDC);

	nHits.resize(Geometry::MAX_TUBE_LAYER);
	nMiss.resize(Geometry::MAX_TUBE_LAYER);

	for (size_t i = 0; i < Geometry::MAX_TUBE_LAYER; ++i)
		{
			nHits[i].resize(Geometry::MAX_TUBE_COLUMN);
			nMiss[i].resize(Geometry::MAX_TUBE_COLUMN);
		}

	// TODO: Set up p_tdc_hit_rate_graph

	for(int tdc = 0; tdc < Geometry::MAX_TDC; ++tdc) {

		// TODO: Eventually restore this check and make sure everything still
		//       lines up well. Adapt canvas divisions to the actual number of
		//       active TDCs, and make sure we're converting TDC number into
		//       the correct index with each plot.
		// TODO: We will want to store some metadata with everything to e.g.
		//       map TDC number to a plot index.

		// if(!Geometry::getInstance()->IsActiveTDC(tdc)) continue;

		plot_name_buffer.Form("tdc_%d_adc_time_spectrum", tdc);
		p_tdc_adc_time.push_back(new TH1F(
					plot_name_buffer,
					plot_name_buffer,
					ADC_HIST_TOTAL_BIN,
					ADC_HIST_LEFT,
					ADC_HIST_RIGHT
					));
		p_tdc_adc_time.back()->GetXaxis()->SetTitle("time/ns");
		p_tdc_adc_time.back()->GetYaxis()->SetTitle("entries");

		string p_tdc_adc_time_corrected_title = "tdc_";
		p_tdc_adc_time_corrected_title += tdc;
		p_tdc_adc_time_corrected_title += "_tdc_time_spectrum_corrected";

		plot_name_buffer.Form("tdc_%d_tdc_time_spectrum_corrected", tdc);
		p_tdc_tdc_time_corrected.push_back(new TH1F(
					plot_name_buffer, 
					plot_name_buffer,
					TDC_HIST_TOTAL_BIN,
					TDC_HIST_LEFT,
					TDC_HIST_RIGHT
					));
		p_tdc_tdc_time_corrected.back()->GetXaxis()->SetTitle("time/ns");
		p_tdc_tdc_time_corrected.back()->GetYaxis()->SetTitle("entries");

		for(int channel = 0; channel < Geometry::MAX_TDC_CHANNEL; ++channel) {

			plot_name_buffer.Form("tdc_%d_channel_%d_tdc_time_spectrum_corrected", tdc, channel);
			p_tdc_time_corrected[tdc].push_back(new TH1F(
						plot_name_buffer,
						plot_name_buffer,
						TDC_HIST_TOTAL_BIN,
						TDC_HIST_LEFT,
						TDC_HIST_RIGHT
						));
			p_tdc_time_corrected[tdc].back()->GetXaxis()->SetTitle("time/ns");
			p_tdc_time_corrected[tdc].back()->GetYaxis()->SetTitle("entries");

			plot_name_buffer.Form("tdc_%d_channel_%d_tdc_time_spectrum", tdc, channel);
			p_tdc_time[tdc].push_back(new TH1F(
						plot_name_buffer,
						plot_name_buffer, 
						TDC_HIST_TOTAL_BIN,
						TDC_HIST_LEFT, 
						TDC_HIST_RIGHT
						));
			p_tdc_time[tdc].back()->GetXaxis()->SetTitle("time/ns");
			p_tdc_time[tdc].back()->GetYaxis()->SetTitle("entries");

			plot_name_buffer.Form("tdc_%d_channel_%d_adc_time_spectrum", tdc, channel);
			p_adc_time[tdc].push_back(new TH1F(
						plot_name_buffer, 
						plot_name_buffer,
						ADC_HIST_TOTAL_BIN, 
						ADC_HIST_LEFT, 
						ADC_HIST_RIGHT
						));
			p_adc_time[tdc].back()->GetXaxis()->SetTitle("time/ns");
			p_adc_time[tdc].back()->GetYaxis()->SetTitle("entries");

		}

	}

	hitByLC = new TH2D(
			"hitByLC", 
			"All hits on tubes (that passed clustering)", 
			54,
			-0.5,
			53.5,
			8,
			-0.5,
			7.5
			);
	hitByLC->SetStats(0);

	badHitByLC = new TH2D(
			"badHitByLC", 
			"Hits on tubes outside window (that passed clustering)", 
			54,
			-0.5,
			53.5,
			8,
			-0.5,
			7.5
			);
	badHitByLC->SetStats(0);

	goodHitByLC = new TH2D(
			"goodHitByLC", 
			"Hits on tubes inside window (that passed clustering)", 
			54,
			-0.5,
			53.5,
			8,
			-0.5,
			7.5
			);
	goodHitByLC->SetStats(0);


	residuals = new TH1D("residuals",
			"Residuals;Residual[#mum];Number of hits/2#mum", 
			500,
			-500, 
			500);
	residuals->SetStats(0);

	tube_efficiency = new TH2D("tube_efficiency",
			";Layer;Column",
			Geometry::MAX_TUBE_COLUMN, -0.5, Geometry::MAX_TUBE_COLUMN - 0.5,
			Geometry::MAX_TUBE_LAYER, -0.5, Geometry::MAX_TUBE_LAYER - 0.5);

	tp = TrackParam();
	tp.SetRT(rtp);
	tp.setVerbose(0);
	tp.setMaxResidual(1000000);

	ed_event_counter = 0;
	ed_event_rate    = 1000;
}

void Plots::binEvent(const Event &e) {

	// TODO: Event display
	// TODO: Go through DAQ.cpp and find everything we need to include
	// TODO: Make sure all plots print
	// TODO: Clean up and redesign UI
	// TODO: Split up binning and drawing so we can delay draw until we're done
	//       binning

	for(const Hit &hit : e.Hits()) {

		p_tdc_tdc_time_corrected[hit.TDC()]->Fill(hit.CorrTime());
		p_tdc_adc_time          [hit.TDC()]->Fill(hit.ADCTime ());

		// TODO: We can wait to make these until the end
		p_tdc_time_corrected[hit.TDC()][hit.Channel()]->Fill(hit.CorrTime ());
		p_tdc_time          [hit.TDC()][hit.Channel()]->Fill(hit.DriftTime());
		p_adc_time          [hit.TDC()][hit.Channel()]->Fill(hit.ADCTime  ());

		int hitL, hitC;
		Geometry::getInstance().GetHitLayerColumn(hit.TDC(), hit.Channel(), &hitL, &hitC);

		hitByLC->Fill(hitC, hitL);
		if(hit.CorrTime() < 0 || hit.CorrTime() > 400) { // TODO: Magic numbers

			badHitByLC->Fill(hitC, hitL);

		} else {

			goodHitByLC->Fill(hitC, hitL);

		}
		p_hits_distribution[hitL]->Fill(hitC);

		int panelIndex = hit.TDC() + 1;

	}
	rtp = new RTParam(Geometry::getInstance());
	tp.SetRT(rtp);
	tp.SetGEO(&Geometry::getInstance());
	TimeCorrection tc = TimeCorrection();

	ed_event = e;
	DoHitFinding(&ed_event, tc, 0, 0);
	DoHitClustering(&ed_event, Geometry::getInstance());
	bool pass_event_check = kTRUE;
	ed_event.SetPassCheck(pass_event_check);
	ed_event.CheckClusterTime();

	if (pass_event_check){
		ed_event_counter+=1;
		TTree *optTree = new TTree("optTree", "optTree");
		optTree->Branch("event", "Event", &ed_event);
		optTree->Fill();
		tp.setTarget(optTree);
		tp.setRangeSingle(0);
		tp.setIgnoreNone();
		tp.optimize();

		for (Cluster c : ed_event.Clusters())
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

						for (Hit h : e.Hits())
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


		delete optTree;
		//rtp->~RTParam();
		//tp.~TrackParam();
		//tc.~TimeCorrection();
		if (ed_event_counter%ed_event_rate==0)
		{
			ed_event.AddTrack(Track(tp.slope(), tp.y_int()));
			ed_events.push_back(ed_event);
		}
	}
	if (ed_event_counter%1000 == 0){
		cout<<"ed_event_counter:"<<ed_event_counter<<endl;
	}
}


void Plots::clear() {

	// TODO: We can also delete and remake everything to reset.

	p_leading_time ->Reset();
	p_trailing_time->Reset();

	for(TH1F *h : p_hits_distribution) { h->Reset(); }

	for(vector<TH1F*> &v : p_tdc_time) {

		for(TH1F *h : v) h->Reset();

	}

	for(vector<TH1F*> &v : p_tdc_time_original) {

		for(TH1F *h : v) h->Reset();

	}

	for(vector<TH1F*> &v : p_tdc_time_corrected) {

		for(TH1F *h : v) h->Reset();

	}

	for(vector<TH1F*> &v : p_tdc_time_selected) {

		for(TH1F *h : v) h->Reset();

	}

	for(vector<TH1F*> &v : p_adc_time) {

		for(TH1F *h : v) h->Reset();

	}            ;


	for(TH1F *h : p_tdc_tdc_time_original ) { h->Reset(); }
	for(TH1F *h : p_tdc_tdc_time_corrected) { h->Reset(); }
	for(TH1F *h : p_tdc_tdc_time_selected ) { h->Reset(); }
	for(TH1F *h : p_tdc_adc_time          ) { h->Reset(); }
	for(TH1F *h : p_tdc_channel           ) { h->Reset(); }

	for(TH2F *h : p_adc_vs_tdc            ) { h->Reset(); }

	delete p_tdc_hit_rate_graph;
	p_tdc_hit_rate_graph = nullptr;
	// TODO: Remake the hit rate graph


	hitByLC    ->Reset();
	badHitByLC ->Reset();
	goodHitByLC->Reset();
	residuals   ->Reset();
	tube_efficiency ->Reset();
}
