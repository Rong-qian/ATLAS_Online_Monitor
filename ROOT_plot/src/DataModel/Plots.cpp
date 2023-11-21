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

using namespace Muon;
using namespace std;

// TODO: Pull out initialization from default constructor, copy constructor,
//       and copy assignment

// NOTE: There should really only be one of those. DAQData holds it.
struct Plots {

	Plots();

	Plots         (const Plots &other);
	void operator=(const Plots &other) = delete; // TODO: It's a little bit
	                                             //       weird to have a
	                                             //       copy constructor but
	                                             //       not copy assignment

	TH1F *                 p_leading_time          ;
	TH1F *                 p_trailing_time         ;

	vector<TH1F*>          p_hits_distribution     ;

	vector<vector<TH1F*>>  p_tdc_time              ;
	vector<vector<TH1F*>>  p_tdc_time_original     ;
	vector<vector<TH1F*>>  p_tdc_time_corrected    ;
	vector<vector<TH1F*>>  p_tdc_time_selected     ;
	vector<vector<TH1F*>>  p_adc_time              ;


	vector<TH1F*>          p_tdc_tdc_time_original ;
	vector<TH1F*>          p_tdc_tdc_time_corrected;
	vector<TH1F*>          p_tdc_tdc_time_selected ;
	vector<TH1F*>          p_tdc_adc_time          ;
	vector<TH1F*>          p_tdc_channel           ;

	vector<TH2F*>          p_adc_vs_tdc            ;

	vector<vector<double>> p_tdc_hit_rate          ;

	TH2D *                 hitByLC                 ;
	TH2D *                 badHitByLC              ;
	TH2D *                 goodHitByLC             ; 

	void binEvent(const Event &e);

	void clear();

};

Plots::Plots(const Plots &other) {

	p_leading_time  = dynamic_cast<TH1F*>(other.p_leading_time ->Clone());
	p_trailing_time = dynamic_cast<TH1F*>(other.p_trailing_time->Clone());

	p_hits_distribution.reserve(Geometry::MAX_TUBE_LAYER);
	for(int layer_id = 0; layer_id < Geometry::MAX_TUBE_LAYER; ++layer_id) {

		p_hits_distribution.push_back(dynamic_cast<TH1F*>(other.p_hits_distribution[layer_id]->Clone()));

	}

	p_tdc_time              .resize (Geometry::MAX_TDC) ;
	p_tdc_time_original     .resize (Geometry::MAX_TDC) ;
	p_tdc_time_corrected    .resize (Geometry::MAX_TDC) ;
	p_tdc_time_selected     .resize (Geometry::MAX_TDC) ;
	p_adc_time              .resize (Geometry::MAX_TDC) ;

	p_tdc_tdc_time_original .reserve(Geometry::MAX_TDC);
	p_tdc_tdc_time_corrected.reserve(Geometry::MAX_TDC);
	p_tdc_tdc_time_selected .reserve(Geometry::MAX_TDC);
	p_tdc_adc_time          .reserve(Geometry::MAX_TDC);
	p_tdc_channel           .reserve(Geometry::MAX_TDC);

	p_adc_vs_tdc            .reserve(Geometry::MAX_TDC);

	p_tdc_hit_rate          .resize (Geometry::MAX_TDC);

	for(vector<double> &vec : p_tdc_hit_rate) {

		vec.resize(Geometry::MAX_TDC_CHANNEL);

	}

	for(int tdc = 0; tdc < Geometry::MAX_TDC; ++tdc) {

		p_tdc_adc_time          .push_back(dynamic_cast<TH1F*>(other.p_tdc_adc_time[tdc]          ->Clone()));
		p_tdc_tdc_time_corrected.push_back(dynamic_cast<TH1F*>(other.p_tdc_tdc_time_corrected[tdc]->Clone()));

		for(int channel = 0; channel < Geometry::MAX_TDC_CHANNEL; ++channel) {

			p_tdc_time_corrected[tdc].push_back(dynamic_cast<TH1F*>(other.p_tdc_time_corrected[tdc][channel]->Clone()));
			p_tdc_time[tdc]          .push_back(dynamic_cast<TH1F*>(other.p_tdc_time[tdc][channel]          ->Clone()));
			p_adc_time[tdc]          .push_back(dynamic_cast<TH1F*>(other.p_adc_time[tdc][channel]          ->Clone()));

		}

	}

	hitByLC     = dynamic_cast<TH2D*>(other.hitByLC    ->Clone());
	badHitByLC  = dynamic_cast<TH2D*>(other.badHitByLC ->Clone());
	goodHitByLC = dynamic_cast<TH2D*>(other.goodHitByLC->Clone());

}

Plots::Plots() {

	TString plot_name_buffer;

	p_leading_time  = new TH1F("leading time spectrum" , "leading time spectrum", 100, 0, 1000);
	p_trailing_time = new TH1F("trailing time spectrum", "trailing time spectrum", 100, 0, 1000);

	p_hits_distribution.reserve(Geometry::MAX_TUBE_LAYER);
	for(int layer_id = 0; layer_id < Geometry::MAX_TUBE_LAYER; ++layer_id) {

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

	p_tdc_time              .resize (Geometry::MAX_TDC) ;
	p_tdc_time_original     .resize (Geometry::MAX_TDC) ;
	p_tdc_time_corrected    .resize (Geometry::MAX_TDC) ;
	p_tdc_time_selected     .resize (Geometry::MAX_TDC) ;
	p_adc_time              .resize (Geometry::MAX_TDC) ;

	p_tdc_tdc_time_original .reserve(Geometry::MAX_TDC);
	p_tdc_tdc_time_corrected.reserve(Geometry::MAX_TDC);
	p_tdc_tdc_time_selected .reserve(Geometry::MAX_TDC);
	p_tdc_adc_time          .reserve(Geometry::MAX_TDC);
	p_tdc_channel           .reserve(Geometry::MAX_TDC);

	p_adc_vs_tdc            .reserve(Geometry::MAX_TDC);

	p_tdc_hit_rate          .resize (Geometry::MAX_TDC);

	for(vector<double> &vec : p_tdc_hit_rate) {

		vec.resize(Geometry::MAX_TDC_CHANNEL);

	}

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

}

void Plots::binEvent(const Event &e) {

	// TODO: Event display
	// TODO: Go through DAQ.cpp and find everything we need to include
	// TODO: Make sure all plots print
	// TODO: Clean up and redesign UI
	// TODO: Noise rate display

	for(const Hit &hit : e.Hits()) {

		p_tdc_tdc_time_corrected[hit.TDC()]->Fill(hit.CorrTime());
		p_tdc_adc_time          [hit.TDC()]->Fill(hit.ADCTime ());

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

	}


	for(TH1F *h : p_tdc_tdc_time_original ) { h->Reset(); }
	for(TH1F *h : p_tdc_tdc_time_corrected) { h->Reset(); }
	for(TH1F *h : p_tdc_tdc_time_selected ) { h->Reset(); }
	for(TH1F *h : p_tdc_adc_time          ) { h->Reset(); }
	for(TH1F *h : p_tdc_channel           ) { h->Reset(); }

	for(TH2F *h : p_adc_vs_tdc            ) { h->Reset(); }

	for(vector<double> &vec : p_tdc_hit_rate) {

		vec.clear();
		vec.resize(Geometry::MAX_TDC_CHANNEL);

	}


	hitByLC    ->Reset();
	badHitByLC ->Reset();
	goodHitByLC->Reset();

}