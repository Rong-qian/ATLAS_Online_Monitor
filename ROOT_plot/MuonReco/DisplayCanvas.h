/*******************************************************************************
  file name: DisplayCanvas.h
  author: Rongqian Qian
  created: 10/31/2023
  last modified: 10/31/2023

  description:
  -Header file for DisplayerCanvas.cpp which create the online monitor canvas and update when receiving data

*******************************************************************************/

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

#include "src/DataModel/DAQData.cpp"

#define SAVE_TRACKS_OUT_OF_ROOT // comment this line if you don't need to save plots out of rootfile
#define ERROR_WORD 29
#define HEADER_WORD 31
#define TRAILER_WORD 30
#define SPEEDFACTOR 1
#define NEWTDC_NUMBER 17
#define WIDTH_RES 1

#ifndef DISPLAYCANVAS_H
#define DISPLAYCANVAS_H

namespace Muon {

  class DisplayCanvas
  {
  public:
      DisplayCanvas();
      void UpdateCanvas(DAQData daq_data);
      void DisplayPrints();
  private:    
      TCanvas*EDCanvas;
      Geometry geo;
      RTParam rtp;
      EventDisplay ed;
      TrackParam tp;
      int total_events_pass;
      TimeCorrection tc;
  };
}
#endif