
/*******************************************************************************
  file name: DecodeRawData.cxx
  author: Zhe Yang
  created: 01/25/2019
  last modified: 04/26/2019

  description:
  -Decode .raw data from HPTDC and save data to ntuple

  remark:
  -Learned basic decode method from Shuzhou Zhang, redeveloped and added new
  function for new HPTDC data format.

*******************************************************************************/
#include "macros/GlobalIncludes.h"

#include <stdio.h>
#include <iostream>
#include <bitset>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TString.h"
// Geometry includes
// #include "src/CheckEvent.cpp"

// Muon Reconstruction includes
#include "src/Signal.cpp"
#include "src/EventID.cpp"
#include "src/Event.cpp"
#include "src/EventDisplay.cpp"
#include "src/Geometry.cpp"
#include "src/Hit.cpp"
#include "src/Cluster.cpp"
#include "src/TimeCorrection.cpp"
#include "src/HitFinder.cpp"
// #include "src/HitFinder_BCID.cpp"
#include "src/HitCluster.cpp"
#include "src/Noisecut.cpp"

// New Reco Includes
#include "src/T0Fit.h"
#include "src/TrackParam.cpp"
#include "src/Observable.cpp"
#include "src/Parameterization.cpp"
#include "src/RTParam.cpp"
#include "src/ResolutionResult.cpp"

// Old Reco Includes
#include "src/ConfigParser.cxx"
#include "src/Track.cpp"
#include "src/IOUtility.cpp"
#include "src/Optimizer.cpp"

#include "src/udp_client_init.c"
#include "src/channel_packet.c"
#include "src/sockhelp.c"

// #define DECODER_DEBUG // comment this line when debugging information is not needed
// #define SET_MAXWORDS // comment this line if you want to decode the whole data words
#define SAVE_TRACKS_OUT_OF_ROOT // comment this line if you don't need to save plots out of rootfile
#define ERROR_WORD 29
#define HEADER_WORD 31
#define TRAILER_WORD 30
#define SPEEDFACTOR 1
#define NEWTDC_NUMBER 17
#define WIDTH_RES 1

#ifndef ONLINE_DAQ_MONITOR
#define ONLINE_DAQ_MONITOR

using namespace std;
using namespace Muon;

class DAQ_monitor
{
public:
  DAQ_monitor(short portno_input, int bisno = 1);
  void error(const char *msg);
  void tcp_server_setup(in_addr_t server_ip_int);
  void DecodeOnline(TString filename);

private:
  TDirectory *event_track[2];
  char track_group_name[128];
  std::vector<TH1F *> p_tdc_adc_time;
  std::vector<TH1F *> p_tdc_tdc_time_corrected;
  std::vector<std::vector<TH1F *>> p_tdc_chnl_adc_time;
  std::vector<std::vector<TH1F *>> p_tdc_chnl_tdc_time_corrected;
  std::vector<std::vector<TH1F *>> p_tdc_chnl_adc_time_raw;
  std::vector<std::vector<TH1F *>> p_tdc_chnl_tdc_time_corrected_raw;
  std::vector<std::vector<double>> p_tdc_hit_rate;
  std::vector<std::vector<double>> nHits;
  std::vector<std::vector<double>> nMiss;
  std::vector<std::vector<Bool_t>> first_signal_flag;
  
  std::vector<double> p_tdc_hit_rate_x;
  std::vector<bool> isDrawn;
  Geometry geo;
  TH2D *tube_efficiency;
  std::vector<TH1F *> p_tdc_hit_rate_graph;
  EventDisplay *ed;
  ResolutionResult *rr;
  TrackParam tp;
  RTParam *rtp;
  TCanvas *adc_canvas, *tdc_canvas, *rate_canvas, *trigger_rate_canvas, *residual_canvas, *EDCanvas, *eff_canvas;
  short tcp_portno;
  int sockfd, newsockfd, udp_sock_fd;
  struct sockaddr_in udp_servaddr;
  in_addr_t tcp_server_ip;
  int pad_num;
  unsigned int buffer[4096];
  int sockReadCount, bytes_recv, total_bytes_recv;
  ofstream oFile;
  FILE *fp_rate_File;
  char oFile_name[30];
  char filename_time[30];
  ifstream data_in_flow;
  TFile *p_output_rootfile;
  TH2D *hitByLC, *badHitByLC, *goodHitByLC;
  unsigned long total_triggers, total_events, total_triggers_pass, total_events_pass;
  unsigned long total_signals, total_signals_pass, total_events_fail, event_print, event_signals;
  int hitL, hitC;
  int current_track_group;
  int temp_track_group;
  bool pass_event_check;
  uint64_t word;
  char readbuff[5] = {0, 0, 0, 0, 0};
  unsigned int header_type;
  vector<Signal> sigVec;
  bitset<4> header;
  bitset<1> bit1;
  bitset<3> bit3;
  bitset<3> bit4;
  bitset<5> bit5;
  // EventID currEventID, prevEventID;
  unsigned int current_event_ID = 0;
  unsigned int previous_event_ID = -1;
  unsigned int tdc_header_eventID_current = 0;
  unsigned int tdc_header_eventID_previous = -1;
  unsigned int header_count_Event_trailer;
  unsigned int trailer_count_Event_trailer;
  unsigned int header_count_TDC_data;
  unsigned int trailer_count_TDC_data;
  int error_word_count;
  int unexpected_data;
  unsigned int error_flag;
  unsigned int error_channel;
  vector<Signal> trigVec;
  Signal sig, header_sig, trailer_sig;
  Event event, event_raw, ed_event;
  TTree *eTree = new TTree("eTree", "eTree");
  stringstream *filar_1, *filar_2;
  time_t start_time, current_time;
  double DAQ_time;
  int status;
  int tracking_evt_no = 0;
  TH1D *residuals;
};

DAQ_monitor::DAQ_monitor(short portno_input, int bisno /*=1*/)
{
  tcp_portno = portno_input;
  geo = Geometry();
  TString filename = "run_20231003_152805.dat";
  TString fn = TString(filename);
  int runN = ((TObjString *)(TString(fn(3, 256)).Tokenize("_")->At(0)))->String().Atoi();
  geo.SetRunN(runN);
  ConfigParser cp;
  if (bisno == 1)
  {
    cp = ConfigParser("/atlas/data18a/rongqian/update/sMDT_MiniDAQ/smdt-reco/smdt-reco/conf/BIS1.conf");
  }
  else
  {
    cp = ConfigParser("/atlas/data18a/rongqian/update/sMDT_MiniDAQ/smdt-reco/smdt-reco/conf/BIS2-6.conf");
  }

  residuals = new TH1D("residuals", "Residuals;Residual[#mum];Number of hits/2#mum", 500, -500, 500);
  rtp = new RTParam(geo);
  ed = new EventDisplay();
  tp = TrackParam();
  tp.SetGEO(&geo);

  current_track_group = 0;
  int temp_track_group;
  word = 0; // sizeof(long)=8, which equals uint64_t.
  // unsigned int word;
  char readbuff[5] = {0, 0, 0, 0, 0};
  int unexpected_data = 0;
  unsigned int header_type;
  unsigned int current_event_ID = 0;
  unsigned int previous_event_ID = -1;
  unsigned int tdc_header_eventID_current = 0;
  unsigned int tdc_header_eventID_previous = -1;

  // EventID currEventID;
  // EventID prevEventID = EventID(0x00000000);

  // ed->SetRT(rtp);  // TODO: GET A WORKING RT ONLINE
  TString h_name;

  p_tdc_adc_time.resize(Geometry::MAX_TDC);
  p_tdc_tdc_time_corrected.resize(Geometry::MAX_TDC);
  p_tdc_hit_rate_x.resize(Geometry::MAX_TDC_CHANNEL);
  isDrawn.resize(Geometry::MAX_TDC);

  p_tdc_chnl_adc_time.resize(Geometry::MAX_TDC);
  p_tdc_chnl_adc_time_raw.resize(Geometry::MAX_TDC);
  p_tdc_chnl_tdc_time_corrected.resize(Geometry::MAX_TDC);
  p_tdc_chnl_tdc_time_corrected_raw.resize(Geometry::MAX_TDC);
  p_tdc_hit_rate.resize(Geometry::MAX_TDC);
  p_tdc_hit_rate_graph.resize(Geometry::MAX_TDC);
  first_signal_flag.resize(Geometry::MAX_TDC);
  nHits.resize(Geometry::MAX_TUBE_LAYER);
  nMiss.resize(Geometry::MAX_TUBE_LAYER);

  for (size_t i = 0; i < Geometry::MAX_TUBE_LAYER; ++i)
  {
    nHits[i].resize(Geometry::MAX_TUBE_COLUMN);
    nMiss[i].resize(Geometry::MAX_TUBE_COLUMN);
  }

  for (size_t i = 0; i < Geometry::MAX_TDC; ++i)
  {
    p_tdc_chnl_adc_time[i].resize(Geometry::MAX_TDC_CHANNEL);
    p_tdc_chnl_adc_time_raw[i].resize(Geometry::MAX_TDC_CHANNEL);
    p_tdc_chnl_tdc_time_corrected[i].resize(Geometry::MAX_TDC_CHANNEL);
    p_tdc_chnl_tdc_time_corrected_raw[i].resize(Geometry::MAX_TDC_CHANNEL);
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
  for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
  {
    int local_tdc_id = tdc_id % 18;
    h_name.Form("tdc_%d_tdc_time_spectrum_corrected", tdc_id);
    p_tdc_tdc_time_corrected[tdc_id] = new TH1F(h_name, TString::Format("tdc_%d_tdc_time_spectrum;time/ns;entries", local_tdc_id), TDC_HIST_TOTAL_BIN, TDC_HIST_LEFT, TDC_HIST_RIGHT);
    if (tdc_id >= 18)
      p_tdc_tdc_time_corrected[tdc_id]->SetLineColor(kRed);

    h_name.Form("tdc_%d_adc_time_spectrum", tdc_id);
    p_tdc_adc_time[tdc_id] = new TH1F(h_name, TString::Format("tdc_%d_adc_time_spectrum;time/ns;entries", local_tdc_id), ADC_HIST_TOTAL_BIN, ADC_HIST_LEFT, ADC_HIST_RIGHT);
    cout << "creating:" << tdc_id << endl;
    if (tdc_id >= 18)
      p_tdc_adc_time[tdc_id]->SetLineColor(kRed);

    h_name.Form("tdc_%d_hit_rate", tdc_id);
    p_tdc_hit_rate_graph[tdc_id] = new TH1F(h_name, TString::Format("tdc_%d_tdc_hit_rate;Channel No.;Rate(Hz)", local_tdc_id), 24, -0.5, 23.5);
    if (tdc_id < 18)
      p_tdc_hit_rate_graph[tdc_id]->SetFillColor(4);
    else
      p_tdc_hit_rate_graph[tdc_id]->SetFillColor(kRed);

    p_tdc_hit_rate_graph[tdc_id]->SetBarWidth(0.4);
    p_tdc_hit_rate_graph[tdc_id]->SetBarOffset(0.1 + 0.4 * (tdc_id < 18));
    p_tdc_hit_rate_graph[tdc_id]->SetStats(0);

    p_tdc_hit_rate_graph[tdc_id]->SetMaximum(1);
    p_tdc_hit_rate_graph[tdc_id]->SetMinimum(0);

    for (int tdc_chnl_id = 0; tdc_chnl_id != Geometry::MAX_TDC_CHANNEL; tdc_chnl_id++)
    {
      h_name.Form("tdc_%d_chnl_%d_adc_time_spectrum", tdc_id, tdc_chnl_id);
      p_tdc_chnl_adc_time[tdc_id][tdc_chnl_id] = new TH1F(h_name, h_name, ADC_HIST_TOTAL_BIN, ADC_HIST_LEFT, ADC_HIST_RIGHT);
      h_name.Form("tdc_%d_chnl_%d_adc_time_raw_spectrum", tdc_id, tdc_chnl_id);
      p_tdc_chnl_adc_time_raw[tdc_id][tdc_chnl_id] = new TH1F(h_name, h_name, ADC_HIST_TOTAL_BIN, ADC_HIST_LEFT, ADC_HIST_RIGHT);

      h_name.Form("tdc_%d_chnl_%d_tdc_time_spectrum_corrected", tdc_id, tdc_chnl_id);
      p_tdc_chnl_tdc_time_corrected[tdc_id][tdc_chnl_id] = new TH1F(h_name, h_name, TDC_HIST_TOTAL_BIN, TDC_HIST_LEFT, TDC_HIST_RIGHT);
      h_name.Form("tdc_%d_chnl_%d_tdc_time_spectrum_raw_corrected", tdc_id, tdc_chnl_id);
      p_tdc_chnl_tdc_time_corrected_raw[tdc_id][tdc_chnl_id] = new TH1F(h_name, h_name, TDC_HIST_TOTAL_BIN, TDC_HIST_LEFT, TDC_HIST_RIGHT);
    }

  } // end for: all TDC

  tube_efficiency = new TH2D("tube_efficiency", ";Layer;Column",
                             Geometry::MAX_TUBE_COLUMN, -0.5, Geometry::MAX_TUBE_COLUMN - 0.5,
                             Geometry::MAX_TUBE_LAYER, -0.5, Geometry::MAX_TUBE_LAYER - 0.5);

  adc_canvas = new TCanvas("c1", "ADC Plots", 0, 0, 2160, 750);
  adc_canvas->Divide(6, 2);
  tdc_canvas = new TCanvas("c2", "TDC Plots", 0, 750, 2160, 750);
  tdc_canvas->Divide(6, 2);
  rate_canvas = new TCanvas("c3", "Hit Rate Plots", 2160, 0, 1800, 750);
  rate_canvas->Divide(6, 2);
  trigger_rate_canvas = new TCanvas("c4", "Trigger Board", 1440, 750, 400, 300);
  residual_canvas = new TCanvas("c5", "Residuals", 2100, 900, 400, 300);
  EDCanvas = new TCanvas("c6", "Event Display", 2700, 900, 800, 800);
  eff_canvas = new TCanvas("C7", "Efficiency", 2300, 900, 400, 300);

  printf("Canvases created and divided.\n");
  for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
  {
    if (geo.IsActiveTDC(tdc_id) || tdc_id == geo.TRIGGER_MEZZ)
    {
      int local_tdc_id = tdc_id % 18;
      if (tdc_id == geo.TRIGGER_MEZZ)
      {

        trigger_rate_canvas->cd();
        p_tdc_hit_rate_graph[tdc_id]->Draw("B");
      }
      else
      {
        pad_num = 6 * ((local_tdc_id + 1) % 2) + (local_tdc_id / 2) + 1;
        TString opts;
        if (isDrawn[tdc_id])
          opts = "same";
        else
          opts = "";

        adc_canvas->cd(pad_num);
        cout << tdc_id << "Setting:"
             << ":" << pad_num << endl;
        p_tdc_adc_time[tdc_id]->Draw(opts);
        tdc_canvas->cd(pad_num);
        p_tdc_tdc_time_corrected[tdc_id]->Draw(opts);
        rate_canvas->cd(pad_num);
        // gPad->SetLogy();
        p_tdc_hit_rate_graph[tdc_id]->Draw(opts + " B");
        // if (gSystem->ProcessEvents()) break;
        isDrawn[tdc_id] = 1;
      }
    }
  }
  pad_num = 12;
  adc_canvas->cd();
  adc_canvas->Modified();
  adc_canvas->Update();
  tdc_canvas->cd();
  tdc_canvas->Modified();
  tdc_canvas->Update();
  rate_canvas->cd();
  rate_canvas->Modified();
  rate_canvas->Update();
  residual_canvas->cd();
  residuals->Draw();
  residual_canvas->Modified();
  residual_canvas->Update();
  eff_canvas->cd();
  tube_efficiency->Draw("colz");
  eff_canvas->Modified();
  eff_canvas->Update();
  gSystem->ProcessEvents();
  printf("Canvases updated.\n");
  /*
  // TCP SERVER SETUP
  tcp_server_ip = INADDR_ANY;
  // tcp_server_ip = inet_addr("141.211.96.10");
  // tcp_server_ip = inet_addr("141.213.133.230");
  tcp_server_setup(tcp_server_ip);
  // TCP SERVER SETUP DONE


  // UDP SERVER SETUP
  udp_sock_fd = udp_client_init(UDP_PORT);

  in_addr_t udp_server_ip= inet_addr("127.0.0.1");

  memset(&udp_servaddr, 0, sizeof(udp_servaddr));
  // Filling server information
  udp_servaddr.sin_family = AF_INET;
  udp_servaddr.sin_port = htons(UDP_PORT);
  udp_servaddr.sin_addr.s_addr = udp_server_ip;
  // UDP SERVER SETUP DONE

  time_t sys_time;
  struct tm * timeinfo;
  sys_time = time(0);
  timeinfo = localtime(&sys_time);
  memset(filename_time, 0, sizeof(filename_time));

  strftime(filename_time, 30, "%Y%m%d_%H%M%S", timeinfo);
  sprintf(oFile_name,"./data/%s.dat",filename_time);
  oFile.open(oFile_name, ios::out | ios::binary);
  printf("File %s opened for raw data recording.\n",oFile_name);
  */
  // data_in_flow.open(oFile_name);
  // data_in_flow.open("../data/run_20231003_152805.dat");

  p_output_rootfile = new TFile("output.root", "RECREATE");

  hitByLC = new TH2D("hitByLC", "All hits on tubes (that passed clustering)", 54, -0.5, 53.5, 8, -0.5, 7.5);
  hitByLC->SetStats(0);
  badHitByLC = new TH2D("badHitByLC", "Hits on tubes outside window (that passed clustering)", 54, -0.5, 53.5, 8, -0.5, 7.5);
  badHitByLC->SetStats(0);
  goodHitByLC = new TH2D("goodHitByLC", "Hits on tubes inside window (that passed clustering)", 54, -0.5, 53.5, 8, -0.5, 7.5);
  goodHitByLC->SetStats(0);

  total_triggers = 0;
  total_events = 0;
  total_triggers_pass = 0;
  total_events_pass = 0;
  total_signals = 0;
  total_signals_pass = 0;
  total_events_fail = 0;
  event_signals = 0;
  event_print = 100;

  event = Event();
  event_raw = Event();
  // prevEventID = EventID(0x00000000);
  eTree->Branch("event", "Event", &event);

  filar_1 = new stringstream;
  filar_2 = new stringstream;
  cout << "Processing..." << endl;
}

void DAQ_monitor::tcp_server_setup(in_addr_t server_ip_int)
{
  socklen_t clilen;
  char server_ip_address[INET_ADDRSTRLEN], client_ip_address[INET_ADDRSTRLEN];
  struct sockaddr_in serv_addr, cli_addr;

  // Creating socket file descriptor
  if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == 0)
  {
    perror("socket failed");
    exit(EXIT_FAILURE);
  }
  // Forcefully attaching socket to the port number from main() argument
  // error for local IP
  int opt = 1;
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)))
  {
    perror("setsockopt");
    exit(EXIT_FAILURE);
  }

  bzero((char *)&serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = server_ip_int;
  serv_addr.sin_port = htons(tcp_portno);
  inet_ntop(AF_INET, &(serv_addr.sin_addr), server_ip_address, INET_ADDRSTRLEN);
  printf("Waiting for client to connect to IP: %s Port: %u\n", server_ip_address, tcp_portno);

  // Forcefully attaching socket
  if (bind(sockfd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0)
  {
    perror("bind failed");
    exit(EXIT_FAILURE);
  }
  if (listen(sockfd, 5) < 0)
  {
    perror("listen");
    exit(EXIT_FAILURE);
  }

  clilen = sizeof(cli_addr);
  if ((newsockfd = accept(sockfd, (struct sockaddr *)&cli_addr, (socklen_t *)&clilen)) < 0)
  {
    perror("accept");
    exit(EXIT_FAILURE);
  }
  time_t local_time;
  local_time = time(0);
  inet_ntop(AF_INET, &(cli_addr.sin_addr), client_ip_address, INET_ADDRSTRLEN);
  printf("Connected to client IP:%s\n at %s\n", client_ip_address, ctime(&local_time));
  return;
}

void DAQ_monitor::error(const char *msg)
{
  perror(msg);
  exit(1);
}

/*
 * Get File Name from a Path with or without extension
 */
std::string getFileName(std::string filePath, bool withExtension = true, char seperator = '/')
{
  // Get last dot position
  std::size_t dotPos = filePath.rfind('.');
  std::size_t sepPos = filePath.rfind(seperator);
  if (sepPos != std::string::npos)
  {
    return filePath.substr(sepPos + 1, filePath.size() - (withExtension || dotPos != std::string::npos ? 1 : dotPos));
  }
  return "";
}

void DAQ_monitor::DecodeOnline(TString filename = "run_20231003_152805.dat")
{

  //gROOT->SetBatch(kTRUE); // set to batch mode to inprove the speed
  int maxEventCount = 100000000;
  // int maxEventCount = 100;
  // gStyle->SetOptStat(10); //only print entries
  gStyle->SetOptStat(1110); // std, mean, entris, name printed
  // gStyle->SetTitleX(999.);//hist no title
  // gStyle->SetTitleY(999.);
  gStyle->SetStatY(0.9);
  // Set y-position (fraction of pad size)
  gStyle->SetStatX(0.9);
  // Set x-position (fraction of pad size)
  gStyle->SetStatW(0.25);
  // Set width of stat-box (fraction of pad size)
  gStyle->SetStatH(0.25);
  // Set height of stat-box (fraction of pad size)

  double match_window = 1.5;
  // open input file
  /*
  total_bytes_recv = 0;
  bzero(buffer,sizeof(buffer));
  bytes_recv = sock_read(newsockfd, (char *) buffer, sizeof(buffer));
  total_bytes_recv += bytes_recv;
  sockReadCount = 1;
  printf("\nReceiving data...\n");
  printf("Received message %i\n",sockReadCount);
  */
  // time(&start_time);
  int iter = 0;
  bool anyTrig = 0;

  int loop_file = 0;

  while (loop_file < 3)
  {

    loop_file += 1;
    gSystem->Sleep(2000);
    cout << "Test loop_file:" << loop_file << endl;
    iter++;
    // std::cout << "CSM NUMBER: " << (buffer[20]&3) << " " << (buffer[21]&3) << " " << (buffer[22]&3) << " " << (buffer[23]&3) << std::endl;
    // std::cout << " bytes:     " << bytes_recv << std::endl;

    // bytes_recv = sock_read(newsockfd, (char *) buffer, sizeof(buffer));

    // total_bytes_recv += bytes_recv;
    // sockReadCount++;

    data_in_flow.tellg(); // Needed to ensure reading continues
    int nloop = 0;
    oFile.write((const char *)buffer, bytes_recv);
    TString input_filename = "/atlas/data18a/rongqian/update/phase2_MiniDAQ_2/ROOT_plot/data/";
    // TString input_filename = "data/";

    // if (loop_file ==2){
    //   filename = "run_20231003_152805_1.dat";
    // }
    // if (loop_file ==3){
    //   filename = "run_20231003_152805_2.dat";
    // }

    TString fn = TString(filename);
    input_filename += filename;
    std::cout << input_filename.Data() << std::endl;

    data_in_flow.open(input_filename.Data());
    data_in_flow.seekg(0, data_in_flow.end);
    int data_in_flow_length = data_in_flow.tellg(); // get file size
    if (data_in_flow_length <= 0)
    {
      printf("file name incorrect!\n");
      printf("file size = %d\n", data_in_flow_length);
      return 1;
    }
    else
    {
      printf("file size = %d\n", data_in_flow_length);
    }
    data_in_flow.seekg(0, data_in_flow.beg);

 
    static TimeCorrection tc = TimeCorrection();
    // static EventDisplay   ed = EventDisplay();

    Noisecut ncut = Noisecut();


    cout << "Processing..." << endl;

    int readbytes = 0;
    while (data_in_flow.read(readbuff, 5) && total_events < maxEventCount)
    {
      // while (data_in_flow.read((char *) &word, 4)  && nloop<maxEventCount) {
      nloop++;
      word = (((uint64_t)readbuff[4]) & 0xff) +
             ((((uint64_t)readbuff[3]) & 0xff) << 8) +
             ((((uint64_t)readbuff[2]) & 0xff) << 16) +
             ((((uint64_t)readbuff[1]) & 0xff) << 24) +
             ((((uint64_t)readbuff[0]) & 0xff) << 32);
      // networking uses big-endian while processor uses small-endian, so a byte-swap is used here.

      // word = readbuff[0] + readbuff[0]*256 + readbuff[0]*65536 + readbuff[0]*16777216 + readbuff[0]*4294967296;
      header = word >> 36; // get the four bits header of this word
      header_type = static_cast<unsigned int>((header.to_ulong()));

      if (header_type == Signal::TRAILER)
      {
        // analyze data if we reached a trailer for an event
        if (total_events % event_print == 0)
        {
          std::cout << "Processing Event " << total_events << std::endl;
          if (TMath::Floor(TMath::Log10(total_events)) > TMath::Floor(TMath::Log10(event_print)))
            event_print *= 10;
        }
        trailer_sig = Signal(word);
        if (trailer_sig.TrailerEID() != header_sig.HeaderEID())
        {
          printf("Warning! Event ID mismatch! HEADER ID= %d, TRAILER ID = %d\n", header_sig.HeaderEID(), trailer_sig.TrailerEID());
          sigVec.clear();
          continue;
        }
        bit4 = word >> 32; // get the TDC header count in Event trailer
        header_count_Event_trailer = static_cast<unsigned int>((bit4.to_ulong()));
        if (header_count_Event_trailer != header_count_TDC_data)
        {
          printf("Warning! %d header(s) found in data, event trailer indicates %d!\n", header_count_TDC_data, header_count_Event_trailer);
        }
        bit4 = word >> 28; // get the TDC trailer count in Event trailer
        trailer_count_Event_trailer = static_cast<unsigned int>((bit4.to_ulong()));
        if (trailer_count_Event_trailer != trailer_count_TDC_data)
        {
          printf("Warning! %d trailer(s) found in data, event trailer indicates %d!\n", header_count_TDC_data, header_count_Event_trailer);
        }
        bit1 = word >> 15; // get the TDC header count error flag in Event trailer
        unsigned int error_flag_tmp;
        error_flag_tmp = static_cast<unsigned int>((bit1.to_ulong()));
        if (error_flag_tmp)
        {
          printf("Warning! %d header count error flag. Got %d header(s)!\n", header_count_Event_trailer);
        }
        bit1 = word >> 14; // get the TDC trailer count error flag in Event trailer
        error_flag_tmp = static_cast<unsigned int>((bit1.to_ulong()));
        if (error_flag_tmp)
        {
          printf("Warning! %d trailer count error flag. Got %d trailer(s)!\n", trailer_count_Event_trailer);
        }
        bit3 = word >> 36; // get the CSM ID
        int isCSM2 = static_cast<unsigned int>((bit3.to_ulong()));
        total_events++;
        total_signals += event_signals;
        if (trailer_sig.HitCount() != (sigVec.size() + error_word_count))
        {
          printf("Warning! Hit count in trailer = %d, real hit count = %d, error word count = %d\n",
                 trailer_sig.HitCount(), sigVec.size(), error_word_count);
        }

        // std::cout << "trailer_sig.HitCount() "<< trailer_sig.HitCount() << std::endl;
        if (trailer_sig.HitCount() != 0)
        {

          event = Event(header_sig, trailer_sig, sigVec);
          // DoHitFinding(&event,tc,ncut,0);
          DoHitFinding(&event, tc, 0, 0);
          // DoHitClustering(&event, geo);
          pass_event_check = kTRUE;
          // pass_event_check = CheckEvent(event, geo);
          event.SetPassCheck(pass_event_check);
          event.CheckClusterTime();
          for (Hit h : event.Hits()){
	          p_tdc_chnl_adc_time_raw          [h.TDC()][h.Channel()]->Fill(h.ADCTime()); 
          }

          if (pass_event_check)
          {
            eTree->Fill();

            // for (Cluster c : event.Clusters()) {
            for (Hit h : event.Hits())
            {

              // for (Hit h : c.Hits()) {
              p_tdc_tdc_time_corrected[h.TDC()]->Fill(h.CorrTime());
              p_tdc_adc_time[h.TDC()]->Fill(h.ADCTime());

              p_tdc_chnl_tdc_time_corrected[h.TDC()][h.Channel()]->Fill(h.CorrTime());
              p_tdc_chnl_adc_time[h.TDC()][h.Channel()]->Fill(h.ADCTime());
              // p_tdc_time_corrected[h.TDC()][h.Channel()]->Fill(h.CorrTime());
              // p_tdc_time          [h.TDC()][h.Channel()]->Fill(h.DriftTime());

              // std::cout << "TDC "<< h.TDC() << " Channel " << h.Channel() << " CorrTime " << h.CorrTime() << " DriftTime " << h.DriftTime() << std::endl;

              // 	      if(h.TDC() == 2){
              // 		      std::cout << " TDC 2 Channel " << h.Channel() << " CorrTime " << h.CorrTime() << " DriftTime " << h.DriftTime() << std::endl;
              // 	      }
              // 	      if(h.TDC() == 3){
              // 		      std::cout << "TDC 3 Channel " << h.Channel() << " CorrTime " << h.CorrTime() << " DriftTime " << h.DriftTime() << std::endl;
              // 	      }

              //    p_adc_time          [h.TDC()][h.Channel()]->Fill(h.ADCTime());

              geo.GetHitLayerColumn(h.TDC(), h.Channel(), &hitL, &hitC);
              hitByLC->Fill(hitC, hitL);
              if (h.CorrTime() < 0 || h.CorrTime() > 400)
                badHitByLC->Fill(hitC, hitL);
              else
                goodHitByLC->Fill(hitC, hitL);
              // p_hits_distribution[hitL]->Fill(hitC);
            } // end for (Hit h : event.Hits())
            // } //end for (Cluster c : event.Clusters())

            TTree *optTree = new TTree("optTree", "optTree");
            optTree->Branch("event", "Event", &event);
            optTree->Fill();

            // TOUPDATE timeprocess fill
            tp.setTarget(optTree);
            tp.setRangeSingle(0);
            tp.setIgnoreNone();
            tp.optimize();

            for (Cluster c : event.Clusters())
            {
              for (Hit h : c.Hits())
              {
                residuals->Fill(tp.Residual(h) * 1000.0);
              }
            }
            ed_event = event;

            if (pass_event_check)
            {
              total_events_pass++;
            }
            else
            {
              total_events_fail++;
            }

            // fill efficiency distribution
            double _hitX, _hitY;
            for (int iL = 0; iL < Geometry::MAX_TUBE_LAYER; iL++)
            {
              for (int iC = 0; iC < Geometry::MAX_TUBE_COLUMN; iC++)
              {
                if (geo.IsActiveTDCChannel(iL, iC))
                {
                  geo.GetHitXY(iL, iC, &_hitX, &_hitY);
                  // get track x position and figure out what tube(s) it may go through
                  double trackDist = tp.Distance(Hit(0, 0, 0, 0, iL, iC));
                  if (trackDist <= Geometry::column_distance / 2)
                  {
                    Bool_t tubeIsHit = kFALSE;
                    for (Hit h : event.Hits())
                    {
                      int *hit_layer;
                      int *hit_column;
                      geo.GetHitLayerColumn(h.TDC(), h.Channel(), hit_layer, hit_column);
                      if (*hit_layer == iL && *hit_column == iC)
                        tubeIsHit = kTRUE;
                    }

                    // TODO: THIS IS A HORRIBLE HACK
                    int col = iC + isCSM2 * Geometry::MAX_TUBE_COLUMN / 2;
                    if (col >= Geometry::MAX_TUBE_COLUMN)
                    {
                      continue;
                    }

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
          } // end if: pass_event_check

          /*
            // plot the event
            // plot the first 10 events meeting and not meeting the pass event check criteria
          if ((pass_event_check && total_events_pass < 10) || (!pass_event_check && total_events_fail < 10)) {

            if (pass_event_check)
              sprintf(track_group_name, "events_passing");
            else
              sprintf(track_group_name, "events_failing");

            if ((pass_event_check && total_events_pass == 1) || (!pass_event_check&& total_events_fail == 1) )
              event_track[(int)pass_event_check] = p_output_rootfile->mkdir(track_group_name);


            #ifdef SAVE_TRACKS_OUT_OF_ROOT
            if (mkdir(track_group_name, 0777) == -1) {
              cerr << strerror(errno) << endl;
            }
            #endif

            event_track[(int)pass_event_check]->cd();
            // chdir(track_group_name);
            // ed.DrawEvent(event, geo, event_track[(int)pass_event_check]);
            // chdir("..");
            // ed.Clear();
          } // end if: pass event check for first 100 events
          */

          for (int i = 0; i != Geometry::MAX_TDC; i++)
          {
            for (int j = 0; j != Geometry::MAX_TDC_CHANNEL; j++)
            {
              first_signal_flag[i][j] = kFALSE;
            }
          }
          if (pass_event_check)
            sprintf(track_group_name, "events_passing");
          else
            sprintf(track_group_name, "events_failing");

          if ((pass_event_check && total_events_pass == 1) || (!pass_event_check && total_events_fail == 1))
          {
            event_track[(int)pass_event_check] = p_output_rootfile->mkdir(track_group_name);
          }

          // Plot every 100th event
          if (pass_event_check && (total_events_pass % (10 * SPEEDFACTOR) == 0))
          {
            // printf("Entered event display loop\n");

            printf("(int) pass_event_check = %i\n", (int)pass_event_check);
            event_track[(int)pass_event_check]->cd();

            // ed->DrawTubeHistAndEvent(event, geo, goodHitByLC);
            gSystem->ProcessEvents();

            ed->Clear();
          }
          //
        }
      } // end if: got event trailer

      // got event header
      else if (header_type == Signal::HEADER)
      {
        header_sig = Signal(word);
        current_event_ID = header_sig.HeaderEID();
        if ((current_event_ID != (previous_event_ID + 1) % 4096) && (previous_event_ID != -1))
          printf("Warning! Event lost! Curr=%d, Pre=%d\n", current_event_ID, previous_event_ID);
        if (current_event_ID == previous_event_ID)
          printf("Warning! Repeating event! Curr=%d, Pre=%d\n", current_event_ID, previous_event_ID);
        previous_event_ID = current_event_ID;
        sigVec.clear();
        event_signals = 0;
        error_word_count = 0;
        tdc_header_eventID_current = 0;
        tdc_header_eventID_previous = -1; // reset tdc header eventID for each event
        header_count_TDC_data = 0;
        trailer_count_TDC_data = 0;
      }

      // got hit data
      else
      {
        sig = Signal(word);
        if (geo.IsActiveTDC(sig.TDC()) && sig.Channel() < Geometry::MAX_TDC_CHANNEL)
        {
          // printf("data from TDC %d, Channel = %d\n",sig.TDC(),sig.Channel());
          sigVec.push_back(sig);
        }
        else if (geo.IsActiveTDC(sig.TDC()) && sig.Channel() == ERROR_WORD)
        {                   // TDC triggerless error word
          bit3 = word >> 5; // get the LSB error flag
          error_flag = static_cast<unsigned int>((bit3.to_ulong()));
          error_word_count++;
          if (error_flag > 0)
          {
            bit5 = word >> 0; // get the LSB error channel
            error_channel = static_cast<unsigned int>((bit5.to_ulong()));
            printf("Warning! TDCID = %d, 	Channel = %d 	overflowed!\n", sig.TDC(), error_channel);
          }
          bit3 = word >> 13; // get the second LSB error flag
          error_flag = static_cast<unsigned int>((bit3.to_ulong()));
          if (error_flag > 0)
          {
            bit5 = word >> 8; // get the second LSB error channel
            error_channel = static_cast<unsigned int>((bit5.to_ulong()));
            printf("Warning! TDCID = %d, 	Channel = %d 	overflowed!\n", sig.TDC(), error_channel);
          }
        }
        else if (geo.IsActiveTDC(sig.TDC()) && sig.Channel() == HEADER_WORD)
        {
          header_count_TDC_data++;
          tdc_header_eventID_current = sig.TDCHeaderEID();
          if ((tdc_header_eventID_current != tdc_header_eventID_previous) && (tdc_header_eventID_previous != -1))
          {
            printf("TDC %d EventID mismatch! Current = %d, previous = %d\n", sig.TDC(), tdc_header_eventID_current, tdc_header_eventID_previous);
          }
          tdc_header_eventID_previous = tdc_header_eventID_current;
        }
        else if (geo.IsActiveTDC(sig.TDC()) && sig.Channel() == TRAILER_WORD)
        {
          trailer_count_TDC_data++;
        }
        else
        {
          printf("Error! At event %d unexpected data TDCID = %d, Channel=%d\n", total_events, sig.TDC(), sig.Channel());
          unexpected_data++;
          if (unexpected_data > 10)
            return (1);
        }
        event_signals++;
        // if(event_signals>1000)return 0;
      } // end got hit data
    }   // end while dataloop

    cout << "Decoding completed !" << endl;
    if (total_events <= 1)
    {
      cout << "No event! Terminated." << endl;
      return 1;
    }

    for (Int_t tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
    {
      if (geo.IsActiveTDC(tdc_id))
      {
        for (int tdc_chnl_id = 0; tdc_chnl_id != Geometry::MAX_TDC_CHANNEL; tdc_chnl_id++)
        {
          p_tdc_hit_rate[tdc_id][tdc_chnl_id] =
              p_tdc_chnl_adc_time_raw[tdc_id][tdc_chnl_id]->GetEntries() / match_window * 1000 / total_events;
        }
      }
    }

    if (iter % SPEEDFACTOR == 0)
    {
      float max = 0.0;
      for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; ++tdc_id)
      {
        isDrawn[tdc_id] = 0;
        max = 0.0;
        for (int iBin = 1; iBin <= p_tdc_hit_rate_graph[tdc_id]->GetNbinsX(); ++iBin)
        {
          //cout<<tdc_id<<" TestiBin:"<<iBin<< p_tdc_hit_rate[tdc_id][iBin - 1];
          p_tdc_hit_rate_graph[tdc_id]->SetBinContent(iBin, p_tdc_hit_rate[tdc_id][iBin - 1]);
          if (p_tdc_hit_rate[tdc_id][iBin - 1] > max)
            max = p_tdc_hit_rate[tdc_id][iBin - 1];
          p_tdc_hit_rate_graph[tdc_id]->SetMaximum(max);
        }
      }

      for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
      {
        string text_content;
        int CSM = tdc_id >= 18;
        int local_tdc_id = tdc_id % 18;

        if (geo.IsActiveTDC(tdc_id) || tdc_id == geo.TRIGGER_MEZZ)
        {

          if (tdc_id == geo.TRIGGER_MEZZ)
          {

            trigger_rate_canvas->cd();
            text_content = "Entries = " + to_string((int)total_triggers);
          }
          else
          {
            rate_canvas->cd(6 * ((local_tdc_id + 1) % 2) + (local_tdc_id / 2) + 1);
            text_content = "Entries = " + to_string((int)p_tdc_adc_time[tdc_id]->GetEntries());
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

          if (tdc_id != geo.TRIGGER_MEZZ)
          {
            if (CSM && isDrawn[local_tdc_id])
              opts = "same";
            else
              opts = "";
            adc_canvas->cd(6 * ((local_tdc_id + 1) % 2) + (local_tdc_id / 2) + 1);
            // cout<<"Drawing"<<6*((local_tdc_id+1)%2)+(local_tdc_id/2)+1<<":"<<tdc_id<<endl;

            p_tdc_adc_time[tdc_id]->Draw(opts);
            tdc_canvas->cd(6 * ((local_tdc_id + 1) % 2) + (local_tdc_id / 2) + 1);
            p_tdc_tdc_time_corrected[tdc_id]->Draw(opts);
          }

          isDrawn[local_tdc_id] = 1;
        }
      }

      for (int i = 1; i != pad_num + 1; i++)
      {

        adc_canvas->cd(i);
        gPad->Modified();
        tdc_canvas->cd(i);
        gPad->Modified();
        rate_canvas->cd(i);
        gPad->Modified();
      }

      // Update plots
      adc_canvas->cd();
      adc_canvas->Modified();
      adc_canvas->Update();
      adc_canvas->Print("/atlas/data18a/rongqian/test1.png");

      tdc_canvas->cd();
      tdc_canvas->Modified();
      tdc_canvas->Update();
      tdc_canvas->Print("/atlas/data18a/rongqian/test2.png");

      rate_canvas->cd();
      rate_canvas->Update();
      trigger_rate_canvas->cd();
      trigger_rate_canvas->Update();

      residual_canvas->cd();
      residuals->Draw();
      residual_canvas->Update();
      EDCanvas->cd();
      ed_event.AddTrack(Track(tp.slope(), tp.y_int()));
      // ed->DrawEvent(ed_event, geo, NULL);
      gSystem->ProcessEvents();
      printf("Plots got updated.\n");
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
      eff_canvas->Modified();
      eff_canvas->Update();
      struct Channel_packet p_chnl;

      for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
      {
        if (geo.IsActiveTDC(tdc_id))
        {
          if (tdc_id == geo.TRIGGER_MEZZ)
          {
            // write p_chnl information
            p_chnl.tdc_id = geo.TRIGGER_MEZZ;
            p_chnl.tdc_chnl_id = geo.TRIGGER_CH;

            // write ADC histogram
            p_chnl.adc_entries_raw = 0;
            for (int bin_index = 0; bin_index < ADC_HIST_TOTAL_BIN; bin_index++)
            {
              p_chnl.adc_hist[bin_index] = 0;
            }
            // write TDC histogram
            p_chnl.adc_entries_raw = 0;
            for (int bin_index = 0; bin_index < TDC_HIST_TOTAL_BIN; bin_index++)
            {
              p_chnl.tdc_hist[bin_index] = 0;
            }
            // write ADC raw histogram
            p_chnl.adc_entries_raw = p_tdc_chnl_adc_time_raw[geo.TRIGGER_MEZZ][geo.TRIGGER_CH]->GetEntries();
            for (int bin_index = 0; bin_index < ADC_HIST_TOTAL_BIN; bin_index++)
            {
              p_chnl.adc_hist_raw[bin_index] = p_tdc_chnl_adc_time_raw[geo.TRIGGER_MEZZ][geo.TRIGGER_CH]->GetBinContent(bin_index + 1);
            }
            // write TDC raw histogram
            p_chnl.tdc_entries_raw = p_tdc_chnl_adc_time_raw[geo.TRIGGER_MEZZ][geo.TRIGGER_CH]->GetEntries();
            for (int bin_index = 0; bin_index < TDC_HIST_TOTAL_BIN; bin_index++)
            {
              p_chnl.tdc_hist_raw[bin_index] = p_tdc_chnl_tdc_time_corrected_raw[geo.TRIGGER_MEZZ][geo.TRIGGER_CH]->GetBinContent(bin_index + 1);
            }
            // sending UDP packet outs
            sendto(udp_sock_fd, (char *)&p_chnl, sizeof(p_chnl),
                   MSG_CONFIRM, (const struct sockaddr *)&udp_servaddr, sizeof(udp_servaddr));
          }

          else
          {
            for (int tdc_chnl_id = 0; tdc_chnl_id != Geometry::MAX_TDC_CHANNEL; tdc_chnl_id++)
            {
              // write p_chnl information
              p_chnl.tdc_id = tdc_id;
              p_chnl.tdc_chnl_id = tdc_chnl_id;

              // write ADC histogram
              p_chnl.adc_entries = p_tdc_chnl_adc_time[tdc_id][tdc_chnl_id]->GetEntries();
              for (int bin_index = 0; bin_index < ADC_HIST_TOTAL_BIN; bin_index++)
              {
                p_chnl.adc_hist[bin_index] = p_tdc_chnl_adc_time[tdc_id][tdc_chnl_id]->GetBinContent(bin_index + 1);
              }

              // write TDC histogram
              p_chnl.tdc_entries = p_tdc_chnl_tdc_time_corrected[tdc_id][tdc_chnl_id]->GetEntries();
              for (int bin_index = 0; bin_index < TDC_HIST_TOTAL_BIN; bin_index++)
              {
                p_chnl.tdc_hist[bin_index] = p_tdc_chnl_tdc_time_corrected[tdc_id][tdc_chnl_id]->GetBinContent(bin_index + 1);
              }

              // write ADC raw histogram
              p_chnl.adc_entries_raw = p_tdc_chnl_adc_time_raw[tdc_id][tdc_chnl_id]->GetEntries();
              for (int bin_index = 0; bin_index < ADC_HIST_TOTAL_BIN; bin_index++)
              {
                p_chnl.adc_hist_raw[bin_index] = p_tdc_chnl_adc_time_raw[tdc_id][tdc_chnl_id]->GetBinContent(bin_index + 1);
              }

              // write TDC raw histogram
              p_chnl.tdc_entries_raw = p_tdc_chnl_tdc_time_corrected_raw[tdc_id][tdc_chnl_id]->GetEntries();
              for (int bin_index = 0; bin_index < TDC_HIST_TOTAL_BIN; bin_index++)
              {
                p_chnl.tdc_hist_raw[bin_index] = p_tdc_chnl_tdc_time_corrected_raw[tdc_id][tdc_chnl_id]->GetBinContent(bin_index + 1);
              }

              // sending UDP packet out
              //  printf("TDC=%d, CHNL=%d, Entries=%d\n",p_chnl.tdc_id,p_chnl.tdc_chnl_id,p_chnl.adc_entries);
              sendto(udp_sock_fd, (char *)&p_chnl, sizeof(p_chnl),
                     MSG_CONFIRM, (const struct sockaddr *)&udp_servaddr, sizeof(udp_servaddr));
            }
          }
        }
      }

      //  if(bytes_recv<=0)break;
    } // end if : iter%SPEEDFACTOR==0
    if (gSystem->ProcessEvents())
    {
      std::cout << "Processing Events" << std::endl;
      break;
    }
    cout << "Test clear and close" << endl;
    data_in_flow.clear();
    data_in_flow.close();

  } // end while 1

  p_output_rootfile->cd();
  for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++) {
    if (geo.IsActiveTDC(tdc_id)) {
      if (tdc_id == geo.TRIGGER_MEZZ) continue;
      p_tdc_adc_time[tdc_id]->Write();
      p_tdc_tdc_time_corrected[tdc_id]->Write();
    }
  }
  
  p_output_rootfile->Write();
  eTree->Write();
  int nEntries = eTree->GetEntries();
  delete p_output_rootfile;
  
  
  oFile.close();
  close(newsockfd);
  close(sockfd);
  data_in_flow.close();

  // create output file

  system("mkdir output");
  chdir("output");
  char output_directoryname[256];
  memset(output_directoryname, 0, sizeof(output_directoryname));
  sprintf(output_directoryname, "mkdir %s", filename_time);
  system(output_directoryname);
  chdir(filename_time);

  char rate_canvas_name[256];
  memset(rate_canvas_name, 0, sizeof(rate_canvas_name));
  sprintf(rate_canvas_name, "%s_rate.png", filename_time);
  rate_canvas->Print(rate_canvas_name);

  char adc_canvas_name[256];
  memset(adc_canvas_name, 0, sizeof(adc_canvas_name));
  sprintf(adc_canvas_name, "%s_adc.png", filename_time);
  adc_canvas->Print(adc_canvas_name);

  char tdc_canvas_name[256];
  memset(tdc_canvas_name, 0, sizeof(tdc_canvas_name));
  sprintf(tdc_canvas_name, "%s_tdc.png", filename_time);
  tdc_canvas->Print(tdc_canvas_name);

  char trigger_canvas_name[256];
  memset(trigger_canvas_name, 0, sizeof(trigger_canvas_name));
  sprintf(trigger_canvas_name, "%s_trigger_rate.png", filename_time);
  trigger_rate_canvas->Print(trigger_canvas_name);

  char fp_rate_File_name[256];
  memset(fp_rate_File_name, 0, sizeof(fp_rate_File_name));
  sprintf(fp_rate_File_name, "%s_rate.csv", filename_time);

  fp_rate_File = fopen(fp_rate_File_name, "w");
  fprintf(fp_rate_File, "tdc_id,");
  for (int tdc_chnl_id = 0; tdc_chnl_id != Geometry::MAX_TDC_CHANNEL; tdc_chnl_id++)
  {
    fprintf(fp_rate_File, "%d,", tdc_chnl_id);
  }

  for (int tdc_id = 0; tdc_id != Geometry::MAX_TDC; tdc_id++)
  {
    if (geo.IsActiveTDC(tdc_id))
    {
      fprintf(fp_rate_File, "\n%d,", tdc_id);
      for (int tdc_chnl_id = 0; tdc_chnl_id != Geometry::MAX_TDC_CHANNEL; tdc_chnl_id++)
      {
        fprintf(fp_rate_File, "%.4f,", p_tdc_hit_rate[tdc_id][tdc_chnl_id]);
      }
    }
  }
  fprintf(fp_rate_File, "%d,", total_events);
  fclose(fp_rate_File);

  cout << endl;
  cout << "Total Triggers: " << total_triggers << endl;
  cout << "Pass  Triggers: " << total_triggers_pass << endl;
  cout << endl;
  cout << "Total Events:   " << total_events << endl;
  cout << "Pass  Events:   " << total_events_pass << endl;
  cout << endl;
  cout << "Total Signals:  " << total_signals << endl;
  cout << "Pass  Signals:  " << total_signals_pass << endl;
  cout << endl;
  cout << "N tree entries: " << nEntries << endl;
  adc_canvas->Print("/atlas/data18a/rongqian/test3.png");

  // gROOT->SetBatch(kFALSE);
  return 0;
}

#endif

int OnlineDAQ_fake(short portno, int bisno = 1)
{
  DAQ_monitor *p_DAQ_monitor = new DAQ_monitor(portno, bisno);
  TString filename = "run_20231003_152805.dat";
  p_DAQ_monitor->DecodeOnline(filename);
  printf("END");
  return 1;

}
