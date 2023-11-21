/**
 * @file PlotSavingOperations.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <thread>

#include "macros/DAQState.cpp"

#include "src/Geometry.cpp"

#include "src/DataModel/DAQData.cpp"

#include "DAQMonitor/PlotSaving/Views/ProgressBar.cpp"

using namespace std;

// TODO: Remove dependency on the progress bar

namespace Muon {
namespace PlotSaving {

	void savePlots();

}
namespace PlotSavingIMPL {

    void makeDirectory(const string &path);

}
}

void Muon::PlotSaving::savePlots() {

    using namespace Muon::PlotSavingIMPL;

    cout << "Saving snapshot..." << endl;

    makeDirectory("output");
    cout << "Created output directory." << endl;

    State::DAQState state = State::DAQState::getState();
    string outputDirName = string("output/") + state.tempState.runLabel;
    makeDirectory(outputDirName);
    cout << "Created directory " << outputDirName << endl;

    // TODO: We need to save the root output too -- see DecodeOffline.cpp:253

    // TODO: Figure out behavior for saving the same run twice
    //         -- I think it should overwrite the previous save

    // Plot saving should happen in a separate thread, and we should have a progress bar
    // associated with it
    //   -- Maybe it's okay to block the UI while saving plots though, as long as we 
    //      provide feedback that makes it clear that the UI is blocked.
    //   -- NOTE: It's a little tricky to do plot saving outside the UI thread -- we need
    //            a mutex or the output canvas overwrites other plots

    DAQData &data = DAQData::getInstance();

    data.lock();
    Plots snapshot(data.plots);
    data.unlock();

    // TODO: We need to:
    //         -- Set up the directory structure
    //         -- Save all the plots

    // NOTE: We can speed things up a lot by setting gROOT->SetBatch(), but
    //       this also breaks all the graphics. Useful if saves only happen
    //       on stopRun, but not otherwise. Setting outputCanvas->SetBatch()
    //       has no effect on performance.

    // TODO: Can we consider multithreading the save operation to speed things
    //       up?

    TRootEmbeddedCanvas *outputCanvas = new TRootEmbeddedCanvas("Output Canvas", gClient->GetRoot());

    // NOTE: We avoid displaying the window by omitting the calls to 
    //       MapSubwindows() and MapWindow()

    outputCanvas->SetWindowName("Output Preview");
    outputCanvas->Resize(700, 700);
    outputCanvas->GetCanvas()->cd();

    ProgressBar *progressBar = new ProgressBar(gClient->GetRoot());

    progressBar->SetWindowName("Save Progress");
    progressBar->MapSubwindows();
    progressBar->Resize(progressBar->GetDefaultSize());
    progressBar->MapWindow();

    int activeTDCs = 0;
    int activeChannels = 0;

    for(int tdc = 0; tdc < Geometry::MAX_TDC; ++tdc) {

        if(Geometry::getInstance().IsActiveTDC(tdc)) {

            ++activeTDCs;

        }

        for(int chnl = 0; chnl < Geometry::MAX_TDC_CHANNEL; ++chnl) {

            if(
                Geometry::getInstance().IsActiveTDCChannel(tdc, chnl) || 
                (tdc == Geometry::getInstance().TRIGGER_MEZZ)
            ) {

                ++activeChannels;

            }

        }

    }

    float incr = 100. / (2 * activeTDCs + 3 * activeChannels);

    for(int tdc = 0; tdc < Geometry::MAX_TDC; ++tdc) {

        if(Geometry::getInstance().IsActiveTDC(tdc)) {

            string dirName = outputDirName
                + "/TDC_" 
                + to_string(tdc) 
                + "_of_" 
                + to_string(Geometry::MAX_TDC) 
                + "_Time_Spectrum";

            makeDirectory(dirName);

            cout << "Created directory " << dirName << endl;

            // TODO: Consider this implementation:
            //       https://root.cern/doc/v608/pad2png_8C.html

            string filename = "";

            snapshot.p_tdc_tdc_time_corrected[tdc]->Draw("colz");
            filename = string("tdc_") + tdc + "_tdc_time_spectrum_corrected.png";
            outputCanvas->GetCanvas()->SaveAs((dirName + "/" + filename).data());

            progressBar->increment(incr);

            snapshot.p_tdc_adc_time[tdc]->Draw("colz");
            filename = string("tdc_") + tdc + "_adc_time_spectrum.png";
            outputCanvas->GetCanvas()->SaveAs((dirName + "/" + filename).data());

            progressBar->increment(incr);

            for(int chnl = 0; chnl < Geometry::MAX_TDC_CHANNEL; ++chnl) {

                if(
                    Geometry::getInstance().IsActiveTDCChannel(tdc, chnl) || 
                    (tdc == Geometry::getInstance().TRIGGER_MEZZ)
                ) {

                    snapshot.p_tdc_time_corrected[tdc][chnl]->Draw("colz");
                    filename = string("tdc_") + tdc + "__channel_" + chnl + "__tdc_time_spectrum_corrected.png";
                    outputCanvas->GetCanvas()->SaveAs((dirName + "/" + filename).data());

                    progressBar->increment(incr);

                    snapshot.p_tdc_time[tdc][chnl]->Draw("colz");
                    filename = string("tdc_") + tdc + "__channel_" + chnl + "__tdc_time_spectrum.png";
                    outputCanvas->GetCanvas()->SaveAs((dirName + "/" + filename).data());

                    progressBar->increment(incr);

                    snapshot.p_adc_time[tdc][chnl]->Draw("colz");
                    filename = string("tdc_") + tdc + "__channel_" + chnl + "__adc_time_spectrum.png";
                    outputCanvas->GetCanvas()->SaveAs((dirName + "/" + filename).data());

                    progressBar->increment(incr);

                }

            }

        }

    }

    delete outputCanvas;
    delete progressBar;

    cout << "Plots saved!" << endl;

}

void Muon::PlotSavingIMPL::makeDirectory(const string &path) {

  if(mkdir(path.data(), 0777) == -1) {
    cerr << strerror(errno) << endl;
  }

}