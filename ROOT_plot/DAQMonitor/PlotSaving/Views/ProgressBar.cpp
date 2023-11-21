/**
 * @file ProgressBar.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

using namespace std;

// TODO: Disable the close window button

class ProgressBar : public TGTransientFrame {

public:

	ProgressBar(const TGWindow *p);

	void increment(float inc);

private:

	// VIEW

	TGHProgressBar *bar;

};

ProgressBar::ProgressBar(const TGWindow *p) : TGTransientFrame(p) {

	bar = new TGHProgressBar(this, TGProgressBar::kFancy, 300);
	AddFrame(bar, new TGLayoutHints(kLHintsExpandX));

	bar->SetBarColor("green");

	bar->ShowPosition();

}

void ProgressBar::increment(float inc) {

	bar->Increment(inc);

}