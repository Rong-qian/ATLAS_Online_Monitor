/**
 * @file RunStats.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <string>

using namespace std;

// TODO: Examine this example:
//       https://root.cern/doc/master/calendar_8C.html

class RunStats : public TGVerticalFrame {

public:

	RunStats(
		const TGWindow *parent
	);

private:

	// VIEW

	TGHorizontalFrame *packetCounter;
		TGLabel *packetCountLabel;
		TGLabel *packetCount;

	TGHorizontalFrame *signalCounter;
		TGLabel *signalCountLabel;
		TGLabel *signalCount;

	TGHorizontalFrame *eventCounter;
		TGLabel *eventCountLabel;
		TGLabel *eventCount;

	TGHorizontalFrame *nonemptyEventCounter;
		TGLabel *nonemptyEventCountLabel;
		TGLabel *nonemptyEventCount;

	TGHorizontalFrame *errorCounter;
		TGLabel *errorCountLabel;
		TGLabel *errorCount;

	TGHtml *table;

};

RunStats::RunStats(
	const TGWindow *parent
) : TGVerticalFrame(parent) {

	/*
	packetCounter = new TGHorizontalFrame(this);
	AddFrame(packetCounter, new TGLayoutHints(kLHintsCenterX));

		packetCountLabel = new TGLabel(packetCounter, "Packets Recorded:  ");
		packetCounter->AddFrame(packetCountLabel, new TGLayoutHints(kLHintsLeft));

		packetCount = new TGLabel(packetCounter, "UNIMPLEMENTED");
		packetCounter->AddFrame(packetCount, new TGLayoutHints(kLHintsRight));

	signalCounter = new TGHorizontalFrame(this);
	AddFrame(signalCounter, new TGLayoutHints(kLHintsCenterX));

		signalCountLabel = new TGLabel(signalCounter, "Processed Signals:  ");
		signalCounter->AddFrame(signalCountLabel, new TGLayoutHints(kLHintsLeft));

		signalCount = new TGLabel(signalCounter, "UNIMPLEMENTED");
		signalCounter->AddFrame(signalCount, new TGLayoutHints(kLHintsRight));

	eventCounter = new TGHorizontalFrame(this);
	AddFrame(eventCounter, new TGLayoutHints(kLHintsCenterX));

		eventCountLabel = new TGLabel(eventCounter, "Processed Events:  ");
		eventCounter->AddFrame(eventCountLabel, new TGLayoutHints(kLHintsLeft));

		eventCount = new TGLabel(eventCounter, "UNIMPLEMENTED");
		eventCounter->AddFrame(eventCount, new TGLayoutHints(kLHintsRight));

	nonemptyEventCounter = new TGHorizontalFrame(this);
	AddFrame(nonemptyEventCounter, new TGLayoutHints(kLHintsCenterX));

		nonemptyEventCountLabel = new TGLabel(nonemptyEventCounter, "Processed Nonempty Events:  ");
		nonemptyEventCounter->AddFrame(nonemptyEventCountLabel, new TGLayoutHints(kLHintsLeft));

		nonemptyEventCount = new TGLabel(nonemptyEventCounter, "UNIMPLEMENTED");
		nonemptyEventCounter->AddFrame(nonemptyEventCount, new TGLayoutHints(kLHintsRight));

	errorCounter = new TGHorizontalFrame(this);
	AddFrame(errorCounter, new TGLayoutHints(kLHintsCenterX));

		errorCountLabel = new TGLabel(errorCounter, "Errors:  ");
		errorCounter->AddFrame(errorCountLabel, new TGLayoutHints(kLHintsLeft));

		errorCount = new TGLabel(errorCounter, "UNIMPLEMENTED");
		errorCounter->AddFrame(errorCount, new TGLayoutHints(kLHintsRight));
	*/
	table = new TGHtml(this, 100, 100);
	AddFrame(table, new TGLayoutHints(/*kLHintsExpandX | kLHintsExpandY*/kLHintsCenterX));

	table->Clear();

	char htmlText[256] = "<TABLE width=100%><TBODY><TR><TH>One</TH><TH>Two</TH></TR></TBODY></TABLE>";
	table->ParseText(htmlText);

	table->Layout();

}