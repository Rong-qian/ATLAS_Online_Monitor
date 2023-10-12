/**
 * @file DeviceSelectorPanel.cpp
 *
 * @brief TODO: Write
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include <string>

#include "src/EthernetCapture/DeviceManager.cpp"

using namespace std;

// TODO: Return a whole PCapDevice rather than just its name

class DeviceSelector : public TGVerticalFrame {

public:

    DeviceSelector(const TGWindow *p);

    virtual ~DeviceSelector() override;

    // METHODS

    void initialize();

    bool ready();

    void enable();
    void disable();

    string getDeviceName();
    void setDeviceName(char *name);

private:

    // VIEW

    ///////////////////////////////////////////////////////////////////////////
    TGLabel    *label   ;
    TGComboBox *dropdown;
    ///////////////////////////////////////////////////////////////////////////

    bool initialized;

};

string DeviceSelector::getDeviceName() {

    if(!dropdown->GetSelectedEntry()) return "";
    if(dropdown->GetSelected() == 0) return "";

    return string(dropdown->GetSelectedEntry()->GetTitle());

}

void DeviceSelector::setDeviceName(char *name) {

    // FIXME: The current entry must be cleared if name is blank
    if(!strcmp(name, "")) {

        dropdown->Select(0);
        return;

    }

    int entryID = dropdown->FindEntry(name)->EntryId();

    if(!entryID) {

        // We can handle this quietly. If user selects a nonexistent device,
        // we've done something very wrong. But the program might try it
        // on settings state load if a previously used device has been 
        // disconnected. In that case, we really just want to clear the default
        // setting.
        dropdown->Select(0);
        return;

    }

    dropdown->Select(entryID, false);

}

DeviceSelector::DeviceSelector(const TGWindow *p) 
    : TGVerticalFrame(p), initialized(false) {

    label = new TGLabel(this, "Select Device:");
    AddFrame(label, new TGLayoutHints(kLHintsLeft));


    dropdown = new TGComboBox(this);
    AddFrame(dropdown, new TGLayoutHints(kLHintsLeft));

    dropdown->Resize(150, 20);

}

DeviceSelector::~DeviceSelector() {

    TGVerticalFrame::~TGVerticalFrame();

}

void DeviceSelector::initialize() {

    DeviceManager devices;
    devices.initialize();

    vector<PCapDevice> deviceList = devices.getAllDevices();

    dropdown->AddEntry("(None)", 0);

    for(int i = 0; i < deviceList.size(); ++i) {

        dropdown->AddEntry(deviceList[i].name(), i + 1);

    }

    initialized = true;

}

bool DeviceSelector::ready() {

    return initialized;

}


void DeviceSelector::enable() {

    dropdown->SetEnabled(true);

    label->Enable();

}

void DeviceSelector::disable() {

    // TODO: Selections in the box are not grayed out

    dropdown->SetEnabled(false);

    label->Disable();

}