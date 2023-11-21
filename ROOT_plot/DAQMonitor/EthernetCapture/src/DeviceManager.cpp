/**
 * @file DeviceManager.cpp
 *
 * @brief Manages network devices.
 *
 * @author Robert Myers
 * Contact: romyers@umich.edu
 */

#pragma once

#include "DAQMonitor/EthernetCapture/src/NetworkDeviceException.cpp"
#include "DAQMonitor/EthernetCapture/src/PCapDevice.cpp"

using namespace std;



class DeviceManager {

public:

    DeviceManager();
    ~DeviceManager();

    DeviceManager (      DeviceManager &other) = delete;
    bool operator=(const DeviceManager &other) = delete;

    void initialize();
    PCapDevice getDevice(int deviceNumber);
    vector<PCapDevice> getAllDevices();

    void printDevices(ostream &out) const;

private:

    pcap_if_t *devices;

};



DeviceManager::DeviceManager() : devices(nullptr) {}

DeviceManager::~DeviceManager() {

    if(devices) pcap_freealldevs(devices);
    devices = nullptr;

}

void DeviceManager::initialize() {

    char errorBuffer[PCAP_ERRBUF_SIZE];

    if(pcap_findalldevs(&devices, errorBuffer) == -1) {

        if(devices) pcap_freealldevs(devices);
        devices = nullptr;

        throw NetworkDeviceException(
            string("Error in pcap_findalldevs_ex: ") + errorBuffer
        );

    }

    if(devices == nullptr) {

        throw NetworkDeviceException(
            "No network devices found. Did you run with root privileges?"
        );

    }

}

PCapDevice DeviceManager::getDevice(int deviceNumber) {

    if(!devices) {

        throw NetworkDeviceException(
            "DeviceManager::getDevice -- DeviceManager must be initialized before use."
        );

    }

    pcap_if_t *d = devices;

    for(int i = 0; i < deviceNumber - 1; ++i) {

        d = d->next;

        // NOTE: This check must come after d = d->next. Otherwise it will 
        //       fail when the user selects a device number
        //       one greater than the total number of devices, which may be
        //       a difficult edge case to find.
        if(!d) {

            throw NetworkDeviceException(
                "Selected device does not exist."
            );

        }
    }

    return PCapDevice(d);

}

vector<PCapDevice> DeviceManager::getAllDevices() {

    if(!devices) {

        throw logic_error(
            "DeviceManager::getDevice -- DeviceManager must be initialized before use."
        );

    }

    vector<PCapDevice> deviceList;

    pcap_if_t *d = devices;

    while(d) {

        deviceList.push_back(PCapDevice(d));
        d = d->next;

    }

    return deviceList;


}

void DeviceManager::printDevices(ostream &out) const {

    if(!devices) {

        throw logic_error(
            "DeviceManager::printDevices -- DeviceManager must be initialized before use."
        );

    }

    int i = 0;
    for(pcap_if_t *d = devices; d != NULL; d = d->next) {

        ++i;

        cout << i << ": " << d->name << "    ";
        if(d->description) {
            cout << d->description << endl;
        } else {
            cout << "(No description available)" << endl;
        }

    }

}