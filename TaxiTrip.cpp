#include <bits/stdc++.h>
using namespace std;

struct TaxiTrip {
    string medallion;
    string hackLicense;
    string vendorID;
    int rateCode;
    bool storeAndFwdFlag;
    string pickUpDateTime;
    string dropOffDateTime;
    int passengerCount;
    int tripTimeSecs;
    double tripDistance;
    pair<double, double> pickupLocation;  // latitud, longitud
    pair<double, double> dropoffLocation; // latitud, longitud

    TaxiTrip(string med, string hack, string vendor, int rate, bool storeFlag, 
    string pickupDT, string dropoffDT, int passengerCount, int tripTime, 
    double tripDist, double pickupLon, double pickupLat, double dropoffLon, double dropoffLat) 
    : medallion(med), hackLicense(hack), vendorID(vendor), rateCode(rate), storeAndFwdFlag(storeFlag), 
    pickUpDateTime(pickupDT), dropOffDateTime(dropoffDT), passengerCount(passengerCount), 
    tripTimeSecs(tripTime), tripDistance(tripDist), 
    pickupLocation(pickupLat, pickupLon), dropoffLocation(dropoffLat, dropoffLon) {}

    void displayInfo() const {
        cout << "Taxi Trip Info: \n";
        cout << "Medallion: " << medallion << endl;
        cout << "Hack License: " << hackLicense << endl;
        cout << "Vendor ID: " << vendorID << endl;
        cout << "Rate Code: " << rateCode << endl;
        cout << "Store and Forward Flag: " << (storeAndFwdFlag ? "Yes" : "No") << endl;
        cout << "Pickup Datetime: " << pickUpDateTime << endl;
        cout << "Dropoff Datetime: " << dropOffDateTime << endl;
        cout << "Passenger Count: " << passengerCount << endl;
        cout << "Trip Time (secs): " << tripTimeSecs << endl;
        cout << "Trip Distance (miles): " << tripDistance << endl;
        cout << "Pickup Location: (" << pickupLocation.first << ", " << pickupLocation.second << ")\n";
        cout << "Dropoff Location: (" << dropoffLocation.first << ", " << dropoffLocation.second << ")\n";
        cout << "--------------------------------------\n";
    }
};