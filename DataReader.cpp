#include <vector>
#include "TaxiTrip.cpp"
using namespace std;

class DataReader {
    public:
    static vector<TaxiTrip> readTaxiTrips(const string& file_name, int maxRows) {
        vector<TaxiTrip> taxiTrips;
        ifstream file(file_name);
        string line;
        int rowCount = 0;

        getline(file, line);

        while (getline(file, line)) {
            if (rowCount >= maxRows) break;
            rowCount++;

            istringstream ss(line);
            string medallion, hackLicense, vendorID, pickupDatetime, dropoffDatetime;
            int rateCode, passengerCount, tripTimeSecs;
            double tripDistance, pickupLongitude, pickupLatitude, dropoffLongitude, dropoffLatitude;
            bool storeAndFwdFlag;

            getline(ss, medallion, ',');
            getline(ss, hackLicense, ',');
            getline(ss, vendorID, ',');

            ss >> rateCode;
            ss.ignore();

            string storeFlagStr;
            getline(ss, storeFlagStr, ',');
            storeAndFwdFlag = (storeFlagStr == "Y");

            // Pickup Datetime
            getline(ss, pickupDatetime, ',');
            // Dropoff Datetime
            getline(ss, dropoffDatetime, ',');

            // Passenger Count
            ss >> passengerCount;
            ss.ignore();

            // Trip Time in Seconds
            ss >> tripTimeSecs;
            ss.ignore();

            // Trip Distance
            ss >> tripDistance;
            ss.ignore();

            // Coordenadas
            ss >> pickupLatitude;
            ss.ignore();
            ss >> pickupLongitude;
            ss.ignore();
            ss >> dropoffLatitude;
            ss.ignore();
            ss >> dropoffLongitude;

            taxiTrips.push_back(TaxiTrip(
                medallion, hackLicense, vendorID, rateCode, storeAndFwdFlag,
                pickupDatetime, dropoffDatetime, passengerCount,
                tripTimeSecs, tripDistance, pickupLongitude, pickupLatitude,
                dropoffLongitude, dropoffLatitude
            ));
        }
        return taxiTrips;
    }
};