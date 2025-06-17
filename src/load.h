#include "algo.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath> //std::log
#include <algorithm> //std::min, std::max
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <set>
#include <unordered_set>
#include <unordered_map> 

// Custom hash function for double keys
struct DoubleHash {
    size_t operator()(const double& key) const {
        // Convert the double to an integer for hashing
        return std::hash<double>{}(key);
    }
};

// Custom equality function for double keys
struct DoubleEqual {
    bool operator()(const double& lhs, const double& rhs) const {
        // Compare doubles within a tolerance epsilon
        double epsilon = 0.0001; // Define your desired precision
        return std::abs(lhs - rhs) < epsilon;
    }
};


class Load {
public:
	
	//Generic method for BSS with commong data format
	void Load_CSV_File(BSS_Data & data, std::string filename_csv, int start_hour, int end_hour);	
	
	//Generic method for BSS with commong data format
	void Load_BSS_Data(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour, std::string type, int max_stations_allowed);
	
	//EcoBici methods
	void Load_Ecobici(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour);
	
	//Toronto methods
	void Load_Toronto(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour);	
	
	//Madrid methods
	void Load_Madrid(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour);	
	
	//Bixi methods
	void Load_Bixi(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour);
	void Load_Bixi(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour, int max_stations_allowed);

	//CapitalBikeShare method
	void Load_CapitalBikeShare(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour, int max_stations_allowed);

	//Maps
	std::unordered_map<std::string, Station> station_name_map;
	std::unordered_map<double, Station, DoubleHash, DoubleEqual> station_lat_map;
	std::unordered_map<double, Station, DoubleHash, DoubleEqual> station_lon_map;
	std::map<std::string,int> bss_id_shortName_map;

};

