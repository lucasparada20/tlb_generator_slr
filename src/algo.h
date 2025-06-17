#ifndef ALGO_2_H
#define ALGO_2_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <unordered_map>
#include <map>
#include "RandomNumbers.h"

#define P_CHAR 'p'
#define Q_CHAR 'q'

// Define a custom hash function for std::pair<int, int>
struct PairHash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2>& p) const {
		auto h1 = std::hash<T1>{}(p.first);
		auto h2 = std::hash<T2>{}(p.second);

		// Combine the hash values of the two elements
		return h1 ^ h2;
	}
};

class Station{
public:
	int id;
	std::string bss_station_id;
    std::string name;
    double lat;
    double lon;
    int cap;
	void Show(){std::cout << "Station " << "id:" << id << " lat:" << lat <<" lon:" << lon << " cap:" << cap << " name:"<< name << std::endl;}
	
};


class TripData {
public:
    std::string start_date_str;
    std::string start_time_str;
    std::string end_date_str;
    std::string end_time_str;
    double duration;
    int start_id;
    int end_id;

    // Default constructor
    TripData() = default;

    // Parameterized constructor
    TripData(
        std::string start_date_str,
        std::string start_time_str,
        std::string end_date_str,
        std::string end_time_str,
        double duration,
        int start_id,
        int end_id
    ) : start_date_str(std::move(start_date_str)),
        start_time_str(std::move(start_time_str)),
        end_date_str(std::move(end_date_str)),
        end_time_str(std::move(end_time_str)),
        duration(duration),
        start_id(start_id),
        end_id(end_id) {}

    void Show() const {
		printf("Trip start:%s end:%s duration:%.1lf from:%d to:%d\n",(start_date_str+'_'+start_time_str).c_str(),(end_date_str+'_'+end_time_str).c_str(),duration,start_id,end_id);
	}
};

class TripDataAllInts
{
public:
	TripDataAllInts() : start_year(-1), start_month(-1), start_day(-1), start_hour(-1), start_minute(-1), start_second(-1), end_year(-1), end_month(-1), end_day(-1), end_hour(-1), end_minute(-1), end_second(-1), start_id(-1), end_id(-1), duration(-1.0) {}

	int8_t start_month, start_day, end_month, end_day, start_hour, end_hour, start_minute, end_minute, start_second, end_second;
	int16_t start_year, end_year;
	int16_t start_id, end_id;
	double duration;
	
    void Show() const
    {
        printf("Trip start:%d/%d/%d %02d:%02d:%02d end:%d/%d/%d %02d:%02d:%02d from:%d to:%d dur:%.1lf\n",
               start_year, start_month, start_day, start_hour, start_minute, start_second,
               end_year, end_month, end_day, end_hour, end_minute, end_second,
               start_id, end_id, duration);
    }
};


class Requests {
public:
    char type;
	int stat_id;
    int time;
    int scenario;
	
	Requests(char _type = 'n', int _stat_id=-1, int _time=-1, int _scenario=-1)
		: type(_type), stat_id(_stat_id), time(_time), scenario(_scenario){}	

    void Show() const {
        printf("req:%c station:%d scenario:%d\n", type, stat_id,scenario);
    }
};


class Trip {
public:
    int origin;
    int destination;
    int start_time;
    int end_time;
    int scenario;

    // Constructor to initialize Trip object
    Trip(int _origin = -1, int _destination = -1, int _start_time = -1, int _end_time = -1, int _scenario = -1)
        : origin(_origin), destination(_destination), start_time(_start_time), end_time(_end_time), scenario(_scenario) {}

    void Show() const {
        printf("o:%d d:%d start_time:%d end_time:%d scenario:%d\n", origin, destination, start_time, end_time, scenario);
    }
};

struct BSS_Data {
		// Copy constructor
    BSS_Data(const BSS_Data& other) {
        trips = other.trips;
        stations = other.stations;
    }
	
	//To count the data days in the load method ...
	int dayCounter;
	int nb_days_of_data;
	std::unordered_map<std::string, int> dateToDayMap;
	
	//Takes the trips in std::vector<TripData> and generates instance data
	void ConvertToTrips();
	//Takes the trips in std::vector<TripData> and generates instance data
	void ConvertToRequestss();
	
    // Default constructor
    BSS_Data() {
		all_int_trips.reserve(14000000);
		real_trips.reserve(200);
		trips.reserve(14000000);
        stations.reserve(5000);
	}
	
	std::map<int,std::vector<Trip>> real_trips_map;
	std::vector<std::vector<Trip>> real_trips;
	
	std::map<int,std::vector<Requests>> real_picks_map;
	std::map<int,std::vector<Requests>> real_dels_map;
	
	std::vector<TripDataAllInts> all_int_trips;
	std::vector<TripData> trips;
	std::vector<Station> stations;
};

class LambdaData {
public:
    LambdaData() : f(1), start_id(-1), end_id(-1){}
	int f = 1;//Multiplier > 1 to make more trips per scenario
    int start_id;
    int end_id;
    std::vector<double> duration_vec;
    std::vector<int> counts;
    
    // Function to modify counts based on the value of f
    LambdaData& ModifyCounts() {
        if (f > 1) {
            for (int& count : counts) {
                count *= f;
            }
        }
        return *this;
    }
    
    void Show() {
        printf("o:%d d:%d Durations[s]:", start_id, end_id);
        for (int i = 0; i < duration_vec.size(); i++)
            printf("%.1lf ", duration_vec[i]);
        printf("\nCounts:");
        for (int i = 0; i < counts.size(); i++)
            printf("%d ", counts[i]);
        printf("\n");
    }
};


class MyException : public std::exception {
public:
    MyException(const std::string& message) : message_(message) {}

    virtual const char* what() const noexcept {
        return message_.c_str();
    }

private:
    std::string message_;
};

//The main class ...
class Algorithm{
public:

    Algorithm() : AvgScenario(0) {
        lambda_object.reserve(1000000);
		lambdaPick_object.reserve(1300000);
		lambdaPick_object.reserve(1300000);
		tripIndex.reserve(14000000);
    }
	
	int randomSeed;
	RandomNumbers rn; 
	
	//Miscellaneous utilities
	void generateHistogramImage(std::vector<TripDataAllInts> & trips, int start_hour, int end_hour, std::string plot_name);	
		
	//To make Trip type instances
	void Make_avg_scenario_map(BSS_Data & BSSData, std::map<int, std::vector<Trip>> & trips_map);
	void MakeLambdaObject(std::vector<TripDataAllInts> & trips, int start_hour, int end_hour);	
	void Make_scenarios_map(int nb_scenarios, int nb_days_of_data, int start_hour, int end_hour, int instance_no, std::map<int, std::vector<Trip>> & trip_map) ;
	
	std::vector<LambdaData> lambda_object;
	std::unordered_map<std::pair<int, int>, std::vector<TripData>, PairHash> tripIndex;
	
	std::vector<LambdaData>* getLambdaPointer() {
        return &lambda_object;
    }
	std::vector<LambdaData>& getLambdaReference() {
		return lambda_object;
	}	
	
	//To make PickDel type instances
	void Make_avg_request_scenario_map(BSS_Data & BSSData, std::map<int, std::vector<Requests>> & picks_map, std::map<int, std::vector<Requests>> & dels_map);
	void MakePickDelLambdaObject(std::vector<TripDataAllInts> & trips, int start_hour, int end_hour);
	void Make_PickDelscenarios_map(int nb_scenarios, int nb_days_of_data, int start_hour, int end_hour, int instance_no, std::map<int, std::vector<Requests>> & pick_map, std::map<int, std::vector<Requests>> & del_map);
	void GetPickDelCounts(char type, std::map<int,std::vector<int>> & map, std::string time_str, int station_int, int start_hour_int, int end_hour_int);
	void GenReqs(char type,int scenario,int nb_days_of_data,int start_hour,int end_hour);
	
	std::vector<LambdaData> lambdaPick_object;
	std::vector<Requests> scenario_picks;
	std::vector<LambdaData> lambdaDel_object;
	std::vector<Requests> scenario_dels;
	
	std::vector<LambdaData>* getLambdaPickPointer() {
        return &lambdaPick_object;
    }
	std::vector<LambdaData>* getLambdaDelPointer() {
        return &lambdaDel_object;
    }
	std::vector<LambdaData>& getLambdaPickReference() {
		return lambdaPick_object;
	}
	std::vector<LambdaData>& getLambdaDelReference() {
		return lambdaDel_object;
	}
	
	int AvgScenario;
};


#endif