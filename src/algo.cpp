#include "algo.h"
#include <cmath> //std::log
#include <algorithm> //std::min, std::max
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <set>
#include "load.h"
#include <unordered_map>
#include <fstream>
#include <cstdlib> //system
#include <numeric> //std::accumulate
#include <locale>
#include <memory>

// Customized numpunct facet for thousands separation
struct separate_thousands : std::numpunct<char> {
    char_type do_thousands_sep() const override { return ','; }  // Separate with commas
    string_type do_grouping() const override { return "\3"; }     // Groups of 3 digits
};

// Function to format integer with custom thousands separation
std::string formatWithThousandsSeparator(int value) {
    std::stringstream ss;
    ss.imbue(std::locale(std::cout.getloc(), new separate_thousands())); // Apply the custom locale
    ss << value;
    return ss.str();
}

//Global variables
int nb_days_of_data;
int dayCounter;
int nb_instances;

std::unordered_map<std::string, int> dateToDayMap;

//kwrg parameters;
std::unordered_map<std::string, std::string> parameters; //The container

void writeTripsToFile(const std::map<int, std::vector<Trip>> & trips, std::string filename);

struct PickDropEvent
{
  int station; int time; int day; int qty;
};

bool pick_drop_event_sorter(const PickDropEvent &a, const PickDropEvent &b) {
    if (a.day != b.day) {
        return a.day < b.day;
    } else if (a.time != b.time) {
        return a.time < b.time;
    } else {
        // If day and time are the same, sort by station ID
        return a.station < b.station;
    }
}

int DayOfYear(int month, int day) {
    static const int daysBeforeMonth[] = {0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    return daysBeforeMonth[month] + day;
}

void CalculateInitialCapacities(BSS_Data & BSSData)
{
	std::vector<TripDataAllInts> & trips = BSSData.all_int_trips;
	std::vector<Station> & stations = BSSData.stations;
	
	std::vector<PickDropEvent> events;
	int min_day = 9999; int max_day = -1;
	for(size_t i=0;i<trips.size();i++)
	{
		PickDropEvent e;
		
		//if(trips[i].start_id == stations.size() || trips[i].end_id == stations.size())
		//{
			//trips[i].Show(); getchar();
		//}
		
		e.station = trips[i].start_id;
		e.time = 60*trips[i].start_hour + trips[i].start_minute;
		//First day is April 1st = Day 91
		e.day = DayOfYear(trips[i].start_month,trips[i].start_day) - 91;
		
		if( e.day < BSSData.dayCounter )
		{
			min_day = min_day > e.day ? e.day : min_day;
			e.qty = -1; //pickup 1 unit is removed
			events.push_back(e);			
		}
		
		e.station = trips[i].end_id;
		e.time = 60*trips[i].end_hour + trips[i].end_minute;
		e.day = DayOfYear(trips[i].end_month,trips[i].end_day) - 91;
		
		if(trips[i].end_month>7 || e.day > BSSData.dayCounter - 1) 
			continue;
		
		max_day = max_day < e.day ? e.day : max_day;				
		e.qty = 1; //drop 1 unit is added
		events.push_back(e);
	}
	std::sort(events.begin(), events.end(), pick_drop_event_sorter);

	int number_of_days = BSSData.dayCounter;
	std::vector<std::vector<int>> initial_qty(stations.size()); // station 0 is the depot and does not have trips.
	for(size_t i=0;i<stations.size();i++)
		initial_qty[i].resize( number_of_days, 0 );
	
	printf("Size of initial_qty rows:%d columns:%d nb_stations:%d min_day:%d max_day:%d difference:%d\n",(int)initial_qty.size(),(int)initial_qty[0].size(),(int)stations.size(),min_day,max_day,max_day-min_day); getchar();
	
	std::vector<std::vector<int>> current_qty = initial_qty;

	for(size_t i=0;i<events.size();i++)
	{
		printf("StatId:%d day:%d\n",events[i].station,events[i].day);
		current_qty[ events[i].station ][ events[i].day ] += events[i].qty;
		if( current_qty[ events[i].station ][ events[i].day ] < 0 )
		{
			initial_qty[ events[i].station ][ events[i].day ] -= current_qty[ events[i].station ][ events[i].day ];
			current_qty[ events[i].station ][ events[i].day ] = 0;
		}
	}

    FILE* outFile = fopen("initial_capacities.txt", "w");
    if (!outFile) {
        std::cerr << "Error opening file: initial_capacities.txt" << std::endl;
        return;
    }
	fprintf(outFile,"%d\n%d\n",(int)stations.size(),number_of_days);
	
    for(int i = 0; i < stations.size(); i++)
    {
        for(int j = 0; j < number_of_days; j++)
        {
            fprintf(outFile, "%d ", initial_qty[i][j]); // Print initial capacity for each day
        }
        fprintf(outFile, "\n");
    }

    fclose(outFile); 
}
	
void Algorithm::generateHistogramImage(std::vector<TripDataAllInts> & trips, int start_hour, int end_hour, std::string plot_name)
{
    printf("Generating image of trips:%d in time interval [%d,%d]\n",(int)trips.size(),start_hour,end_hour);
	
	// Create a map to count occurrences of hours
    std::map<int, int> hourCounts;
	std::map<int, int> dayCounts;
	std::vector<std::string> weekDays = {"Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"};

    // Iterate through the vector and count occurrences of each hour
	for(const auto& trip : trips)
	{
		hourCounts[trip.start_hour]++;
		if(trip.start_day == 1) // Saturday 1st, July 2023
			dayCounts[5]++;
		else if(trip.start_day == 2) // Sunday 2nd, July 2023
			dayCounts[6]++;
		else // Assuming trip.start_day == 3 for Monday, ..., 9 for Sunday
			dayCounts[(trip.start_day + 4) % 7]++; // dayCounts index based on the start day of the week
	}
	
	for(auto it = hourCounts.begin(); it!=hourCounts.end(); ++it)
	{
		//printf("Hour:%d Counts:%d\n",it->first,it->second);
		hourCounts[ it->first ] /= hourCounts.size(); //Average trip count per hour
		printf("Hour:%d AvgCounts:%d sizeOfMap:%d\n",it->first,it->second,(int)hourCounts.size());
	}
	for(auto it = dayCounts.begin(); it!=dayCounts.end(); ++it)
	{
		//printf("Hour:%d Counts:%d\n",it->first,it->second);
		dayCounts[ it->first ] /= dayCounts.size(); //Average trip count per hour
		printf("Day:%d AvgCounts:%d weekDay:%s sizeOfMap:%d\n",it->first,it->second,weekDays[it->first].c_str(),(int)dayCounts.size());
	}
	// Find the maximum value in the dayCounts map to scale
	int maxCount = 0;
	for (const auto& pair : dayCounts) {
		if (pair.second > maxCount) {
			maxCount = pair.second;
		}
	}	
		
    std::ofstream gnuplotFile("histogram.plt");
    //gnuplotFile << "set terminal png\n";
	gnuplotFile << "set terminal png size 800,600 enhanced font 'Arial,11'\n";
	//gnuplotFile << "set size ratio 0.75\n";
	gnuplotFile << "show terminal\n";
    gnuplotFile << "set output 'histogram.png'\n";
	//gnuplotFile << "set xlabel 'Day'\n";
    //gnuplotFile << "set xlabel 'Time of Day [h]'\n";
    gnuplotFile << "set ylabel 'Average Trip Start Time Count'\n";
    //gnuplotFile << "set title '" << plot_name.c_str() << "'\n";
    //Hour of day
	//gnuplotFile << "set xrange [" << start_hour-1 << ":" << end_hour << "]\n";
	//Day of the week
	gnuplotFile << "set xrange [-0.5:6.5]\n";
    gnuplotFile << "set yrange [0:" << maxCount * 1.3 << "]\n"; // 1.1 is a buffer to prevent bars from touching the top
	//gnuplotFile << "set yrange [10000:*]\n";
    //Day of the week
	//gnuplotFile << "set xtics 1\n";
    gnuplotFile << "set xtics (";
	for (size_t i = 0; i < weekDays.size(); ++i) {
		gnuplotFile << "'" << weekDays[i] << "' " << i << (i < weekDays.size() - 1 ? ", " : ")\n");
	}
	gnuplotFile << "set boxwidth 0.7\n";
	//Hour of the day
	//gnuplotFile << "plot '-' using 1:2 with boxes fill pattern 5 fc 'black' notitle\n";
	//Day of the week
	gnuplotFile << "plot '-' using 0:2:xtic(1) with boxes fill pattern 5 fc 'black' notitle\n";

	
	//Print xlabel: hours of the day
    /*for (int t = start_hour; t <= end_hour; t++) {
        if (hourCounts.find(t) != hourCounts.end()) {
            gnuplotFile << t << " " << hourCounts[t] << "\n";
        } else {
            gnuplotFile << t << " 0\n";
        }
    }*/
	//Print xlabel: days of the week
	for(auto it = dayCounts.begin(); it!=dayCounts.end(); ++it) {
		//std::cout << weekDays[it->first] << " " << it->second << std::endl;
		gnuplotFile << weekDays[it->first] << " " << it->second << "\n";
	}
    gnuplotFile << "e\n";
    gnuplotFile.close();
    int ret = system("gnuplot histogram.plt");
	
}

std::pair<double, double> computeMeanAndVariance(std::vector<double> duration_vec) {
    int n = duration_vec.size();

    if (n == 0) {
        // Handle the case of an empty vector
        throw MyException("Empty vector - cannot compute mean and variance.");
    }

    double sum = 0.0;
    double sum_of_squares = 0.0;

    // Calculate the sum and sum of squares
    for (double x : duration_vec) {
        sum += x;
        sum_of_squares += x * x;
    }

    double mean = sum / (double)n;
    double variance = (sum_of_squares / n) - (mean * mean);

    return std::make_pair(mean, variance);
}

std::pair<double,double> computeAlphaAndBeta(double mean, double variance)
{
	double alpha = 0.0; double beta=0.0;
	alpha = mean*( ( mean*(1-mean) )/variance -1 );
	beta = (1-mean)*( ( mean*(1-mean) )/variance -1 );
	
	if (alpha < 0.0 || beta < 0.0) {
		std::string errorMsg;
		errorMsg = "Could not compute Alpha, Beta -> Alpha: " + std::to_string(alpha) + " Beta: " + std::to_string(beta);

		throw MyException(errorMsg);
	}
			
	return std::make_pair(alpha,beta);
}

std::vector<double> scale_to_01(const std::vector<double>& vec) {
    double min_val = *std::min_element(vec.begin(), vec.end());
    double max_val = *std::max_element(vec.begin(), vec.end());

    // Check for the case where all values in the vector are the same
    if (min_val == max_val) {
        if (vec.size() > 1) {
            // Instead of throwing an exception, return a vector of zeros
            return std::vector<double>(vec.size(), 0.0);
        } else {
            // If there's only one element, return a vector with a single 0.0
            return std::vector<double>(1, 0.0);
        }
    }

    std::vector<double> scaled_vec;
    scaled_vec.reserve(vec.size());  // Reserve memory for performance
    for (double x : vec) {
        double scaled_value = (x - min_val) / (max_val - min_val);
        scaled_vec.push_back(scaled_value);
    }

    return scaled_vec;
}

/*
std::vector<double> scale_to_01(std::vector<double> vec) {
    double min_val = *std::min_element(vec.begin(), vec.end());
    double max_val = *std::max_element(vec.begin(), vec.end());

    if (min_val == max_val && vec.size() > 1) {
		std::string errorMsg;
		errorMsg = "Could not scale duration_vec. Min_element:" + std::to_string(min_val) + " is equal to Max_element:" + std::to_string(max_val)	;

		throw MyException(errorMsg);
    }
	
	
    std::vector<double> scaled_vec;
    for (double x : vec) {
        double scaled_value = (x - min_val) / (max_val - min_val);
        scaled_vec.push_back(scaled_value);
    }

    return scaled_vec;
}
*/

double descale_from_01(double scaled_value, double min_val, double max_val ){
    double descaled_value = scaled_value * (max_val - min_val) + min_val;
    return descaled_value;
}

void Algorithm::MakeLambdaObject(std::vector<TripDataAllInts> & trips, int start_hour, int end_hour) 
{
	int T = end_hour - start_hour;
	printf("Making lambda_object ...\n");
	printf("Size of trip object:%lu Sample:\n",sizeof(trips[0])); trips[0].Show();
	lambda_object.clear();
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
	clock_t currentTime = clock();
	
	// Initialize the maps with unique pairs
	//std::pair<int, int>: start_id, end_id of trip
	//std::pair<std::vector<double>: duration in seconds of the trip
	//std::vector<int>: counts of the trips in the interval hours of 8am,9am,...,10pm
	std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<int>>> unique_pairs;
	double trips_processed = 0.0;

	for (const auto& trip : trips) {
		++trips_processed;
		if (trips_processed / trips.size() * 10 == 0) {
			double percentage = (trips_processed / trips.size()) * 100.0;
			currentTime = clock();
			double elapsedSeconds = (double)(currentTime - startTime) / CLOCKS_PER_SEC;
			std::cout << "Processed " << (int)(percentage + 0.5) << "% of trips. ElapsedTime:" << std::fixed << std::setprecision(2) << elapsedSeconds << std::endl;
		}

		std::pair<int, int> pair(trip.start_id, trip.end_id);
		auto& unique_pair = unique_pairs[pair];
		unique_pair.first.push_back(trip.duration);
		
		// Initialize unique_pair.second only if it's not already initialized
		if (unique_pair.second.empty()) {
			unique_pair.second = std::vector<int>(T, 0);
		}
		
		int hour = trip.start_hour;
		if (hour >= start_hour && hour < end_hour) {
			unique_pair.second[hour - start_hour]++;
		} else {
			trip.Show();
			throw MyException("Wrong hour extracted at the following trip. Did you forget to filter the trips in relation to the times...");
		}		
	}

	currentTime = clock();
	printf("Unique pairs:%d Time taken[s]:%.1lf\n",(int)unique_pairs.size(),(double)((currentTime - startTime) / CLOCKS_PER_SEC));
	//exit(1);
	
	for (auto it = unique_pairs.begin(); it != unique_pairs.end(); ++it) {
		
		LambdaData lambda;
		const auto& unique_pair = it->first; // Access the key pair

		lambda.start_id = unique_pair.first; // Access the first integer
		lambda.end_id = unique_pair.second; // Access the second integer

		// Accessing vectors from the map using the key
		const auto& vectors_pair = it->second; // Access the pair of vectors
		lambda.duration_vec = vectors_pair.first; // Access the vector of durations
		lambda.counts = vectors_pair.second; // Access the vector of counts
		
		lambda_object.push_back(lambda);
	}

	
	//To check
	int count=0;
	for(int i=0; i<lambda_object.size(); i++)
		for(int t=0; t<T; t++)
		{
			count+=lambda_object[i].counts[t];
		}
	
	if(count!=trips.size())
	{
		printf("Lambda_object has different trips than trips object. Exiting ...\n");
		exit(1);
	}else{printf("Correct count of trips from lambda_object equals nb trips! count from lambda_object:%d nb_trips:%d\n",count,(int)trips.size());}
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Total Elapsed Time for lambda_object: " << elapsedSeconds << " seconds" << std::endl;
	printf("Size of lambda_object:%d\n",(int)lambda_object.size());

}

void Algorithm::MakePickDelLambdaObject(std::vector<TripDataAllInts> & trips, int start_hour, int end_hour)
{
	int T = end_hour - start_hour;
	printf("Making PickDel lambda_object in the interval [%dh,%dh]...\n",start_hour,end_hour);
	printf("Size of trip:%lu Sample trip object:\n",sizeof(trips[0])); trips[0].Show();
	lambda_object.clear();
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
	clock_t currentTime = clock();
	
	// Initialize the maps with unique pairs
	//std::pair<int, int>: start_id, end_id of trip
	//std::pair<std::vector<double>: duration in seconds of the trip
	//std::vector<int>: counts of the trips in the interval hours of 8am,9am,...,10pm
	std::map<int,std::vector<int>> start_stat_map;
	std::map<int,std::vector<int>> end_stat_map;
	
	for (int i = 0; i < trips.size(); i++) 
	{
		if ((i + 1) % (trips.size() / 10) == 0) 
		{			
			// Print your output here, e.g., percentage completion
			double percentage = (double)(i + 1) / trips.size() * 100.0;
			currentTime = clock();
			double elapsedSeconds = (double)(currentTime - startTime) / CLOCKS_PER_SEC;
			std::cout << "Processed " << (int)(percentage + 0.5) << "% of trips. ElapsedTime:" << std::fixed << std::setprecision(2) << elapsedSeconds << std::endl;
		}
		
		
		
		auto &start_stat_ref = start_stat_map[ trips[i].start_id ];
		auto &end_stat_ref = end_stat_map[ trips[i].end_id ];

		if (start_stat_ref.empty()) start_stat_ref = std::vector<int>(T, 0);
		if (end_stat_ref.empty()) end_stat_ref = std::vector<int>(T, 0);

		int trip_start_hour = trips[i].start_hour;
		int trip_end_hour = trips[i].end_hour;

		if (start_hour < end_hour && trip_start_hour >= start_hour && trip_start_hour < end_hour && trip_end_hour >= start_hour && trip_end_hour < end_hour) {
			start_stat_ref[trip_start_hour - start_hour]++;
			end_stat_ref[trip_end_hour - start_hour]++;
		} else {
			trips[i].Show();
			throw MyException("Wrong hour extracted at the following trip. Did you forget to filter the trips in relation to the times...");
		}			
	}

	currentTime = clock();
	printf("Picks:%d Dels:%d Time taken[s]:%.1lf\n",(int)start_stat_map.size(),(int)end_stat_map.size(),(double)((currentTime - startTime) / CLOCKS_PER_SEC));
	//exit(1);
	
	//Always preincrement (never postincrement)
	for(auto it=start_stat_map.begin(); it!=start_stat_map.end(); ++it)
	{
		LambdaData lambda;
		lambda.start_id = it->first;
		lambda.counts = it->second;
		lambdaPick_object.push_back(lambda);
		//lambda.Show();
	}
	for(auto it=end_stat_map.begin(); it!=end_stat_map.end(); ++it)
	{
		LambdaData lambda;
		lambda.start_id = it->first;
		lambda.counts = it->second;
		lambdaDel_object.push_back(lambda);
	}

	//To check
	int pick_count=0; int del_count=0;
	for(int i=0; i<lambdaPick_object.size(); i++)
		for(int t=0; t<T; t++)
		{
			pick_count+=lambdaPick_object[i].counts[t];
		}
	for(int i=0; i<lambdaDel_object.size(); i++)
		for(int t=0; t<T; t++)
		{
			del_count+=(std::abs(lambdaDel_object[i].counts[t]));
		}
	
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Total Elapsed Time for lambda_object: " << elapsedSeconds << " seconds" << std::endl;
	printf("Size of lambdaPick_object:%d lambdaDel_object:%d\n",(int)lambdaPick_object.size(),(int)lambdaDel_object.size());
	
}

void Algorithm::Make_PickDelscenarios_map(int nb_scenarios, int nb_days_of_data, int start_hour, int end_hour, int instance_no, std::map<int, std::vector<Requests>> & pick_map, std::map<int, std::vector<Requests>> & del_map)
{	
	std::vector<int> start_hours;
	for(int i=start_hour; i<end_hour; i++) start_hours.push_back(i);
	int nb_time_intervals = (int)start_hours.size();
	int end_of_day_time = 60*nb_time_intervals;
		
	printf("In PickDel Make_scenarios with nb_time_intervals:%d nb_days_of_data:%d start_of_day_hour (t=0):%02dh end_of_day_time[minutes]:%d seed:%d\n",nb_time_intervals,nb_days_of_data,start_hours[0],end_of_day_time,randomSeed);

	for(int e=0;e<nb_scenarios;e++)
	{
		scenario_picks.clear();
		scenario_dels.clear();
		
		GenReqs(P_CHAR,e,nb_days_of_data,start_hour,end_hour);
		GenReqs(Q_CHAR,e,nb_days_of_data,start_hour,end_hour);
		
		//printf("scenario:%d picks:%d dels:%d\n",e,(int)scenario_picks.size(),(int)scenario_dels.size());
		
		pick_map[e]=scenario_picks;
		del_map[e]=scenario_dels;			
	}
}

void Algorithm::GenReqs(char type, int scenario, int nb_days_of_data, int start_hour, int end_hour)
{
	std::vector<LambdaData>& lambdaRef = (type == P_CHAR)? getLambdaPickReference() : getLambdaDelReference();
	int T = end_hour - start_hour;
	for(int i=0; i<lambdaRef.size(); i++)
	{
		if(i%50000 == 0)
		{
			char type = 'n';
			for(int t=0;t<lambdaRef[i].counts.size();t++)
			{
				if(lambdaRef[i].counts[t]>0)
				{ type=P_CHAR; break;}
				if(lambdaRef[i].counts[t]<0)
				{ type=Q_CHAR; break;}
			}
				
			//printf("type:%c sce:%d looping at stat:%d ...\n",type,scenario,i);
			//getchar();
		}
		if(lambdaRef[i].counts.size() != T)
		{ printf("The following lambda object is erroneous with a different number of counts:%d than T:%d. Exiting ...\n",(int)lambdaRef[i].counts.size(),T); exit(1); }
	
		LambdaData lambda = lambdaRef[i];
		for(int start_h=0; start_h<T; start_h++)
		{
			if(lambda.counts[start_h] == 0) continue;
			double lambda_val = std::abs(lambda.counts[start_h])/(60.0*(double)nb_days_of_data);
			//printf("lambda_val:%.2lf count:%d\n",lambda_val,lambda.counts[start_h]);
			//getchar();
			
			int t=0;
			while(true)
			{
				double u = rn.rand01();
				
				int t_next = t + std::ceil(-std::log(1 - u) / lambda_val);
				//printf("u:%.1lf t_next:%d\n",u,t_next);
				//getchar();
				
				if(t_next > 60 ) break;
				
				Requests req;
				if(type == P_CHAR) req.type = P_CHAR;
				if(type == Q_CHAR) req.type = Q_CHAR;
				req.time = (start_h*60)+t_next;
				req.scenario = scenario;
				req.stat_id = lambdaRef[i].start_id; //Pick, Del requests have only start_id (not end_id)
				
				//req.Show(); //getchar();
				
				if(type == P_CHAR)
				{
					scenario_picks.push_back(req);
					if((int)scenario_picks.size()%5000 == 0)
						printf("type:%c sce:%d generated requests:%d ...\n",type,scenario,(int)scenario_picks.size());
				}					
				if(type == Q_CHAR)
				{
					scenario_dels.push_back(req);
					if((int)scenario_dels.size()%5000 == 0)
						printf("type:%c sce:%d generated requests:%d ...\n",type,scenario,(int)scenario_dels.size());
				}			
				t=t_next;
			}			
		}

	}
}

void Algorithm::Make_avg_request_scenario_map(BSS_Data & BSSData, std::map<int, std::vector<Requests>> & picks_map, std::map<int, std::vector<Requests>> & dels_map)
{
	printf("In Make_avg_request_scenario ...\n");
	
	int AvgNbTrips = BSSData.all_int_trips.size() / BSSData.dayCounter;

	for(int i=0;i<AvgNbTrips;i++)
	{
		int index = rn.randInt(0, BSSData.all_int_trips.size() - 1);

		Requests pick;
		pick.type = 'p'; 
		pick.stat_id = BSSData.all_int_trips[ index ].start_id; 
		pick.time = 60*BSSData.all_int_trips[ index ].start_hour + BSSData.all_int_trips[ index ].start_minute; 
		pick.scenario = 0;
		
		Requests del;
		del.type = 'd'; 
		del.stat_id = BSSData.all_int_trips[ index ].end_id; 
		del.time = 60*BSSData.all_int_trips[ index ].end_hour + BSSData.all_int_trips[ index ].end_minute; 
		del.scenario = 0;		
		
		picks_map[ 0 ].push_back( pick );
		dels_map[ 0 ].push_back( del );
		//pick.Show();
		//del.Show();
	}
	printf("The average scenario has %d picks and %d dels\n",(int)picks_map[0].size(),(int)dels_map[0].size());
}

void Algorithm::Make_avg_scenario_map(BSS_Data & BSSData, std::map<int, std::vector<Trip>> & trips_map)
{
	printf("In Make_avg_scenario ...\n");
	
	int AvgNbTrips = BSSData.all_int_trips.size() / BSSData.dayCounter;

	for(int i=0;i<AvgNbTrips;i++)
	{
		int index = rn.randInt(0, BSSData.all_int_trips.size() - 1);
		Trip trip; 
		trip.origin = BSSData.all_int_trips[ index ].start_id;
		trip.destination = BSSData.all_int_trips[ index ].end_id;
		trip.scenario = 0;
		
		//if(BSSData.all_int_trips[ index ].duration>60*14) continue;
		
		trip.start_time = 60*BSSData.all_int_trips[ index ].start_hour + BSSData.all_int_trips[ index ].start_minute;
		trip.end_time = 60*BSSData.all_int_trips[ index ].end_hour + BSSData.all_int_trips[ index ].end_minute;
		
		trips_map[ 0 ].push_back( trip );
		//trip.Show();
	}
	printf("The average scenario has %d trips.\n",(int)trips_map[0].size());
}


void Algorithm::Make_scenarios_map(int nb_scenarios, int nb_days_of_data, int start_hour, int end_hour, int instance_no, std::map<int, std::vector<Trip>> & trips_map)
{
	std::vector<int> start_hours;
	for(int i=start_hour; i<end_hour; i++) start_hours.push_back(i);
	int nb_time_intervals = (int)start_hours.size();
	int end_of_day_time = 60*nb_time_intervals;
	
	printf("In Make_scenarios with nb_time_intervals:%d start_of_day_hour (t=0):%02dh end_of_day_time[minutes]:%d seed:%d\n",nb_time_intervals,start_hours[0],end_of_day_time,randomSeed);	
	
	double avg_trips = 0.0; int no_beta_params = 0; int out_of_time = 0; int max_end_time = 0;
	for(int e=0; e<nb_scenarios; e++)
	{
		int nbScenarioTrips=0;
		int nbScenarioNoBeta=0;
		for(int i=0; i<lambda_object.size(); i++)
		{
			if(i%50000 == 0)
				printf("sce:%d looping at pair:%d ...\n",e,i);
			
			LambdaData od = lambda_object[i];
			//if(od.duration_vec.size()==1) continue;
			if(od.counts.size() < start_hours.size())
			{
				printf("The following OD pair has less lambdas than the amount of specified time intervals. Exiting ..."); od.Show(); exit(1);
			}
			for(int start_h=0; start_h<start_hours.size(); start_h++)
			{
				double lambda_val = od.counts[start_h]/(60.0*(double)nb_days_of_data);
				if(od.counts[start_h] == 0) continue;
							
				int t=0;
				while(true)
				{
					double u = rn.rand01();
					//printf("u:%.1lf\n",u);
					//getchar();
					int t_next = t + std::ceil(-std::log(1 - u) / lambda_val);
					if(t_next > 60) break;
					
					std::vector<double> duration_vec = od.duration_vec;
					
					double alpha = 0.0; double beta = 0.0;
					std::pair<double, double> result;
					std::pair<double,double> parameter;
					try{
				
						std::vector<double> duration_vec_01 = scale_to_01(duration_vec);
						// Fit the beta distribution using the method of moments
						result = computeMeanAndVariance(duration_vec_01);
						parameter = computeAlphaAndBeta(result.first,result.second);
						alpha = parameter.first; beta = parameter.second; // Condition: μ(1−μ)>σ^2 to get positive alpha, beta
						
					}catch (MyException& ex) {
						std::cout << "Caught exception: " << ex.what() << " for the following OD pair:" <<std::endl;
						od.Show();
						//getchar();
						break;
					}
					
					double min_val = *std::min_element(duration_vec.begin(), duration_vec.end());
					double max_val = *std::max_element(duration_vec.begin(), duration_vec.end());
					
					//y is in minutes now
					int y = -1;
					double sum = -1.0; double average = 0.0;
					if(alpha < 0.0001 || beta < 0.0001 && max_val - min_val > 1.0)//You need a alpha, beta, but could not compute them
					{
						/*od.Show(); 
						std::cout << "Mean:" << result.first << ", Variance:" << result.second << std::endl;
						std::cout << "Alpha:" << parameter.first << ", Beta:" << parameter.second << std::endl;
						getchar();*/
						nbScenarioNoBeta++;
						no_beta_params++;
						//if(nbScenarioNoBeta%2==0)//We accept half of the NoParameter situation ...
						// We accept about 40% of the NoParameter situation, skip 60%
						//if (nbScenarioNoBeta % 10 < 6) {
							//break;
						//}
						// We accept about 55% of the NoParameter situation, skip 45%
						if (nbScenarioNoBeta % 20 < 11) {
							break;
						}						
						sum = std::accumulate(duration_vec.begin(), duration_vec.end(), 0.0);
						average = sum / (duration_vec.size() * 60);
						 //break;
					}					
						
					double y_01 = y < 0 ? rn.randBeta(alpha,beta) : -1.0;
					if(max_val - min_val > 1.0 && duration_vec.size() > 1 && y < 0)
						y = std::ceil(descale_from_01(y_01, min_val, max_val)) / 60;
					if(max_val - min_val <= 1.0 && y < 0)
						y = duration_vec[0] / 60;
					
					if(y<0)
					{
						printf("Negative y:%d y_01:%.2lf alpha:%.6lf beta:%.6lf\nYou made a mistake somewhere! Exiting ...",y,y_01,alpha,beta); od.Show(); exit(1);
					}
					
					//Uncomment if you don't want to accept all generated trips
					//if( y>=1 && t_next + y <= end_of_day_time )
					//if( y>=1 && start_h*60 + t_next + y <= end_of_day_time )
					if( y>=1 && start_h*60 + t_next + y <= 60*16 ) //Before the start of next operational day
					{
						//Generate OD trip
						Trip trip;
						//trip.origin=od.start_id; trip.destination=od.end_id; trip.start_time=t_next; trip.end_time=t_next + y; trip.scenario=e;
						trip.origin=od.start_id; trip.destination=od.end_id; trip.start_time=start_h*60 + t_next; trip.end_time=start_h*60 + t_next + y; trip.scenario=e;
						//scenario_trips.push_back(trip);
						trips_map[e].push_back(trip);
						nbScenarioTrips++;
						
						if(nbScenarioTrips%5000 == 0)
							printf("sce:%d generated trips:%d ...\n",e,nbScenarioTrips);
						
						//if((int)scenario_trips.size()%5000 == 0)
							//printf("sce:%d generated trips:%d ...\n",e,(int)scenario_trips.size());
						max_end_time = std::max( max_end_time, start_h*60 + t_next + y );
					}
					else{
						out_of_time++;
					}
					t=t_next;
				}
				
			}
		}
		avg_trips += trips_map[e].size();
		/*if (e % 200 == 199) {
			// Write trips to file after every 200 scenarios
			char filename[100];
			snprintf(filename, sizeof(filename), "../instances/%d_trips_batch_%d.txt", instance_no, e / 200);
			printf("Writing batch %d scenarios (from %d to %d) to %s\n", e / 200, e - 199, e, filename);
			writeTripsToFile(trips_map, filename);

			// Clear the map after writing to the file
			trips_map.clear();
		}*/
		
	}
	
	/*if(nb_scenarios > 200 && !trips_map.empty())
	{
		// Write any remaining scenarios
		if (!trips_map.empty()) {
			char filename[100];
			snprintf(filename, sizeof(filename), "../instances/%d_trips_batch_last.txt",instance_no);
			printf("Writing remaining scenarios to %s\n", filename);
			writeTripsToFile(trips_map, filename);
		}		
	}*/
	printf("In MakeSce AvgTrips:%.1lf AvgNoBetaParams:%.1lf AvgOOT:%.1lf MaxEndTime:%d\n",avg_trips/(double)nb_scenarios,no_beta_params/(double)nb_scenarios,out_of_time/(double)nb_scenarios,max_end_time); //getchar();
}

// Helper function to convert time string to minutes since midnight
int TimeStringToMinutes(const std::string& time_str) 
{
    std::istringstream stream(time_str);
    int hours, minutes;
    char colon;
    stream >> hours >> colon >> minutes;

    return hours * 60 + minutes;
}

void BSS_Data::ConvertToTrips()
{
    printf("Making the real trips ...\n");
	int realTripCounter=0;
	for (const auto& trip : trips) 
	{	
		Trip generated_trip;
        generated_trip.origin = trip.start_id;
        generated_trip.destination = trip.end_id;

        // Convert start_time_str and end_time_str to minutes
        generated_trip.start_time = TimeStringToMinutes(trip.start_time_str);
        generated_trip.end_time = TimeStringToMinutes(trip.end_time_str);
		//Discard the trips that start on one day an end on the next ...
		if(generated_trip.end_time>(60*24)) continue;
		
		generated_trip.scenario = dateToDayMap[trip.start_date_str];
		real_trips_map[generated_trip.scenario].push_back(generated_trip);
		realTripCounter++;
	}
	
	int max_scenario = 0;
	if (!real_trips_map.empty()) {
		max_scenario = real_trips_map.rbegin()->first;
	}
	printf("Nb real trip:%d Maximum scenario in real trips: %d\n", realTripCounter,max_scenario);

}

void BSS_Data::ConvertToRequestss()
{
    printf("Making the real pickups and deliveries ...\n");
	for (const auto& trip : trips) 
	{	
		Requests pickup;
		pickup.type = 'p';
		pickup.stat_id = trip.start_id;
		pickup.time = TimeStringToMinutes(trip.start_time_str);
		pickup.scenario=dateToDayMap[trip.start_date_str];
		
		Requests delivery;
		delivery.type = 'd';
		delivery.stat_id = trip.end_id;
		delivery.time = TimeStringToMinutes(trip.end_time_str);
		delivery.scenario=dateToDayMap[trip.end_date_str];

		if(pickup.time<=(60*24))
			real_picks_map[pickup.scenario].push_back(pickup);
		if(delivery.time<=(60*24))
			real_dels_map[delivery.scenario].push_back(delivery);		
    }
	
	int max_scenario_picks = 0; int max_scenario_dels =0; int max_scenario=0;
	if (!real_picks_map.empty()) {
		max_scenario_picks = real_picks_map.rbegin()->first;
	}
	if (!real_dels_map.empty()) {
		max_scenario_dels = real_dels_map.rbegin()->first;
	}
	max_scenario=std::max(max_scenario_picks,max_scenario_dels);
	
	printf("Maximum scenario in real trips: %d\n", max_scenario);

}

void writeRequestsAndStationsToFile(const std::map<int,std::vector<Requests>> & picks, const std::map<int,std::vector<Requests>> & dels, const std::vector<Station> & stations, int Qtot, std::string filename)
{
	FILE* output_file = fopen(filename.c_str(), "w");

    if (!output_file) {
        std::cerr << "Error: Unable to open the requests and stations output file." << std::endl;
        return;
    }
	printf("First ... Printing output stations to file\n");
    fprintf(output_file, "%d\n", (int)stations.size());
    fprintf(output_file, "%d\n",Qtot);

    // Initialize random number generator
    std::srand((unsigned int)(std::time(nullptr)));

    for (const Station& station : stations)
        fprintf(output_file, "%d ", station.cap);
	fprintf(output_file,"\n");
	printf("Next ... Printing output requests ...\n");

    for (auto it = picks.begin(); it != picks.end(); ++it) {
        int scenario = it->first;
        const std::vector<Requests>& requests_vec = it->second;

        for (const auto& req : requests_vec) {
			fprintf(output_file, "%c %d %d %d\n",
                    req.type, req.stat_id, req.time, scenario);
        }
    }
    for (auto it = dels.begin(); it != dels.end(); ++it) {
        int scenario = it->first;
        const std::vector<Requests>& requests_vec = it->second;

        for (const auto& req : requests_vec) {
			fprintf(output_file, "%c %d %d %d\n",
                    req.type, req.stat_id, req.time, scenario);
        }
    }	
    fclose(output_file);	
}

void WriteTripsAndStationsToFile(const std::map<int,std::vector<Trip>> & trips, const std::vector<Station> & stations, int Qtot, std::string filename)
{
	FILE* output_file = fopen(filename.c_str(), "w");

    if (!output_file) {
        std::cerr << "Error: Unable to open the trips and stations output file." << std::endl;
        return;
    }
	printf("First ... Printing output stations to file\n");
    fprintf(output_file, "%d\n", (int)stations.size());
    fprintf(output_file, "%d\n",Qtot);

    // Initialize random number generator
    std::srand((unsigned int)(std::time(nullptr)));

    for (const Station& station : stations)
        fprintf(output_file, "%d ", station.cap);
	fprintf(output_file,"\n");
	printf("Next ... Printing output trips\n");

    for (auto it = trips.begin(); it != trips.end(); ++it) {
        int scenario = it->first;
        const std::vector<Trip>& trips_vec = it->second;

        for (const auto& trip : trips_vec) {
			fprintf(output_file, "%d %d %d %d %d\n",
                    trip.origin, trip.destination, trip.start_time, trip.end_time, scenario);
        }
    }
    fclose(output_file);	
}

void writeRequestsToFile(const std::map<int,std::vector<Requests>> & picks, const std::map<int,std::vector<Requests>> & dels, std::string filename)
{
    FILE* file = fopen(filename.c_str(), "w");
    if (!file) {
        std::cerr << "Failed to open requests output file: " << filename << std::endl;
        return;
    }
    printf("Printing output requests to file:%s\n", filename.c_str());

	for (auto it = picks.begin(); it != picks.end(); ++it) {
        int scenario = it->first;
        const std::vector<Requests>& requests_vec = it->second;

    for (const auto& req : requests_vec) {
			fprintf(file, "%c %d %d %d\n",
                    req.type, req.stat_id, req.time, scenario);
        }
    }
    for (auto it = dels.begin(); it != dels.end(); ++it) {
        int scenario = it->first;
        const std::vector<Requests>& requests_vec = it->second;

        for (const auto& req : requests_vec) {
			fprintf(file, "%c %d %d %d\n",
                    req.type, req.stat_id, req.time, scenario);
        }
    }
	fclose(file);
}

void writeTripsToFile(const std::map<int, std::vector<Trip>> & trips, std::string filename) {
    FILE* file = fopen(filename.c_str(), "w");
    if (!file) {
        std::cerr << "Failed to open trips output file: " << filename << std::endl;
        return;
    }
    printf("Printing output trips to file:%s\n", filename.c_str());

    for (auto it = trips.begin(); it != trips.end(); ++it) {
        int scenario = it->first;
        const std::vector<Trip>& trips_vec = it->second;

        for (const auto& trip : trips_vec) {
            fprintf(file, "%d %d %d %d %d\n",
                    trip.origin, trip.destination, trip.start_time, trip.end_time, scenario);
        }
    }
    fclose(file);
}

void writeStationsToFileTargets(const std::vector<Station>& stations, const std::string& filename, const int& Qtot)
{
    FILE* output_file = fopen(filename.c_str(), "w");

    if (!output_file) {
        std::cerr << "Error: Unable to open the single station output file." << std::endl;
        return;
    }
	printf("Printing output stations to file:%s for Target problem instances\n",filename.c_str());
    fprintf(output_file, "%d\n", (int)stations.size());
    fprintf(output_file, "%d\n",Qtot);

    // Initialize random number generator
    std::srand((unsigned int)(std::time(nullptr)));

    for (const Station& station : stations)
        fprintf(output_file, "%d ", station.cap);
	fprintf(output_file,"\n");
	
    fclose(output_file);	
}

void writeStationsToFileRebalancing( std::vector<Station>& stations, const std::string& filename, const int& Qtot )
{
	    FILE* output_file = fopen(filename.c_str(), "w");

    if (!output_file) {
        std::cerr << "Error: Unable to open the single station output file." << std::endl;
        return;
    }
	printf("Printing output stations to file:%s for Rebalancing problem instances\n",filename.c_str());
    fprintf(output_file, "%d\n", (int)stations.size());
    fprintf(output_file, "%d\n",Qtot);

    for (const Station& station : stations)
	{
		if(station.name == "depot")
			fprintf(output_file, "%d %.4lf %.4lf\n",station.cap,station.lat,station.lon);
		else fprintf(output_file, "%s %d %.6lf %.6lf\n",station.bss_station_id.c_str(),station.cap,station.lat,station.lon);
	}
    fclose(output_file);	
}

void ProcessKeywordParameters(int argc, char* argv[]) {
    // Loop through the arguments starting from index 1
    for (int i = 1; i < argc; ++i) {
        std::string argument = argv[i];

        // Check if the argument contains "="
        size_t pos = argument.find('=');
        if (pos != std::string::npos) {
            // Split the argument into parameter name and value
            std::string parameter_name = argument.substr(0, pos);
            std::string parameter_value = argument.substr(pos + 1);

            parameters[parameter_name] = parameter_value;
        }
    }

    // Print the parsed parameters
    for (const auto& param : parameters) {
        std::cout << "Parameter Name: " << param.first << ", Value: " << param.second << std::endl;
    }
}


int main(int argc, char* argv[])
{
	ProcessKeywordParameters(argc, argv);
	
	if (parameters.find("trips_data_file") == parameters.end() ||
		parameters.find("stations_data_file") == parameters.end() ||
		parameters.find("instance_type") == parameters.end() ||
		parameters.find("output_trips") == parameters.end() || parameters.find("output_file")==parameters.end() || parameters.find("instance_format")==parameters.end()) {
		
		printf("Wrong usage: You did not specify the trips, station data and/or the output trips to be random or the real ones. Exiting. Hit `make usage' for options ... \n");
		exit(1);
	}
	auto tripsData = parameters.find("trips_data_file");
	auto stationsData = parameters.find("stations_data_file");
	auto instanceType = parameters.find("instance_type");
	auto outputTrips = parameters.find("output_trips");
	auto outputFile = parameters.find("output_file");
	auto instanceFormat = parameters.find("instance_format");
	auto seed = parameters.find("seed");
	auto rebalancing = parameters.find("rebalancing");
	auto avgScenario = parameters.find("AvgScenario");
	
	if(instanceFormat->second == "trips")
		printf("Executing generator for %s and %s from BSS of type %s ... Output trips:%s  and the format will be:%s\n",tripsData->second.c_str(),stationsData->second.c_str(),instanceType->second.c_str(),outputTrips->second.c_str(),instanceFormat->second.c_str());
	if(instanceFormat->second == "requests")
		printf("Executing generator for %s and %s from BSS of type %s ... Output trips:%s and the format will be:%s\n",tripsData->second.c_str(),stationsData->second.c_str(),instanceType->second.c_str(),outputTrips->second.c_str(),instanceFormat->second.c_str());

	Algorithm algo;
	if(parameters.find("avg_scenarion") != parameters.end())
		if(std::stoi( avgScenario->second ) == 1)
		{
			printf("You set the AvgScenario to 1. The generator will compute it and the exit ...\n");
			algo.AvgScenario = 1;
		}
		
	algo.randomSeed = std::stoi(seed->second);
	printf("Random seed given:%d\nRemember that the seed needs to always change!\n",algo.randomSeed);
	algo.rn.init(algo.randomSeed);
	
	clock_t startTime = clock();
	int start_hour = std::stoi(parameters["start_hour"]); int end_hour = std::stoi(parameters["end_hour"]); 
	
	int nb_scenarios=0; nb_instances=0; int max_stations_allowed=0;
	if (parameters.find("nb_scenarios") != parameters.end()) {
        nb_scenarios = std::stoi(parameters["nb_scenarios"]);
    } else {
        std::cout << "Parameter 'nb_scenarios' not provided." << std::endl;
    }
	if (parameters.find("nb_instances") != parameters.end()) {
        nb_instances = std::stoi(parameters["nb_instances"]);
    } else {
        std::cout << "Parameter 'nb_instances' not provided." << std::endl;
    }
	if (parameters.find("max_stations_allowed") != parameters.end()) {
        max_stations_allowed = std::stoi(parameters["max_stations_allowed"]);
    } else {
        std::cout << "Parameter 'max_stations_allowed' not provided." << std::endl;
    }

	nb_days_of_data = 0;
	int Qtot = 0;
	BSS_Data BSSData; 
	Load load_object;
	
	if(max_stations_allowed==0)
	{
		printf("Setting nb_instances:%d nb_scenarios:%d start_hour:%d end_hour:%d\n",nb_instances,nb_scenarios,start_hour,end_hour);
		
		if(instanceType->second=="washington" || instanceType->second=="newyork" || instanceType->second=="boston" || instanceType->second=="chicago" || instanceType->second=="sanfrancisco")
			load_object.Load_BSS_Data(BSSData,tripsData->second,stationsData->second,start_hour,end_hour,instanceType->second,0); // 0 is for the cap on stations allowed		
		if(instanceType->second=="montreal")
			load_object.Load_Bixi(BSSData,tripsData->second,stationsData->second,start_hour,end_hour);
		if(instanceType->second=="mexicocity")
			load_object.Load_Ecobici(BSSData,tripsData->second,stationsData->second,start_hour,end_hour);
		if(instanceType->second=="toronto")
			load_object.Load_Toronto(BSSData,tripsData->second,stationsData->second,start_hour,end_hour);
		if(instanceType->second=="madrid")
			load_object.Load_Madrid(BSSData,tripsData->second,stationsData->second,start_hour,end_hour);
	}
	else{ 
		printf("Setting nb_instances:%d nb_scenarios:%d max_stations_allowed:%d start_hour:%d end_hour:%d\n",nb_instances,nb_scenarios,max_stations_allowed,start_hour,end_hour);
		
		if(instanceType->second=="montreal")
			load_object.Load_Bixi(BSSData,tripsData->second,stationsData->second,start_hour,end_hour,max_stations_allowed);
		if(instanceType->second=="washington")
			load_object.Load_CapitalBikeShare(BSSData,tripsData->second,stationsData->second,start_hour,end_hour,max_stations_allowed);
		if(instanceType->second=="washington" || instanceType->second=="newyork" || instanceType->second=="boston" || instanceType->second=="chicago" || instanceType->second=="sanfrancisco")
					load_object.Load_BSS_Data(BSSData,tripsData->second,stationsData->second,start_hour,end_hour,instanceType->second,max_stations_allowed);	
	}
	
	printf("back in the main ...\n");	
	//getchar();
	
	if(instanceFormat->second == "trips")
		printf("With %d trips and %d days of data then each scenario should have on average %d trips ...\n",(int)BSSData.all_int_trips.size(),BSSData.dayCounter,(int)BSSData.all_int_trips.size()/BSSData.dayCounter);
	if(instanceFormat->second == "requests")
		printf("With %d trips and %d days of data then each scenario should have on average %d pickups and deliveries ...\n",(int)BSSData.all_int_trips.size(),BSSData.dayCounter,(2*(int)BSSData.all_int_trips.size())/BSSData.dayCounter);	

	//printf("Calculating and printing the initial capacities of the stations ...\n");
	//CalculateInitialCapacities(BSSData);
	
	//Uncomment below if you want to plot the trips.
	/*std::string plot_name;
	std::stringstream ss;
	if (instanceType->second == "montreal")
	{
        std::string formattedTripsCount = formatWithThousandsSeparator((int)BSSData.all_int_trips.size());
        std::stringstream ss;
        ss << "Histogram of Bixi Start Times in 04/2023-07/2023: " << formattedTripsCount << " Trips.";
        std::cout << ss.str() << '\n';
        plot_name = ss.str();			
	}		
	if(instanceType->second=="washington")
	{
        std::string formattedTripsCount = formatWithThousandsSeparator((int)BSSData.all_int_trips.size());
        std::stringstream ss;
        ss << "Histogram of Capital Bike Share Start Times in 04/2023-07/2023: " << formattedTripsCount << " Trips.";
        std::cout << ss.str() << '\n';
        plot_name = ss.str();			
	}	
	if(instanceType->second=="boston")
	{
        std::string formattedTripsCount = formatWithThousandsSeparator((int)BSSData.all_int_trips.size());
        std::stringstream ss;
        ss << "Histogram of Blue Bikes Start Times in 04/2023-07/2023: " << formattedTripsCount << " Trips.";
        std::cout << ss.str() << '\n';
        plot_name = ss.str();			
	}	
	if(instanceType->second=="newyork")
	{
        std::string formattedTripsCount = formatWithThousandsSeparator((int)BSSData.all_int_trips.size());
        std::stringstream ss;
        ss << "Histogram of Citi Bike Start Times in 04/2023-07/2023: " << formattedTripsCount << " Trips.";
        std::cout << ss.str() << '\n';
        plot_name = ss.str();			
	}		
	algo.generateHistogramImage(BSSData.all_int_trips,start_hour,end_hour,plot_name);
	exit(1);*/
	
	if(instanceFormat->second == "trips")
	{
		algo.MakeLambdaObject(BSSData.all_int_trips,start_hour,end_hour);
	}
	if(instanceFormat->second == "requests")
	{
		algo.MakePickDelLambdaObject(BSSData.all_int_trips,start_hour,end_hour);
	}
	
	// countBikes.sh calculates the average of bikes in May-June 2025
	if(instanceType->second=="boston") Qtot = 4000; //countBikes.sh: 4497.42 ~ 4500
	if(instanceType->second=="washington") Qtot = 6000; //countBikes.sh: 5125.21 ~ 5100
	if(instanceType->second=="montreal") Qtot = 10000; //countBikes.sh: 10139.68 ~ 10100
	if(instanceType->second=="newyork") Qtot = 32000; //countBikes.sh: 32361.94 ~ 32400
	if(instanceType->second=="chicago") Qtot = 16500; //countBikes.sh: 8642.50 ~ 8600
	if(instanceType->second=="mexicocity") Qtot = 6800; //countBikes.sh: 3418.26 ~ 3400
	if(instanceType->second=="toronto") Qtot = 9000; //countBikes.sh: 5954.26 ~ 6000
	if(instanceType->second=="madrid") Qtot = 7500; //countBikes.sh: 4786.02 ~ 4800
	if(instanceType->second=="sanfrancisco") Qtot = 7000; //countBikes.sh: 7462.26 ~ 7500
	
	int sum_stat_cap=0;
	for(Station s : BSSData.stations)
		sum_stat_cap += s.cap;
	
	printf("Qtot:%d SumStatCap:%d\n",Qtot,sum_stat_cap);
	//printf("exiting ...\n"); exit(1);
	
	if (rebalancing->second == "true") {
		//getchar();
		char filename_stations[100];
		std::map<std::string, std::string> instanceFiles = {
			{"montreal", "../instances/montreal_reb_n%d.txt"},
			{"washington", "../instances/washington_reb_n%d.txt"},
			{"newyork", "../instances/newyork_reb_n%d.txt"},
			{"boston", "../instances/boston_reb_n%d.txt"},
			{"chicago", "../instances/chicago_reb_n%d.txt"},
			{"mexicocity", "../instances/mexicocity_reb_n%d.txt"},
			{"toronto", "../instances/toronto_reb_n%d.txt"},
			{"madrid", "../instances/madrid_reb_n%d.txt"},
			{"sanfrancisco","../instances/sanfrancisco_reb_n%d.txt"}
		};
		auto it = instanceFiles.find(instanceType->second);
		if (it != instanceFiles.end()) {
			snprintf(filename_stations, sizeof(filename_stations), it->second.c_str(), (int)(BSSData.stations.size()));
			writeStationsToFileRebalancing(BSSData.stations, filename_stations, Qtot);
		}
	}

	if (outputFile->second == "split") {
		char filename_stations[100];
		std::map<std::string, std::string> instanceFiles = {
			{"montreal", "../instances/montreal_n%d.txt"},
			{"washington", "../instances/washington_n%d.txt"},
			{"newyork", "../instances/newyork_n%d.txt"},
			{"boston", "../instances/boston_n%d.txt"},
			{"chicago", "../instances/chicago_n%d.txt"},
			{"mexicocity", "../instances/mexicocity_n%d.txt"},
			{"toronto", "../instances/toronto_n%d.txt"},
			{"madrid", "../instances/madrid_n%d.txt"},
			{"sanfrancisco","../instances/sanfrancisco_n%d.txt"}
		};
		auto it = instanceFiles.find(instanceType->second);
		if (it != instanceFiles.end()) {
			snprintf(filename_stations, sizeof(filename_stations), it->second.c_str(), static_cast<int>(BSSData.stations.size()));
			writeStationsToFileTargets(BSSData.stations, filename_stations, Qtot);
		}
	}
	if(outputTrips->second == "random")
	{
		std::vector<LambdaData> & lambdaObject = algo.getLambdaReference();
		std::vector<LambdaData> & lambdaPickObject = algo.getLambdaPickReference();
		std::vector<LambdaData> & lambdaDelObject = algo.getLambdaPickReference();
		//ModifyCounts() is hard coded in the header file. If > 1, then it increases the number of generated random trips or requests (pickups,deliveries)
		if(instanceFormat->second == "trips")
			for(int i=0;i<lambdaObject.size();i++)
				lambdaObject[i].ModifyCounts();
		if(instanceFormat->second == "requests")
		{
			for(int i=0;i<lambdaPickObject.size();i++)
				lambdaPickObject[i].ModifyCounts();
			for(int i=0;i<lambdaDelObject.size();i++)
				lambdaDelObject[i].ModifyCounts();			
		}
		for(int i=0;i<nb_instances;i++)
		{
			//Create a matrix where first index is the scenario and second index is an std::vector of random trips or pickups, deliveries
			//Current version stores the trips, requests in a map
			std::map<int,std::vector<Trip>> trips;
			std::map<int,std::vector<Requests>> picks;
			std::map<int,std::vector<Requests>> dels;
			
			if(algo.AvgScenario==1 && instanceFormat->second == "trips")
				algo.Make_avg_scenario_map(BSSData,trips);
			if(algo.AvgScenario==1 && instanceFormat->second == "requests")
				algo.Make_avg_request_scenario_map(BSSData,picks,dels);
			
			if(instanceFormat->second == "trips" && algo.AvgScenario==0)
			{
				algo.Make_scenarios_map(nb_scenarios,BSSData.nb_days_of_data,start_hour,end_hour,i,trips); 
				printf("Instance:%d Generated scenarios:\n",i);
				for(int e=0;e<nb_scenarios;e++)
					printf("Scenario %d:%d\n",e,(int)trips[e].size());
			}
			if(instanceFormat->second == "requests" && algo.AvgScenario==0)
			{
				algo.Make_PickDelscenarios_map(nb_scenarios,BSSData.nb_days_of_data,start_hour,end_hour,i,picks,dels); 
				for(int e=0;e<nb_scenarios;e++)
					printf("Picks scenario %d:%d\n",e,(int)picks[e].size());
				for(int e=0;e<nb_scenarios;e++)
					printf("Dels scenario %d:%d\n",e,(int)dels[e].size());
			}
			
			char filename[100];
			if(algo.AvgScenario == 1 || ( outputFile->second == "split" && instanceFormat->second == "trips" ))
			{	
				if(algo.AvgScenario == 1)
				{
					if(instanceType->second=="montreal")
						snprintf(filename,sizeof(filename),"../instances/montreal_avg_scenario.txt");
					if(instanceType->second=="washington")
						snprintf(filename,sizeof(filename),"../instances/washington_avg_scenario.txt");
					if(instanceType->second=="boston")
						snprintf(filename,sizeof(filename),"../instances/boston_avg_scenario.txt");
					if(instanceType->second=="newyork")
						snprintf(filename,sizeof(filename),"../instances/newyork_avg_scenario.txt");
					if(instanceType->second=="chicago")	
						snprintf(filename,sizeof(filename),"../instances/chicago_avg_scenario.txt");
					if(instanceType->second=="mexicocity")
						snprintf(filename,sizeof(filename),"../instances/mexicocity_avg_scenario.txt");
					if(instanceType->second=="toronto")
						snprintf(filename,sizeof(filename),"../instances/toronto_avg_scenario.txt");
					if(instanceType->second=="madrid")
						snprintf(filename,sizeof(filename),"../instances/madrid_avg_scenario.txt");
					if(instanceType->second=="sanfrancisco")
						snprintf(filename,sizeof(filename),"../instances/sanfrancisco_avg_scenario.txt");
				}
				else{ 
					if(instanceType->second=="montreal")
						snprintf(filename,sizeof(filename),"../instances/montreal_e%d.txt",nb_scenarios);
					if(instanceType->second=="washington")
						snprintf(filename,sizeof(filename),"../instances/washington_e%d.txt",nb_scenarios);
					if(instanceType->second=="boston")
						snprintf(filename,sizeof(filename),"../instances/boston_e%d.txt",nb_scenarios);
					if(instanceType->second=="newyork")
						snprintf(filename,sizeof(filename),"../instances/newyork_e%d.txt",nb_scenarios);
					if(instanceType->second=="chicago")	
						snprintf(filename,sizeof(filename),"../instances/chicago_e%d.txt",nb_scenarios);
					if(instanceType->second=="mexicocity")
						snprintf(filename,sizeof(filename),"../instances/mexicocity_e%d.txt",nb_scenarios);
					if(instanceType->second=="toronto")
						snprintf(filename,sizeof(filename),"../instances/toronto_e%d.txt",nb_scenarios);
					if(instanceType->second=="madrid")
						snprintf(filename,sizeof(filename),"../instances/madrid_e%d.txt",nb_scenarios);
					if(instanceType->second=="sanfrancisco")
						snprintf(filename,sizeof(filename),"../instances/sanfrancisco_e%d.txt",nb_scenarios);
				}

				printf("Instance:%d Random trips written to %s.\n",i,filename);
				writeTripsToFile(trips, filename);
				//if(algo.AvgScenario == 1) exit(1);
			}
			if(outputFile->second == "split" && instanceFormat->second == "requests")
			{
				if(algo.AvgScenario == 1)
				{
					if(instanceType->second=="montreal")
						snprintf(filename,sizeof(filename),"../instances/montreal_req_avg_scenario.txt");
					if(instanceType->second=="washington")
						snprintf(filename,sizeof(filename),"../instances/washington_req_avg_scenario.txt");
					if(instanceType->second=="boston")
						snprintf(filename,sizeof(filename),"../instances/boston_req_avg_scenario.txt");
					if(instanceType->second=="newyork")
						snprintf(filename,sizeof(filename),"../instances/newyork_req_avg_scenario.txt");
					if(instanceType->second=="chicago")	
						snprintf(filename,sizeof(filename),"../instances/chicago_req_avg_scenario.txt");
					if(instanceType->second=="mexicocity")
						snprintf(filename,sizeof(filename),"../instances/mexicocity_req_avg_scenario.txt");
					if(instanceType->second=="toronto")
						snprintf(filename,sizeof(filename),"../instances/toronto_req_avg_scenario.txt");
					if(instanceType->second=="madrid")
						snprintf(filename,sizeof(filename),"../instances/madrid_req_avg_scenario.txt");
					if(instanceType->second=="sanfrancisco")
						snprintf(filename,sizeof(filename),"../instances/sanfrancisco_req_avg_scenario.txt");					
				}
				else 
				{
					if(instanceType->second=="montreal")
						snprintf(filename,sizeof(filename),"../instances/montreal_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="washington")
						snprintf(filename,sizeof(filename),"../instances/washington_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="boston")
						snprintf(filename,sizeof(filename),"../instances/boston_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="newyork")
						snprintf(filename,sizeof(filename),"../instances/newyork_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="chicago")	
						snprintf(filename,sizeof(filename),"../instances/chicago_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="mexicocity")
						snprintf(filename,sizeof(filename),"../instances/mexicocity_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="toronto")
						snprintf(filename,sizeof(filename),"../instances/toronto_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="madrid")
						snprintf(filename,sizeof(filename),"../instances/madrid_req_e%d.txt",nb_scenarios);
					if(instanceType->second=="sanfrancisco")
						snprintf(filename,sizeof(filename),"../instances/sanfrancisco_req_e%d.txt",nb_scenarios);
				}

				printf("Instance:%d Random requests written to %s.\n",i,filename);
				writeRequestsToFile(picks,dels, filename);
			}			
			if(outputFile->second == "single" && algo.AvgScenario == 0)
			{
				if(instanceType->second=="montreal" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/montreal_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="washington" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/washington_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="boston" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/boston_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="newyork" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/newyork_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="chicago" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/chicago_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="mexicocity" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/mexicocity_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="toronto" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/toronto_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="madrid" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/madrid_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="sanfrancisco" && instanceFormat->second == "trips")
					snprintf(filename,sizeof(filename),"../instances/sanfrancisco_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);

				if(instanceType->second=="montreal" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/montreal_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="washington" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/washington_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="boston" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/boston_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="newyork" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/newyork_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="chicago" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/chicago_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="mexicocity" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/mexicocity_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="toronto" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/toronto_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="madrid" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/madrid_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);
				if(instanceType->second=="sanfrancisco" && instanceFormat->second == "requests")
					snprintf(filename,sizeof(filename),"../instances/sanfrancisco_req_n%d_e%d_%d.txt",(int)BSSData.stations.size(),nb_scenarios,i);				
				
				printf("Printing instance:%d\n",i);
				if(instanceFormat->second == "trips")
					WriteTripsAndStationsToFile(trips, BSSData.stations,Qtot,filename);
				if(instanceFormat->second == "requests")
					writeRequestsAndStationsToFile(picks,dels,BSSData.stations,Qtot,filename);
				
				printf("Single output_file for instance %d writtent to:%s.\n",i,filename);
				
			}
			
		}
		
	}

	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Total generator time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;
	
	return 0;
}


