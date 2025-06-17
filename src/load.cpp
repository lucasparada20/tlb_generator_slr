#include "load.h"
#include "json.hpp"
#include "csv.h"

bool compare_tm(const std::tm& a, const std::tm& b) {
	std::time_t timestamp_a = std::mktime(const_cast<std::tm*>(&a));
	std::time_t timestamp_b = std::mktime(const_cast<std::tm*>(&b));
	return std::difftime(timestamp_a, timestamp_b) < 0;
}

// Function to convert TripDataAllInts to std::tm
std::tm to_tm(const TripDataAllInts& trip) {
    std::tm date = {0};
    date.tm_year = trip.start_year - 1900;
    date.tm_mon = trip.start_month - 1;
    date.tm_mday = trip.start_day;
    return date;
}

void extractIdentifier(const std::string& input, std::string & identifier) {
    for (char c : input) {
        if (c == '-') {
            break; // Stop when hyphen is encountered
        }
        if (isdigit(c) || isalpha(c)) {
            identifier += c; // Append digits and letters to identifier
        }
    }
}

void Load::Load_Madrid(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour)
{
	printf("In Load_Madrid reading %s and %s\n", filename_csv.c_str(), filename_json.c_str());
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
    clock_t lastPrintTime = startTime;
	
	// Open files
    std::ifstream file_csv(filename_csv); std::ifstream file_json(filename_json);
    if (!file_csv.is_open() || !file_json.is_open()) {
        std::cerr << "Error: Could not open files " << filename_csv << " " << filename_json << ". Exiting ..." <<std::endl;
        exit(1);
    }
	
	// Parse the JSON data; 
	nlohmann::json jsonData;
	file_json >> jsonData;
	nlohmann::json json_stations = jsonData["data"]["stations"];
	
	station_name_map.clear();
	station_lat_map.clear();
	station_lon_map.clear();
	Station depot; depot.name = std::string("depot"); depot.id = 0; depot.cap = 0;
	//HQ of Empresa Municipal de Transportes de Madrid, operator of BiciMad
	depot.lat = 40.39705;
	depot.lon = -3.67147;
	data.stations.push_back(depot);
	int station_counter=1;
	std::map<std::string,int> bss_id_map;
	for (const auto& json_station : json_stations)
	{
		//std::cout << json_station << std::endl;
		Station station;
		std::string station_name_str = json_station["name"]; station.name = station_name_str;
	
		station.id = station_counter;
		
		//Here, annoyingly, the stations have an identifier inside the string of the name ....
		std::string identifier;
		extractIdentifier(station_name_str,identifier);
		bss_id_map[ identifier ] = station_counter;
		
		station.lon = json_station["lon"]; station.lat = json_station["lat"]; station.cap = json_station["capacity"];
		//station.Show();
		//getchar();
		station_name_map.insert(std::pair<std::string,Station>(station.name,station));
		station_lat_map.insert(std::pair<double,Station>(station.lat,station));
		station_lon_map.insert(std::pair<double,Station>(station.lon,station));
		
		//For the rebalancing instances
		station.bss_station_id = json_station["station_id"];
		
		data.stations.push_back(station);
		station_counter++;
	}
	file_json.close();
	printf("Stations:%d ObjectSize:%lu Sample station object:\n",(int)data.stations.size(),sizeof(data.stations[0])); data.stations[0].Show(); //getchar();

	int line_counter=0; int exception_counter=0; int no_id_counter =0; int out_of_time =0; int empty_line =0; int no_time_difference;
	std::string line;
	std::getline(file_csv, line); //the headers of the columns
	
    while (std::getline(file_csv, line))
	{		
		line_counter++;
		if (line_counter % 500000 == 0)
			std::cout << "Read " << line_counter << " lines from .csv ..." << std::endl;
		
		//std::cout << line << std::endl;
		
		std::vector<std::string> tokens;
		size_t startPos = 0;
		// Find the positions of commas and extract substrings
		for (size_t pos = line.find(';'); pos != std::string::npos; pos = line.find(';', startPos)) {
			tokens.push_back(line.substr(startPos, pos - startPos));
			startPos = pos + 1;
		}
		// Add the last token after the final comma
		tokens.push_back(line.substr(startPos));
		
		if(tokens.size() == 19 || tokens.size() == 18)
		{			
				
			//for(int i=0;i<tokens.size();i++)
			//	printf("i:%d tok:%s ",i,tokens[i].c_str());
			//printf("\n");
			//getchar();
			
			int start_id = -1; int end_id = -1;
			std::string unlock_station_name = tokens.size() == 19 ? tokens[15] : tokens[14];
			std::string lock_station_name = tokens.size() == 19 ? tokens[18] : tokens[17];
			
			std::string start_identifier; std::string end_identifier;
			extractIdentifier(unlock_station_name,start_identifier);
			extractIdentifier(lock_station_name,end_identifier);			

			auto it_start = bss_id_map.find(start_identifier);
			auto it_end = bss_id_map.find(end_identifier);			

			if( it_start != bss_id_map.end())
				start_id =  it_start->second;
			if( it_end != bss_id_map.end())
				end_id =  it_end->second;
			
			//std::cout << line << std::endl << std::endl;			
			//std::cout << " start_id: " << start_id << " startName: " << unlock_station_name << " end_id: " << end_id << " endName: " << lock_station_name << std::endl;
			//getchar();
			
			if(start_id==-1 || end_id ==-1)
			{
				//std::cout << line << std::endl << std::endl;			
				//std::cout << " start_id: " << start_id << " startName: " << unlock_station_name << " end_id: " << end_id << " endName: " << lock_station_name << std::endl;
				//getchar();
				no_id_counter++; continue;
			}
			
			TripDataAllInts all_ints_trip;
			all_ints_trip.start_id = start_id;
			all_ints_trip.end_id = end_id;
					
			//Start time is token[7]
			sscanf(tokens.size() == 19 ? tokens[7].c_str() : tokens[6].c_str(),
				"%hd-%hhd-%hhdT%hhd:%hhd:%hhd",
				&all_ints_trip.start_year,&all_ints_trip.start_month,&all_ints_trip.start_day,&all_ints_trip.start_hour,&all_ints_trip.start_minute,&all_ints_trip.start_second);
			
			//End time is token[12]
			sscanf(tokens.size() == 19 ? tokens[12].c_str() : tokens[11].c_str(),
				"%hd-%hhd-%hhdT%hhd:%hhd:%hhd",
				&all_ints_trip.end_year,&all_ints_trip.end_month,&all_ints_trip.end_day,&all_ints_trip.end_hour,&all_ints_trip.end_minute,&all_ints_trip.end_second);
			
			if (all_ints_trip.start_hour < start_hour || all_ints_trip.start_hour > end_hour-1 || all_ints_trip.end_hour < start_hour || all_ints_trip.end_hour > end_hour -1) 
			{
				out_of_time++; continue;
			}
			
			//all_ints_trip.Show(); getchar();
			if(all_ints_trip.start_month<7 || all_ints_trip.end_month>8 || all_ints_trip.end_year>2023 || all_ints_trip.start_day > all_ints_trip.end_day)
			{
				all_ints_trip.Show();
				//getchar();
				out_of_time++;
				continue;
			}			

			int duration_seconds = 0;
			int end_total_seconds = all_ints_trip.end_hour * 3600 + all_ints_trip.end_minute * 60 + all_ints_trip.end_second;
			int start_total_seconds = all_ints_trip.start_hour * 3600 + all_ints_trip.start_minute * 60 + all_ints_trip.start_second;
			
			duration_seconds = end_total_seconds - start_total_seconds;
			if(duration_seconds < 60) 
			{	
				no_time_difference++; continue;
			}
			
			all_ints_trip.duration = duration_seconds;
			if(duration_seconds == -1)
			{
				printf("start:%d end:%d duration:%d\n",start_total_seconds,end_total_seconds,duration_seconds); getchar();
			}
				
			data.all_int_trips.push_back(all_ints_trip);

		}
		else {
			printf("failed to parse csv line:%s\n",line.c_str()); //getchar();
		}
	}//End while
	
	if(line_counter != out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size())
	{
		int sum = out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size();
		printf("line_counter:%d sum:%d out_of_time:%d no_time_difference:%d no_id_counter:%d loaded_trips:%d\n",line_counter,sum,out_of_time,no_time_difference,no_id_counter,(int)data.all_int_trips.size()); getchar();
	}
	
	
	file_csv.close();
	printf("Trips:%d ObjectSize:%lu Sample station object:\n",(int)data.all_int_trips.size(),sizeof(data.all_int_trips[0])); data.all_int_trips[0].Show();

	printf("The .csv has %d rows ...\n",line_counter);
	printf("There were %d empty lines ...\n",empty_line);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... Continuing ... \n",out_of_time,start_hour,end_hour);
	printf("There were %d trips with no time difference.\n",no_time_difference);
	printf("Sum of missing data for exel:%d\n",no_id_counter+no_time_difference);
	
	printf("Loaded %d trips and %d stations.\n",(int)data.all_int_trips.size(),(int)data.stations.size());
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;

	auto minmax_dates = std::minmax_element(data.all_int_trips.begin(), data.all_int_trips.end(), [this](const TripDataAllInts& a, const TripDataAllInts& b) {
		return compare_tm(to_tm(a), to_tm(b));
	});

    if (minmax_dates.first != data.all_int_trips.end() && minmax_dates.second != data.all_int_trips.end()) {
        // Convert start date of earliest trip to std::tm
        std::tm earliest_date = to_tm(*minmax_dates.first);

        // Convert end date of latest trip to std::tm
        std::tm latest_date = to_tm(*minmax_dates.second);

        // Calculate the difference in days
        std::time_t earliest_time = std::mktime(&earliest_date);
        std::time_t latest_time = std::mktime(&latest_date);
        int difference_in_days = std::difftime(latest_time, earliest_time) / (60 * 60 * 24);
		
		// Print the earliest and latest dates
		std::cout << "Earliest date: " << std::put_time(&earliest_date, "%c") << std::endl;
		std::cout << "Latest date: " << std::put_time(&latest_date, "%c") << std::endl;

        // Print the difference in number of days
        std::cout << "Difference between the start date of the earliest trip and the end date of the latest trip in terms of days:" << difference_in_days << " dayCounter:" <<  difference_in_days+1 << std::endl;
		data.nb_days_of_data = difference_in_days;
		data.dayCounter = difference_in_days+1;
    } else {
        std::cout << "No trips recorded." << std::endl;
    }		

	
}

void Load::Load_Toronto(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour)
{
	printf("In Load_Toronto reading %s and %s\n", filename_csv.c_str(), filename_json.c_str());
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
    clock_t lastPrintTime = startTime;
	
	// Open files
    std::ifstream file_csv(filename_csv); std::ifstream file_json(filename_json);
    if (!file_csv.is_open() || !file_json.is_open()) {
        std::cerr << "Error: Could not open files " << filename_csv << " " << filename_json << ". Exiting ..." <<std::endl;
        exit(1);
    }
	
	// Parse the JSON data; 
	nlohmann::json jsonData;
	file_json >> jsonData;
	nlohmann::json json_stations = jsonData["data"]["stations"];
	
	station_name_map.clear();
	station_lat_map.clear();
	station_lon_map.clear();
	Station depot; depot.name = std::string("depot"); depot.id = 0; depot.cap = 0;
	//HQ of Toronto Parking Authority, operator of Toronto Bike Share
	depot.lat = 43.65282;
	depot.lon = -79.37729;
	data.stations.push_back(depot);
	
	std::map<int,int> bss_id_map;
	int station_counter=1;
	for (const auto& json_station : json_stations)
	{
		//std::cout << json_station << std::endl;
		Station station;
		std::string station_name_str = json_station["name"]; station.name = station_name_str;
		
		station.id = station_counter;
		
		std::string station_id_str = json_station["station_id"];
		int bss_id = std::stoi( station_id_str );
		
		bss_id_map[ bss_id ] = station_counter;
		
		station.lon = json_station["lon"]; station.lat = json_station["lat"]; station.cap = json_station["capacity"];
		//station.Show();
		//getchar();
		station_name_map.insert(std::pair<std::string,Station>(station.name,station));
		station_lat_map.insert(std::pair<double,Station>(station.lat,station));
		station_lon_map.insert(std::pair<double,Station>(station.lon,station));
		
		//For the rebalancing instances
		station.bss_station_id = json_station["station_id"];
		
		data.stations.push_back(station);
		station_counter++;
	}
	
	//for(auto it=bss_id_map.begin();it!=bss_id_map.end();++it)
		//printf("BssId:%d Id:%d\n",it->first,it->second);
	
	file_json.close();
	printf("Stations:%d ObjectSize:%lu Sample station object:\n",(int)data.stations.size(),sizeof(data.stations[0])); data.stations[0].Show();

	int line_counter=0; int exception_counter=0; int no_id_counter =0; int out_of_time =0; int no_time_difference=0;
	std::string line;
	std::getline(file_csv, line); //the headers of the columns

    while (std::getline(file_csv, line))
	{		
		line_counter++;
		if (line_counter % 500000 == 0)
			std::cout << "Read " << line_counter << " lines from .csv ..." << std::endl;
		
		//std::cout << line << std::endl;
		
		std::vector<std::string> tokens;
		size_t startPos = 0;
		// Find the positions of commas and extract substrings
		for (size_t pos = line.find(','); pos != std::string::npos; pos = line.find(',', startPos)) {
			tokens.push_back(line.substr(startPos, pos - startPos));
			startPos = pos + 1;
		}
		// Add the last token after the final comma
		tokens.push_back(line.substr(startPos));
		
		if(tokens.size() == 10)
		{
			//for(int i=0;i<tokens.size();i++)
				//printf("i:%d tok:%s ",i,tokens[i].c_str());
			//printf("\n");
			
			TripDataAllInts all_ints_trip;
			
			auto it_start_id = bss_id_map.find(std::atoi(tokens[2].c_str()));
			auto it_end_id = bss_id_map.find(std::atoi(tokens[5].c_str()));
			
			if(it_start_id==bss_id_map.end() && it_end_id==bss_id_map.end())
			{
				no_id_counter++; continue;
			}
			
			if(it_start_id->second == 0 || it_end_id->second == 0)
			{
				no_id_counter++; continue;
				//std::cout << "Tokens[2]:" << tokens[2] << " Tokens[5]:" << tokens[5] << std::endl;
				//std::cout << "Tokens[2]Int:" << std::atoi(tokens[2].c_str()) << " Tokens[5]Int:" << std::atoi(tokens[5].c_str()) << std::endl;
				//printf("Start BssId:%d Id:%d\n",it_start_id->first,it_start_id->second);
				//printf("End BssId:%d Id:%d\n",it_end_id->first,it_end_id->second);
				//std::cout << line << std::endl;
				//getchar();
			}
			
			
			all_ints_trip.start_id = it_start_id->second;
			all_ints_trip.end_id = it_end_id->second;
			
			all_ints_trip.duration = std::stoi( tokens[1] );
			if( all_ints_trip.duration < 60 ) 
			{
				no_time_difference++; continue;
			}
			
			sscanf(tokens[3].c_str(), "%hhd/%hhd/%hd %hhd:%hhd", &all_ints_trip.start_month, &all_ints_trip.start_day, &all_ints_trip.start_year,&all_ints_trip.start_hour,&all_ints_trip.start_minute);
			
			sscanf(tokens[6].c_str(), "%hhd/%hhd/%hd %hhd:%hhd", &all_ints_trip.end_month, &all_ints_trip.end_day, &all_ints_trip.end_year,&all_ints_trip.end_hour,&all_ints_trip.end_minute);
			
			if (all_ints_trip.start_hour < start_hour || all_ints_trip.start_hour > end_hour-1 || all_ints_trip.end_hour < start_hour || all_ints_trip.end_hour > end_hour -1) 
			{
				out_of_time++; continue;
			}

			//all_ints_trip.Show(); 
				
			if(all_ints_trip.start_month<7 || all_ints_trip.end_month>7 || all_ints_trip.start_year<2023 || all_ints_trip.end_year>2023 || all_ints_trip.start_day > all_ints_trip.end_day)
			{
				//all_ints_trip.Show();
				//getchar();
				out_of_time++;
				continue;
			}
				
			data.all_int_trips.push_back(all_ints_trip);

		}
		else {
			printf("failed to parse csv line:%s\n",line.c_str());
		}
	}//End while
	file_csv.close();
	
	if(line_counter != out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size())
	{
		int sum = out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size();
		printf("line_counter:%d sum:%d out_of_time:%d no_time_difference:%d no_id_counter:%d loaded_trips:%d\n",line_counter,sum,out_of_time,no_time_difference,no_id_counter,(int)data.all_int_trips.size()); getchar();
	}
	
	printf("Trips:%d ObjectSize:%lu Sample station object:\n",(int)data.all_int_trips.size(),sizeof(data.all_int_trips[0])); data.all_int_trips[0].Show();

	printf("The .csv has %d rows ...\n",(int)data.all_int_trips.size()+no_id_counter+out_of_time);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d trips with no time difference.\n",no_time_difference);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... Continuing ... \n",out_of_time,start_hour,end_hour);
	printf("Total trips with missing data for exel:%d\n",no_time_difference+no_id_counter);
	
	printf("Loaded %d trips and %d stations.\n",(int)data.all_int_trips.size(),(int)data.stations.size());
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;

	auto minmax_dates = std::minmax_element(data.all_int_trips.begin(), data.all_int_trips.end(), [this](const TripDataAllInts& a, const TripDataAllInts& b) {
		return compare_tm(to_tm(a), to_tm(b));
	});

    if (minmax_dates.first != data.all_int_trips.end() && minmax_dates.second != data.all_int_trips.end()) {
        // Convert start date of earliest trip to std::tm
        std::tm earliest_date = to_tm(*minmax_dates.first);

        // Convert end date of latest trip to std::tm
        std::tm latest_date = to_tm(*minmax_dates.second);

        // Calculate the difference in days
        std::time_t earliest_time = std::mktime(&earliest_date);
        std::time_t latest_time = std::mktime(&latest_date);
        int difference_in_days = std::difftime(latest_time, earliest_time) / (60 * 60 * 24);
		
		// Print the earliest and latest dates
		std::cout << "Earliest date: " << std::put_time(&earliest_date, "%c") << std::endl;
		std::cout << "Latest date: " << std::put_time(&latest_date, "%c") << std::endl;

        // Print the difference in number of days
        std::cout << "Difference between the start date of the earliest trip and the end date of the latest trip in terms of days:" << difference_in_days << " dayCounter:" <<  difference_in_days+1 << std::endl;
		data.nb_days_of_data = difference_in_days;
		data.dayCounter = difference_in_days+1;
    } else {
        std::cout << "No trips recorded." << std::endl;
    }		
		
	
}

void Load::Load_Ecobici(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour)
{
	printf("In Load_Ecobici reading %s and %s\n", filename_csv.c_str(), filename_json.c_str());
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
    clock_t lastPrintTime = startTime;
	
	// Open files
    std::ifstream file_csv(filename_csv); std::ifstream file_json(filename_json);
    if (!file_csv.is_open() || !file_json.is_open()) {
        std::cerr << "Error: Could not open files " << filename_csv << " " << filename_json << ". Exiting ..." <<std::endl;
        exit(1);
    }
	
	// Parse the JSON data; 
	nlohmann::json jsonData;
	file_json >> jsonData;
	nlohmann::json json_stations = jsonData["data"]["stations"];
	
	station_name_map.clear();
	station_lat_map.clear();
	station_lon_map.clear();
	Station depot; depot.name = std::string("depot"); depot.id = 0; depot.cap = 0;
	//Oficinas ecobici mx df
	depot.lat = 19.42997;
	depot.lon = -99.17846;
	data.stations.push_back(depot);
	int station_counter=1;
	for (const auto& json_station : json_stations)
	{
		//std::cout << json_station << std::endl;
		Station station;
		std::string station_name_str = json_station["name"]; station.name = station_name_str;
		//std::string station_id_str = json_station["station_id"]; station.id = atoi(station_id_str.c_str());
		station.id = station_counter;
		station.lon = json_station["lon"];
		station.lat = json_station["lat"];
		station.cap = json_station["capacity"];
		//station.Show();
		//getchar();
		station_name_map.insert(std::pair<std::string,Station>(station.name,station));
		station_lat_map.insert(std::pair<double,Station>(station.lat,station));
		station_lon_map.insert(std::pair<double,Station>(station.lon,station));
		
		//For the rebalancing instances
		station.bss_station_id = json_station["station_id"];
		
		data.stations.push_back(station);
		station_counter++;
	}
	file_json.close();	
	printf("Stations:%d ObjectSize:%lu Sample station object:\n",(int)data.stations.size(),sizeof(data.stations[0])); data.stations[0].Show();

	int line_counter=0; int exception_counter=0; int no_id_counter =0; int out_of_time =0; int no_time_difference=0;
	std::string line;
	std::getline(file_csv, line); //the headers of the columns

    while (std::getline(file_csv, line))
	{		
		line_counter++;
		if (line_counter % 500000 == 0)
			std::cout << "Read " << line_counter << " lines from .csv ..." << std::endl;
		
		std::vector<std::string> tokens;
		size_t startPos = 0;
		// Find the positions of commas and extract substrings
		for (size_t pos = line.find(','); pos != std::string::npos; pos = line.find(',', startPos)) {
			tokens.push_back(line.substr(startPos, pos - startPos));
			startPos = pos + 1;
		}
		// Add the last token after the final comma
		tokens.push_back(line.substr(startPos));
		
		if(tokens.size() == 9)
		{
			//for(int i=0;i<tokens.size();i++)
				//printf("i:%d tok:%s ",i,tokens[i].c_str());
			//printf("\n");
			
			int start_id = -1; int end_id = -1;
			try {
				start_id = std::stoi(tokens[3]);
				end_id = std::stoi(tokens[6]);
			} catch (const std::invalid_argument& e) {
				no_id_counter++;
				continue; // Move to the next iteration of the loop
			}
			
			if(start_id > data.stations.size() || end_id > data.stations.size() || start_id < 0 || end_id < 0)
			{
				no_id_counter++; continue;
			}
			
			TripDataAllInts all_ints_trip;
			sscanf(tokens[4].c_str(), "%hhd/%hhd/%hd", &all_ints_trip.start_day, &all_ints_trip.start_month, &all_ints_trip.start_year);
			sscanf(tokens[7].c_str(), "%hhd/%hhd/%hd", &all_ints_trip.end_day, &all_ints_trip.end_month, &all_ints_trip.end_year);
			sscanf(tokens[5].c_str(), "%hhd:%hhd:%hhd", &all_ints_trip.start_hour, &all_ints_trip.start_minute, &all_ints_trip.start_second);
			sscanf(tokens[8].c_str(), "%hhd:%hhd:%hhd", &all_ints_trip.end_hour, &all_ints_trip.end_minute, &all_ints_trip.end_second);
			
			if (all_ints_trip.start_hour < start_hour || all_ints_trip.start_hour > end_hour-1 || all_ints_trip.end_hour < start_hour || all_ints_trip.end_hour > end_hour -1) 
			{
				out_of_time++; continue;
			}

			double duration_seconds = 0;
			int end_total_seconds = all_ints_trip.end_hour * 3600 + all_ints_trip.end_minute * 60 + all_ints_trip.end_second;
			int start_total_seconds = all_ints_trip.start_hour * 3600 + all_ints_trip.start_minute * 60 + all_ints_trip.start_second;
			
			duration_seconds = end_total_seconds - start_total_seconds;
			if(duration_seconds < 60) 
			{
				no_time_difference++; continue;
			}	
			
			all_ints_trip.duration = duration_seconds;
			all_ints_trip.start_id = start_id;
			all_ints_trip.end_id = end_id;
			//all_ints_trip.Show();
			
			if(all_ints_trip.start_month<7 || all_ints_trip.end_month>7 || all_ints_trip.start_day > all_ints_trip.end_day)
			{
				//all_ints_trip.Show();
				//getchar();
				out_of_time++;
				continue;
			}
				
			data.all_int_trips.push_back(all_ints_trip);

		}
		else {
			printf("failed to parse csv line:%s\n",line.c_str());
		}
	}//End while
	file_csv.close();
	printf("Trips:%d ObjectSize:%lu Sample station object:\n",(int)data.all_int_trips.size(),sizeof(data.all_int_trips[0])); data.all_int_trips[0].Show();

	
	if(line_counter != out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size())
	{
		int sum = out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size();
		printf("line_counter:%d sum:%d out_of_time:%d no_time_difference:%d no_id_counter:%d loaded_trips:%d\n",line_counter,sum,out_of_time,no_time_difference,no_id_counter,(int)data.all_int_trips.size()); getchar();
	}


	printf("The .csv has %d rows ...\n",(int)data.all_int_trips.size()+exception_counter+no_id_counter+out_of_time);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... Continuing ... \n",out_of_time,start_hour,end_hour);
	printf("There were %d trips with no time difference.\n",no_time_difference);
	printf("Total trips with missing data for exel:%d\n",no_time_difference+no_id_counter);
	
	printf("Loaded %d trips and %d stations.\n",(int)data.all_int_trips.size(),(int)data.stations.size());
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;

	auto minmax_dates = std::minmax_element(data.all_int_trips.begin(), data.all_int_trips.end(), [this](const TripDataAllInts& a, const TripDataAllInts& b) {
		return compare_tm(to_tm(a), to_tm(b));
	});

    if (minmax_dates.first != data.all_int_trips.end() && minmax_dates.second != data.all_int_trips.end()) {
        // Convert start date of earliest trip to std::tm
        std::tm earliest_date = to_tm(*minmax_dates.first);

        // Convert end date of latest trip to std::tm
        std::tm latest_date = to_tm(*minmax_dates.second);

        // Calculate the difference in days
        std::time_t earliest_time = std::mktime(&earliest_date);
        std::time_t latest_time = std::mktime(&latest_date);
        int difference_in_days = std::difftime(latest_time, earliest_time) / (60 * 60 * 24);
		
		// Print the earliest and latest dates
		std::cout << "Earliest date: " << std::put_time(&earliest_date, "%c") << std::endl;
		std::cout << "Latest date: " << std::put_time(&latest_date, "%c") << std::endl;

        // Print the difference in number of days
        std::cout << "Difference between the start date of the earliest trip and the end date of the latest trip in terms of days:" << difference_in_days << " dayCounter:" <<  difference_in_days+1 << std::endl;
		data.nb_days_of_data = difference_in_days;
		data.dayCounter = difference_in_days+1;
    } else {
        std::cout << "No trips recorded." << std::endl;
    }		
		
}

void Load::Load_BSS_Data(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour, std::string type, int max_stations_allowed)
{

	printf("In Load_BSS_Data reading %s and %s of type %s...\n", filename_csv.c_str(), filename_json.c_str(),type.c_str());
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
    clock_t lastPrintTime = startTime;
	
	// Open files
    std::ifstream file_csv(filename_csv); std::ifstream file_json(filename_json);
    if (!file_csv.is_open() || !file_json.is_open()) {
        std::cerr << "Error: Could not open files " << filename_csv << " " << filename_json << ". Exiting ..." <<std::endl;
        exit(1);
    }
	
	// Parse the JSON data; 
	nlohmann::json jsonData;
	file_json >> jsonData;
	nlohmann::json json_stations = jsonData["data"]["stations"];
	
	station_name_map.clear();
	station_lat_map.clear();
	station_lon_map.clear();
	
	Station depot; depot.name = std::string("depot"); depot.id = 0; depot.cap = 0;
	if(type=="boston")
	{
		// Blue Bikes headquarters defined at 100 Northern Avenue, Boston, MA 02210
		//42.3517° N, 71.0405° W. 
		depot.lon = -71.0405; depot.lat = 42.3517; 
	}
	if(type=="washington")
	{
		// Capital Bikeshare's headquarters are located at 1501 Wilson Blvd Ste 1100, Arlington, Virginia, 22209
		//38.89538947941347, -77.07424659999998	
		depot.lon = -77.0742465; depot.lat = 38.895389; depot.cap = 0;
	}
	if(type=="newyork")
	{
		// Citi Bike's headquarters are located at Brooklyn, 3A 220 36th St, United States
		//40.65686421320091, -74.00812771723888	
		depot.lon = -74.00812771723888; depot.lat = 40.65686421320091; 
	}
	if(type=="chicago")
	{
		// Divvy Service Warehouse : 41.89695600489008, -87.67771845224394
		depot.lat = 41.89695600489008; depot.lon = -87.67771845224394;
	}
	if(type=="sanfrancisco")
	{
		//Metropolitan Transportation Commission in SF.
		depot.lat = 37.78826; depot.lon = -122.39155;
	}
	data.stations.push_back(depot);
	
	bss_id_shortName_map.clear();
	
	int station_counter=1;
	for (const auto& json_station : json_stations)
	{
		//std::cout << json_station << std::endl;
		Station station;
		std::string station_name_str = json_station["name"];
		station.name = station_name_str;
		//std::string station_id_str = json_station["station_id"]; station.id = atoi(station_id_str.c_str());
		station.id = station_counter;
		
		// Check if "short_name" exists in the JSON object before accessing it
		if (json_station.contains("short_name")) {
			bss_id_shortName_map[json_station["short_name"]] = station_counter;
		} else {
			// Handle case where "short_name" doesn't exist
			// You can skip this station or handle it differently based on your requirements
			// For example:
			std::cerr << "Warning: 'short_name' key not found in JSON station" << std::endl;
			//continue; // Skip this station and move on to the next one
		}
		
		station.lon = json_station["lon"];
		station.lat = json_station["lat"];
		station.cap = json_station["capacity"];
		if(max_stations_allowed > 0)
			station.Show();
		//getchar();
		station_name_map.insert(std::pair<std::string,Station>(station.name,station));
		station_lat_map.insert(std::pair<double,Station>(station.lat,station));
		station_lon_map.insert(std::pair<double,Station>(station.lon,station));
		
		//For the rebalancing instances
		station.bss_station_id = json_station["station_id"];
		
		if(max_stations_allowed > 0 && station_counter == max_stations_allowed) break;
		
		data.stations.push_back(station);
		station_counter++;
	}

	file_json.close();
	printf("Stations:%d Size:%lu Sample station object:\n",(int)data.stations.size(),sizeof(data.stations[0])); data.stations[0].Show();
	//exit(1);

	Load_CSV_File(data, filename_csv, start_hour, end_hour);
	
	printf("Loaded %d trips and %d stations.\n",(int)data.all_int_trips.size(),(int)data.stations.size());
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;	
	
	//for(int i=0;i<data.trips.size();i++)
	//{ printf("i:%d\n",i); data.trips[i].Show(); }

}

void Load::Load_CSV_File(BSS_Data & data, std::string filename_csv, int start_hour, int end_hour)
{
	printf("Loading using csv.h library ...\n");
	
	//to count the data days ... Initialize variables to track earliest and latest dates
    std::tm earliestDate = {0}; TripData earliestTrip;
    std::tm latestDate = {0}; TripData latestTrip;
    bool firstDate = true;
	data.dayCounter = 0; data.nb_days_of_data = 0;
	
	// Open files
    std::ifstream file_csv(filename_csv); 
    if (!file_csv.is_open()) {
        std::cerr << "Error: Could not open csv file " << filename_csv << " " << ". Exiting ..." <<std::endl;
        exit(1);
    }
	//std::string line;
	//std::getline(file_csv, line); //the headers of the columns
	//std::cout << line <<std::endl;
	//std::getline(file_csv, line); //the next line 
	//std::cout << line <<std::endl;
	//file_csv.close();
	
	std::vector<TripData> trips;
	int line_counter=0; int exception_counter=0; int no_id_counter =0; int out_of_time =0; int no_time_difference=0; int trip_with_missing_data=0; int next_day_counter = 0;
	
	int nb_found_by_name, nb_found_from_coord = 0;
	
	try{		
		//io::CSVReader<13> reader(filename_csv); // Create a CSV reader with 13 columns
		io::CSVReader<13, io::trim_chars<>, io::double_quote_escape<',','\"'>, io::throw_on_overflow, io::no_comment> reader(filename_csv);

		reader.read_header(io::ignore_extra_column,
						   "ride_id","rideable_type","started_at","ended_at","start_station_name","start_station_id",
						   "end_station_name","end_station_id","start_lat","start_lng","end_lat","end_lng","member_casual");

		std::string ride_id, rideable_type, started_at, ended_at, start_station_name, end_station_name, member_casual;
		std::string start_station_id_str, end_station_id_str; 
		std::string start_lat_str; std::string start_lng_str; std::string end_lat_str; std::string end_lng_str;
		
		while (reader.read_row(ride_id,rideable_type,started_at,ended_at,start_station_name,start_station_id_str,end_station_name,end_station_id_str,start_lat_str,start_lng_str,end_lat_str,end_lng_str,member_casual)) 
		{
			line_counter++;

			//std::cout << "Row: ";
				//std::cout << ride_id << ", " << rideable_type << ", " << started_at << ", " << ended_at << ", "
              //<< start_station_name << ", " << start_station_id_str << ", " << end_station_name << ", "
              //<< end_station_id_str << ", " << start_lat_str << ", " << start_lng_str << ", "
              //<< end_lat_str << ", " << end_lng_str << ", " << member_casual << std::endl;
			//getchar();
			if (line_counter % 500000 == 0) 
				std::cout << "Read " << line_counter << " lines from .csv ..." << std::endl;
			
			// Check if any of the required columns are empty
			if (started_at.empty() || ended_at.empty() || start_station_name.empty() || end_station_name.empty() || start_lat_str.empty() || start_lng_str.empty() || end_lat_str.empty() || end_lng_str.empty()) {
				trip_with_missing_data++; continue; // Skip this row if any required column is empty
			}

			double start_lat_double; double start_lon_double; double end_lat_double; double end_lon_double;
			try
			{
				start_lat_double = std::stod(start_lat_str); start_lon_double = std::stod(start_lng_str);
				end_lat_double = std::stod(end_lat_str); end_lon_double = std::stod(end_lng_str);				
			} catch (const std::exception& e)
			{
				exception_counter++;
				//std::cout << line << std::endl;
				//for(int i=0; i<tokens.size(); i++)
					//std::cout << tokens[i] << " ";
				//printf("\n");
				//getchar();
				//std::cerr << "Exception encountered: " << e.what() << std::endl;
				continue;
			}
						
			std::string start_date_str, start_time_str;
			std::string end_date_str, end_time_str;

			// Find the position of the space character separating date and time
			size_t spacePos = started_at.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				start_date_str = started_at.substr(0, spacePos);
				start_time_str = started_at.substr(spacePos + 1);
			}
			spacePos = ended_at.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				end_date_str = ended_at.substr(0, spacePos);
				end_time_str = ended_at.substr(spacePos + 1);
			}

			int start_hour_int, start_minute_int, start_second_int; int16_t start_year_int; int8_t start_month_int, start_day_int;
			int end_hour_int, end_minute_int, end_second_int; int16_t end_year_int; int8_t end_month_int, end_day_int;
			if (sscanf(start_time_str.c_str(), "%d:%d:%d", &start_hour_int, &start_minute_int, &start_second_int) == 3) 
			{
				if (start_hour_int < start_hour || start_hour_int > end_hour-1) 
				{
					out_of_time++; continue;
				}
			}
			if (sscanf(end_time_str.c_str(), "%d:%d:%d", &end_hour_int, &end_minute_int, &end_second_int) == 3) 
			{
				if (end_hour_int < start_hour || end_hour_int > end_hour-1) 
				{
					out_of_time++; continue;
				}
			}
			if (sscanf(start_date_str.c_str(), "%hd-%hhd-%hhd", &start_year_int, &start_month_int, &start_day_int) == 3 &&
				sscanf(end_date_str.c_str(), "%hd-%hhd-%hhd", &end_year_int, &end_month_int, &end_day_int) == 3) {
				if (end_day_int > start_day_int) {
					out_of_time++; continue;
				}
			}
			
			double duration_seconds = 0.0;
			// Parse the date and time strings
			struct std::tm time1 = {};
			struct std::tm time2 = {};
			if (strptime((start_date_str + " " + start_time_str).c_str(), "%Y-%m-%d %H:%M:%S", &time1) &&
				strptime((end_date_str + " " + end_time_str).c_str(), "%Y-%m-%d %H:%M:%S", &time2)) {
				// Calculate the time difference in seconds
				duration_seconds = std::difftime(std::mktime(&time2), std::mktime(&time1));

				//std::cout << "Duration in seconds: " << duration_seconds << std::endl;
			} else {
				std::cerr << "Failed to parse time strings." << std::endl;
				no_time_difference++; std::cout << start_time_str << " " << end_time_str << std::endl; continue;
			}
			if(duration_seconds<60.0) 
			{
				no_time_difference++; continue;
			}

						
			auto it_start = station_name_map.find(start_station_name);
			auto it_end = station_name_map.find(end_station_name);
				
			int start_id = -1; int end_id = -1;
			auto it_start_shortName = bss_id_shortName_map.find( start_station_id_str );
			auto it_end_shortName = bss_id_shortName_map.find( end_station_id_str );

			if (it_start_shortName != bss_id_shortName_map.end() && it_end_shortName != bss_id_shortName_map.end()) {
				nb_found_by_name++;
				start_id = it_start_shortName->second;
				end_id = it_end_shortName->second;
			}
			else if (it_start != station_name_map.end() && it_end != station_name_map.end()) {
				nb_found_by_name++;
				start_id = it_start->second.id;
				end_id = it_end->second.id;
			}
			else {
				auto it_start_lat = station_lat_map.find(start_lat_double);
				auto it_start_lon = station_lon_map.find(start_lon_double);
				auto it_end_lat = station_lat_map.find(end_lat_double);
				auto it_end_lon = station_lon_map.find(end_lon_double);

				if (it_start_lat != station_lat_map.end() && it_start_lon != station_lon_map.end() &&
					it_end_lat != station_lat_map.end() && it_end_lon != station_lon_map.end()) {
					//nb_found_from_coord++;
				}
			}
			
			if(start_id == -1 && end_id == -1)
			{
				no_id_counter++;
				//trip.Show();
				//std::cout << "start_lat:" << start_lat << " start_lon:" << start_lon << std::endl; 
				//std::cout << "end_lat:" << end_lat << " end_lon:" << end_lon << std::endl;
				//getchar();
				continue; //printf("Could not find station ids ... \n"); //getchar();
			}
			
			TripDataAllInts all_ints_trip;
			all_ints_trip.start_year = start_year_int; all_ints_trip.start_month = start_month_int; all_ints_trip.start_day = start_day_int;
			all_ints_trip.start_hour = start_hour_int; all_ints_trip.start_minute = start_minute_int; all_ints_trip.start_second = start_second_int;
			
			all_ints_trip.end_year = end_year_int; all_ints_trip.end_month = end_month_int; all_ints_trip.end_day = end_day_int;
			all_ints_trip.end_hour = end_hour_int; all_ints_trip.end_minute = end_minute_int; all_ints_trip.end_second = end_second_int;
			
			all_ints_trip.start_id = start_id;
			all_ints_trip.end_id = end_id;
			all_ints_trip.duration = duration_seconds;
			data.all_int_trips.push_back(all_ints_trip);
			//all_ints_trip.Show(); getchar();
		}

	}catch (const io::error::no_digit& e) {
        std::cerr << "Error: Couldn't read the CSV file. " << e.what() << std::endl; 
    } catch (MyException& e) {
        std::cerr << "Error: An unknown error occurred while reading the CSV file." << e.what() << std::endl; 
    }

	// Inside your member function where you use the lambda
	auto minmax_dates = std::minmax_element(data.all_int_trips.begin(), data.all_int_trips.end(), [this](const TripDataAllInts& a, const TripDataAllInts& b) {
		return compare_tm(to_tm(a), to_tm(b));
	});


    if (minmax_dates.first != data.all_int_trips.end() && minmax_dates.second != data.all_int_trips.end()) {
        // Convert start date of earliest trip to std::tm
        std::tm earliest_date = to_tm(*minmax_dates.first);

        // Convert end date of latest trip to std::tm
        std::tm latest_date = to_tm(*minmax_dates.second);

        // Calculate the difference in days
        std::time_t earliest_time = std::mktime(&earliest_date);
        std::time_t latest_time = std::mktime(&latest_date);
        int difference_in_days = std::difftime(latest_time, earliest_time) / (60 * 60 * 24);

		// Print the earliest and latest dates
        std::cout << "Earliest date: " << (int)minmax_dates.first->start_year << "/" << (int)minmax_dates.first->start_month << "/" << (int)minmax_dates.first->start_day << std::endl;
        std::cout << "Latest date: " << (int)minmax_dates.second->start_year << "/" << (int)minmax_dates.second->start_month << "/" << (int)minmax_dates.second->start_day << std::endl;

        // Print the difference in number of days
        std::cout << "Difference between the start date of the earliest trip and the end date of the latest trip in terms of days:" << difference_in_days << " day counter:" << difference_in_days+1 << std::endl;
		data.nb_days_of_data = difference_in_days;
		data.dayCounter = difference_in_days+1;
    } else {
        std::cout << "No trips recorded. Exiting." << std::endl; exit(1);
    }
	
	if(line_counter != nb_found_by_name+nb_found_from_coord+exception_counter+no_id_counter+no_time_difference+trip_with_missing_data+out_of_time+next_day_counter)
	{
		int sum = nb_found_by_name+nb_found_from_coord+exception_counter+no_id_counter+no_time_difference+trip_with_missing_data+out_of_time;
		printf("Line_counter:%d sum:%d found_by_name:%d found_from_coord:%d exception:%d no_id:%d no_time_diff:%d missing_data:%d out_of_time:%d\nPausing program ...\n",line_counter,sum,nb_found_by_name,nb_found_from_coord,exception_counter,no_id_counter,no_time_difference,trip_with_missing_data,out_of_time);
	printf("The .csv has %d rows ...\n",line_counter);
	printf("Found %d by station name and %d by coordinates ... \n",nb_found_by_name,nb_found_from_coord);
	printf("There were %d rows with exceptions in long long ints or floats (trips that didn't have ending data)... \n",exception_counter);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d rows where the time difference could not be computed from the .csv ... \n",no_time_difference);
	printf("There were %d rows with some missing data ... \n",trip_with_missing_data);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... \n",out_of_time,start_hour,end_hour);
	printf("Total trip_with_missing_data for exel:%d\n",exception_counter+no_id_counter+no_time_difference+trip_with_missing_data);		
		getchar();
	}
	
	printf("The .csv has %d rows ...\n",line_counter);
	printf("Found %d by station name and %d by coordinates ... \n",nb_found_by_name,nb_found_from_coord);
	printf("There were %d rows with exceptions in long long ints or floats (trips that didn't have ending data)... \n",exception_counter);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d rows where the time difference could not be computed from the .csv ... \n",no_time_difference);
	printf("There were %d rows with some missing data ... \n",trip_with_missing_data);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... \n",out_of_time,start_hour,end_hour);
	printf("Total trip_with_missing_data for exel:%d\n",exception_counter+no_id_counter+no_time_difference+trip_with_missing_data);
	
	printf("Eliminating std::string and parsing all trips to pure integers of type int8_t ...\n");

}

void Load::Load_Bixi(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour)
{
	printf("In Load_Bixi reading %s and %s ...\n", filename_csv.c_str(), filename_json.c_str());
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
    clock_t lastPrintTime = startTime;

	data.dayCounter = 0; data.nb_days_of_data = 0;
	
	int nb_found_by_name = 0; int nb_found_from_coord = 0;
	
	// Open files
    std::ifstream file_csv(filename_csv); std::ifstream file_json(filename_json);
    if (!file_csv.is_open() || !file_json.is_open()) {
        std::cerr << "Error: Could not open files " << filename_csv << " " << filename_json << ". Exiting ..." <<std::endl;
        exit(1);
    }
	
	// Parse the JSON data;

	station_name_map.clear();
	station_lat_map.clear();
	station_lon_map.clear();
	
	nlohmann::json jsonData;
	file_json >> jsonData;
	nlohmann::json json_stations = jsonData["data"]["stations"];
	Station depot; depot.name = std::string("depot"); depot.id = 0; depot.lon = -73.6017523; depot.lat = 45.5298971; depot.cap = 0;
	data.stations.push_back(depot);
	//depot_row = {'station_id': 0, 'name': 'depot', 'lat': 45.5298971, 'lon': -73.6017523, 'capacity': 0}
	int cntr = 1;
	for (const auto& json_station : json_stations)
	{
		Station station;
		std::string station_name_str = json_station["name"]; station.name = station_name_str;
		std::string station_id_str = json_station["station_id"]; //station.id = atoi(station_id_str.c_str());
		station.id = cntr; cntr++; 
		station.lon = json_station["lon"];
		station.lat = json_station["lat"];
		station.cap = json_station["capacity"];
		//station.Show();
		//getchar();
		
		station_name_map.insert(std::pair<std::string,Station>(station.name,station));
		station_lat_map.insert(std::pair<double,Station>(station.lat,station));
		station_lon_map.insert(std::pair<double,Station>(station.lon,station));	

		//For the rebalancing instances
		station.bss_station_id = json_station["station_id"];
		
		data.stations.push_back(station);
	}
	file_json.close();
	printf("Stations:%d Sample station object:\n",(int)data.stations.size()); data.stations[0].Show();
	// --------------------------------------------------- //
	// Parse the CSV data
	
	int line_counter=0; int exception_counter=0; int no_id_counter =0; int out_of_time =0; int no_time_difference=0;
	std::string line;
	std::getline(file_csv, line); //the headers of the columns

    while (std::getline(file_csv, line))
	{		
		line_counter++;
		if (line_counter % 500000 == 0)
			std::cout << "Read " << line_counter << " lines from .csv ..." << std::endl;
		
		std::vector<std::string> tokens;
		size_t startPos = 0;
		// Find the positions of commas and extract substrings
		for (size_t pos = line.find(','); pos != std::string::npos; pos = line.find(',', startPos)) {
			tokens.push_back(line.substr(startPos, pos - startPos));
			startPos = pos + 1;
		}
		// Add the last token after the final comma
		tokens.push_back(line.substr(startPos));
		
		if(tokens.size() == 10)
		{
			long long start_time_ms; long long end_time_ms;
			double start_lat_double; double start_lon_double; double end_lat_double; double end_lon_double;
			try
			{
				start_time_ms = std::stoll(tokens[8]);
				end_time_ms = std::stoll(tokens[9]);
				start_lat_double = std::stod(tokens[2]); start_lon_double = std::stod(tokens[3]);
				end_lat_double = std::stod(tokens[6]); end_lon_double = std::stod(tokens[7]);				
			} catch (const std::exception& e)
			{
				exception_counter++;
				//std::cout << line << std::endl;
				//for(int i=0; i<tokens.size(); i++)
					//std::cout << tokens[i] << " ";
				//printf("\n");
				//getchar();
				//std::cerr << "Exception encountered: " << e.what() << std::endl;
				continue;
			}
			

            // Convert milliseconds to seconds
            time_t start_time_posix = start_time_ms / 1000;
            time_t end_time_posix = end_time_ms / 1000;
			
			long long duration_ms = end_time_ms - start_time_ms;
            long long duration_seconds = duration_ms / 1000;

			//std::cout << start_time_posix << " " << tokens[8] << std::endl;
			if(0<duration_seconds && duration_seconds<60) 
			{
				no_time_difference++; continue;
			}	

			//Local time returns time in time zone of my machine. Type in bash:date "+%Z"
            char start_buffer[80]; char end_buffer[80];
			strftime(start_buffer, sizeof(start_buffer), "%Y-%m-%d %H:%M:%S", localtime(&start_time_posix));
			strftime(end_buffer, sizeof(end_buffer), "%Y-%m-%d %H:%M:%S", localtime(&end_time_posix));
			//std::cout << "Posix:" << start_time_posix << " Formatted Start Date and Time: " << start_buffer << std::endl;
			//std::cout << "Posix:" << end_time_posix << " Formatted End Date and Time: " << end_buffer << std::endl;
			std::string start_buffer_str(start_buffer); std::string end_buffer_str(end_buffer);
			std::string start_date_str, start_time_str;
			std::string end_date_str, end_time_str;

			// Find the position of the space character separating date and time
			size_t spacePos = start_buffer_str.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				start_date_str = start_buffer_str.substr(0, spacePos);
				start_time_str = start_buffer_str.substr(spacePos + 1);
			}
			spacePos = end_buffer_str.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				end_date_str = end_buffer_str.substr(0, spacePos);
				end_time_str = end_buffer_str.substr(spacePos + 1);
			}
			
			int start_hour_int, start_minute_int, start_second_int;
			int end_hour_int, end_minute_int, end_second_int;
			if (sscanf(start_time_str.c_str(), "%d:%d:%d", &start_hour_int, &start_minute_int, &start_second_int) == 3) 
			{
				if (start_hour_int < start_hour || start_hour_int > end_hour-1) 
				{
					out_of_time++; continue;
				}
			}
			if (sscanf(end_time_str.c_str(), "%d:%d:%d", &end_hour_int, &end_minute_int, &end_second_int) == 3) 
			{
				if (end_hour_int < start_hour || end_hour_int > end_hour-1) 
				{
					out_of_time++; continue;
				}
			}
			
			std::string start_station_name = tokens[0];
			std::string end_station_name = tokens[4];
			
			auto it_start = station_name_map.find(start_station_name);
			auto it_end = station_name_map.find(end_station_name);

			int start_id = -1; int end_id = -1;
			if(it_start != station_name_map.end() && it_end != station_name_map.end())
			{
				nb_found_by_name++; 
				start_id =  it_start->second.id;
				end_id =  it_end->second.id;
			} else {
				auto it_start_lat = station_lat_map.find(start_lat_double);
				auto it_start_lon = station_lon_map.find(start_lon_double);
				auto it_end_lat = station_lat_map.find(end_lat_double);
				auto it_end_lon = station_lon_map.find(end_lon_double);
				
				
				if(it_start_lat != station_lat_map.end() && it_start_lon != station_lon_map.end() && it_end_lat != station_lat_map.end() && it_end_lon != station_lon_map.end() )
				{
					//nb_found_from_coord++; 
				}
			}

			if(start_id==-1 || end_id==-1)
			{
				no_id_counter++;
				//trip.Show();
				//std::cout << "start_lat:" << start_lat << " start_lon:" << start_lon << std::endl; 
				//std::cout << "end_lat:" << end_lat << " end_lon:" << end_lon << std::endl;
				//getchar();
				continue; //printf("Could not find station ids ... \n"); //getchar();
			}
			
			
			TripDataAllInts all_ints_trip;
			sscanf(start_date_str.c_str(), "%hd-%hhd-%hhd", &all_ints_trip.start_year, &all_ints_trip.start_month, &all_ints_trip.start_day);
			sscanf(end_date_str.c_str(), "%hd-%hhd-%hhd", &all_ints_trip.end_year, &all_ints_trip.end_month, &all_ints_trip.end_day);
			sscanf(start_time_str.c_str(), "%hhd:%hhd:%hhd", &all_ints_trip.start_hour, &all_ints_trip.start_minute, &all_ints_trip.start_second);
			sscanf(end_time_str.c_str(), "%hhd:%hhd:%hhd", &all_ints_trip.end_hour, &all_ints_trip.end_minute, &all_ints_trip.end_second);
			
			if(all_ints_trip.end_day > all_ints_trip.start_day)
			{
				out_of_time++; continue;
			}
			
			all_ints_trip.duration = duration_seconds;
			all_ints_trip.start_id = start_id;
			all_ints_trip.end_id = end_id;
			data.all_int_trips.push_back(all_ints_trip);		
			
		}else {
			std::cerr << "Failed to parse the .csv line." << std::endl;
		}
		
	}
	file_csv.close();

	if(line_counter != out_of_time+no_time_difference+no_id_counter+exception_counter+(int)data.all_int_trips.size())
	{
		int sum = out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size()+exception_counter;
		printf("line_counter:%d sum:%d out_of_time:%d no_time_difference:%d no_id_counter:%d exception_counter:%d loaded_trips:%d\n",line_counter,sum,out_of_time,no_time_difference,no_id_counter,exception_counter,(int)data.all_int_trips.size()); getchar();
	}

	printf("The .csv has %d rows ...\n",(int)data.all_int_trips.size()+exception_counter+no_id_counter+out_of_time);
	printf("Found %d by station name and %d by coordinates ... \n",nb_found_by_name,nb_found_from_coord);
	printf("There were %d rows with exceptions in long long ints or floats (trips that didn't have ending data)... \n",exception_counter);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... Continuing ... \n",out_of_time,start_hour,end_hour);
	printf("There were %d trips with no time difference.\n",no_time_difference);
	printf("Total trips with missing data for exel:%d\n",no_time_difference+no_id_counter+exception_counter);
	
	printf("Loaded %d trips and %d stations.\n",(int)data.all_int_trips.size(),(int)data.stations.size());
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;

	auto minmax_dates = std::minmax_element(data.all_int_trips.begin(), data.all_int_trips.end(), [this](const TripDataAllInts& a, const TripDataAllInts& b) {
		return compare_tm(to_tm(a), to_tm(b));
	});

    if (minmax_dates.first != data.all_int_trips.end() && minmax_dates.second != data.all_int_trips.end()) {
        // Convert start date of earliest trip to std::tm
        std::tm earliest_date = to_tm(*minmax_dates.first);

        // Convert end date of latest trip to std::tm
        std::tm latest_date = to_tm(*minmax_dates.second);

        // Calculate the difference in days
        std::time_t earliest_time = std::mktime(&earliest_date);
        std::time_t latest_time = std::mktime(&latest_date);
        int difference_in_days = std::difftime(latest_time, earliest_time) / (60 * 60 * 24);

        // Print the difference in number of days
        std::cout << "Difference between the start date of the earliest trip and the end date of the latest trip in terms of days: " << difference_in_days << " dayCounter:" <<  difference_in_days+1 << std::endl;
		data.nb_days_of_data = difference_in_days;
		data.dayCounter = difference_in_days+1;
    } else {
        std::cout << "No trips recorded." << std::endl;
    }
}	

void Load::Load_Bixi(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour, int max_stations_allowed)
{
	printf("In Load_Bixi with MaxNbStations reading %s and %s ...\n", filename_csv.c_str(), filename_json.c_str());
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
    clock_t lastPrintTime = startTime;
	
	data.dayCounter = 0; data.nb_days_of_data = 0;
	int nb_found_by_name = 0; int nb_found_from_coord = 0;
	
	//to count the data days ... Initialize variables to track earliest and latest dates
    std::tm earliestDate = {0};
    std::tm latestDate = {0};
    bool firstDate = true;
	data.dayCounter = 0; data.nb_days_of_data = 0;
	
    std::ifstream file_csv(filename_csv); std::ifstream file_json(filename_json);
    if (!file_csv.is_open() || !file_json.is_open()) {
        std::cerr << "Error: Could not open files " << filename_csv << " " << filename_json << ". Exiting ..." <<std::endl;
        exit(1);
    }
	
	// Parse the JSON data; 
	nlohmann::json jsonData;
	file_json >> jsonData;
	nlohmann::json json_stations = jsonData["data"]["stations"];
	Station depot; depot.name = std::string("depot"); depot.id = 0; depot.lon = -73.6017523; depot.lat = 45.5298971; depot.cap = 0;
	data.stations.push_back(depot);
	//depot_row = {'station_id': 0, 'name': 'depot', 'lat': 45.5298971, 'lon': -73.6017523, 'capacity': 0}
	int nb_added=1; int cntr=0;
	for (const auto& json_station : json_stations)
	{
		Station station;
		std::string station_name_str = json_station["name"]; station.name = station_name_str;
		std::string station_id_str = json_station["station_id"]; //station.id = atoi(station_id_str.c_str());
		station.id = cntr; cntr++; 
		station.lon = json_station["lon"];
		station.lat = json_station["lat"];
		station.cap = json_station["capacity"];
		//station.Show();
		//getchar();
		
		station_name_map.insert(std::pair<std::string,Station>(station.name,station));
		station_lat_map.insert(std::pair<double,Station>(station.lat,station));
		station_lon_map.insert(std::pair<double,Station>(station.lon,station));	

		//For the rebalancing instances
		station.bss_station_id = json_station["station_id"];
		
		data.stations.push_back(station);
		nb_added++;
		if(nb_added>max_stations_allowed)
			break;
	}
	file_json.close();
	printf("Loaded Stations:%d:\n",(int)data.stations.size()); 
	for(auto s : data.stations)
		s.Show();
	
	// --------------------------------------------------- //
	// Parse the CSV data
	
	int line_counter=0; int exception_counter=0; int no_id_counter =0; int out_of_time =0; int no_time_difference=0;
	std::string line;
	std::getline(file_csv, line); //the headers of the columns
	//std::cout << line << std::endl;

    while (std::getline(file_csv, line))
	{		
		line_counter++;
		if (line_counter % 500000 == 0)
			std::cout << "Read " << line_counter << " lines from .csv ..." << std::endl;
		
		std::vector<std::string> tokens;
		size_t startPos = 0;
		// Find the positions of commas and extract substrings
		for (size_t pos = line.find(','); pos != std::string::npos; pos = line.find(',', startPos)) {
			tokens.push_back(line.substr(startPos, pos - startPos));
			startPos = pos + 1;
		}
		// Add the last token after the final comma
		tokens.push_back(line.substr(startPos));
		
		if(tokens.size() == 10)
		{
			long long start_time_ms; long long end_time_ms;
			double start_lat_double; double start_lon_double; double end_lat_double; double end_lon_double;
			try
			{
				start_time_ms = std::stoll(tokens[8]);
				end_time_ms = std::stoll(tokens[9]);
				start_lat_double = std::stod(tokens[2]); start_lon_double = std::stod(tokens[3]);
				end_lat_double = std::stod(tokens[6]); end_lon_double = std::stod(tokens[7]);				
			} catch (const std::exception& e)
			{
				exception_counter++;
				//std::cout << line << std::endl;
				//for(int i=0; i<tokens.size(); i++)
					//std::cout << tokens[i] << " ";
				//printf("\n");
				//getchar();
				//std::cerr << "Exception encountered: " << e.what() << std::endl;
				continue;
			}
			

            // Convert milliseconds to seconds
            time_t start_time_posix = start_time_ms / 1000;
            time_t end_time_posix = end_time_ms / 1000;
			
			long long duration_ms = end_time_ms - start_time_ms;
            long long duration_seconds = duration_ms / 1000;

			//std::cout << start_time_posix << " " << tokens[8] << std::endl;
			if(0<duration_seconds && duration_seconds<60) 
			{
				no_time_difference++; continue;
			}	

			//Local time returns time in time zone of my machine. Type in bash:date "+%Z"
            char start_buffer[80]; char end_buffer[80];
			strftime(start_buffer, sizeof(start_buffer), "%Y-%m-%d %H:%M:%S", localtime(&start_time_posix));
			strftime(end_buffer, sizeof(end_buffer), "%Y-%m-%d %H:%M:%S", localtime(&end_time_posix));
			//std::cout << "Posix:" << start_time_posix << " Formatted Start Date and Time: " << start_buffer << std::endl;
			//std::cout << "Posix:" << end_time_posix << " Formatted End Date and Time: " << end_buffer << std::endl;
			std::string start_buffer_str(start_buffer); std::string end_buffer_str(end_buffer);
			std::string start_date_str, start_time_str;
			std::string end_date_str, end_time_str;

			// Find the position of the space character separating date and time
			size_t spacePos = start_buffer_str.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				start_date_str = start_buffer_str.substr(0, spacePos);
				start_time_str = start_buffer_str.substr(spacePos + 1);
			}
			spacePos = end_buffer_str.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				end_date_str = end_buffer_str.substr(0, spacePos);
				end_time_str = end_buffer_str.substr(spacePos + 1);
			}
			
			int start_hour_int, start_minute_int, start_second_int;
			int end_hour_int, end_minute_int, end_second_int;
			if (sscanf(start_time_str.c_str(), "%d:%d:%d", &start_hour_int, &start_minute_int, &start_second_int) == 3) 
			{
				if (start_hour_int < start_hour || start_hour_int > end_hour-1) 
				{
					out_of_time++; continue;
				}
			}
			if (sscanf(end_time_str.c_str(), "%d:%d:%d", &end_hour_int, &end_minute_int, &end_second_int) == 3) 
			{
				if (end_hour_int < start_hour || end_hour_int > end_hour-1) 
				{
					out_of_time++; continue;
				}
			}
			
			std::string start_station_name = tokens[0];
			std::string end_station_name = tokens[4];
			
			auto it_start = station_name_map.find(start_station_name);
			auto it_end = station_name_map.find(end_station_name);

			int start_id = -1; int end_id = -1;
			if(it_start != station_name_map.end() && it_end != station_name_map.end())
			{
				nb_found_by_name++; 
				start_id =  it_start->second.id;
				end_id =  it_end->second.id;
			} else {
				auto it_start_lat = station_lat_map.find(start_lat_double);
				auto it_start_lon = station_lon_map.find(start_lon_double);
				auto it_end_lat = station_lat_map.find(end_lat_double);
				auto it_end_lon = station_lon_map.find(end_lon_double);
				
				
				if(it_start_lat != station_lat_map.end() && it_start_lon != station_lon_map.end() && it_end_lat != station_lat_map.end() && it_end_lon != station_lon_map.end() )
				{
					//nb_found_from_coord++; 
				}
			}

			if(start_id==-1 || end_id==-1)
			{
				no_id_counter++;
				//trip.Show();
				//std::cout << "start_lat:" << start_lat << " start_lon:" << start_lon << std::endl; 
				//std::cout << "end_lat:" << end_lat << " end_lon:" << end_lon << std::endl;
				//getchar();
				continue; //printf("Could not find station ids ... \n"); //getchar();
			}
			
			
			TripDataAllInts all_ints_trip;
			sscanf(start_date_str.c_str(), "%hd-%hhd-%hhd", &all_ints_trip.start_year, &all_ints_trip.start_month, &all_ints_trip.start_day);
			sscanf(end_date_str.c_str(), "%hd-%hhd-%hhd", &all_ints_trip.end_year, &all_ints_trip.end_month, &all_ints_trip.end_day);
			sscanf(start_time_str.c_str(), "%hhd:%hhd:%hhd", &all_ints_trip.start_hour, &all_ints_trip.start_minute, &all_ints_trip.start_second);
			sscanf(end_time_str.c_str(), "%hhd:%hhd:%hhd", &all_ints_trip.end_hour, &all_ints_trip.end_minute, &all_ints_trip.end_second);
			
			if(all_ints_trip.end_day > all_ints_trip.start_day)
			{
				out_of_time++; continue;
			}
			
			all_ints_trip.duration = duration_seconds;
			all_ints_trip.start_id = start_id;
			all_ints_trip.end_id = end_id;
			data.all_int_trips.push_back(all_ints_trip);		
			
		}else {
			std::cerr << "Failed to parse the .csv line." << std::endl;
		}
		
	}
	file_csv.close();

	if(line_counter != out_of_time+no_time_difference+no_id_counter+exception_counter+(int)data.all_int_trips.size())
	{
		int sum = out_of_time+no_time_difference+no_id_counter+(int)data.all_int_trips.size()+exception_counter;
		printf("line_counter:%d sum:%d out_of_time:%d no_time_difference:%d no_id_counter:%d exception_counter:%d loaded_trips:%d\n",line_counter,sum,out_of_time,no_time_difference,no_id_counter,exception_counter,(int)data.all_int_trips.size()); getchar();
	}

	printf("The .csv has %d rows ...\n",(int)data.all_int_trips.size()+exception_counter+no_id_counter+out_of_time);
	printf("Found %d by station name and %d by coordinates ... \n",nb_found_by_name,nb_found_from_coord);
	printf("There were %d rows with exceptions in long long ints or floats (trips that didn't have ending data)... \n",exception_counter);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... Continuing ... \n",out_of_time,start_hour,end_hour);
	printf("There were %d trips with no time difference.\n",no_time_difference);
	printf("Total trips with missing data for exel:%d\n",no_time_difference+no_id_counter+exception_counter);
	
	printf("Loaded %d trips and %d stations.\n",(int)data.all_int_trips.size(),(int)data.stations.size());
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;

	auto minmax_dates = std::minmax_element(data.all_int_trips.begin(), data.all_int_trips.end(), [this](const TripDataAllInts& a, const TripDataAllInts& b) {
		return compare_tm(to_tm(a), to_tm(b));
	});

    if (minmax_dates.first != data.all_int_trips.end() && minmax_dates.second != data.all_int_trips.end()) {
        // Convert start date of earliest trip to std::tm
        std::tm earliest_date = to_tm(*minmax_dates.first);

        // Convert end date of latest trip to std::tm
        std::tm latest_date = to_tm(*minmax_dates.second);

        // Calculate the difference in days
        std::time_t earliest_time = std::mktime(&earliest_date);
        std::time_t latest_time = std::mktime(&latest_date);
        int difference_in_days = std::difftime(latest_time, earliest_time) / (60 * 60 * 24);

        // Print the difference in number of days
        std::cout << "Difference between the start date of the earliest trip and the end date of the latest trip in terms of days: " << difference_in_days << " dayCounter:" <<  difference_in_days+1 << std::endl;
		data.nb_days_of_data = difference_in_days;
		data.dayCounter = difference_in_days+1;
    } else {
        std::cout << "No trips recorded." << std::endl;
    }
}

void Load::Load_CapitalBikeShare(BSS_Data & data, std::string filename_csv, std::string filename_json, int start_hour, int end_hour, int max_stations_allowed)
{
	printf("In Load_CapitalBikeShare reading %s and %s ...\n", filename_csv.c_str(), filename_json.c_str());
	// Initialize elapsed time tracking variables
    clock_t startTime = clock();
    clock_t lastPrintTime = startTime;
	
	//to count the data days ... Initialize variables to track earliest and latest dates
    std::tm earliestDate = {0};
    std::tm latestDate = {0};
    bool firstDate = true;
	data.dayCounter = 0; data.nb_days_of_data = 0;	
	
	// Open files
    std::ifstream file_csv(filename_csv); std::ifstream file_json(filename_json);
    if (!file_csv.is_open() || !file_json.is_open()) {
        std::cerr << "Error: Could not open files " << filename_csv << " " << filename_json << ". Exiting ..." <<std::endl;
        exit(1);
    }
	
	// Parse the JSON data; 
	nlohmann::json jsonData;
	file_json >> jsonData;
	nlohmann::json json_stations = jsonData["data"]["stations"];
	
	Station depot; depot.name = std::string("depot"); depot.id = 0; depot.lon = -77.0742465; depot.lat = 38.895389; depot.cap = 0;
	data.stations.push_back(depot);
	//depot.Show();
	// Capital Bikeshare's headquarters are located at 1501 Wilson Blvd Ste 1100, Arlington, Virginia, 22209
	//38.89538947941347, -77.07424659999998
	int station_counter=1;
	for (const auto& json_station : json_stations)
	{
		//std::cout << json_station << std::endl;
		Station station;
		std::string station_name_str = json_station["name"]; station.name = station_name_str;
		//std::string station_id_str = json_station["station_id"]; station.id = atoi(station_id_str.c_str());
		station.id = station_counter;
		station.lon = json_station["lon"];
		station.lat = json_station["lat"];
		station.cap = json_station["capacity"];
		//station.Show();
		//getchar();
		data.stations.push_back(station);
		station_counter++;
		if(station_counter>max_stations_allowed)
			break;
	}
	file_json.close();
	printf("Sample station object:\n"); data.stations[0].Show();
	//exit(1);

	// --------------------------------------------------- //
	// Parse the CSV data
	int line_counter=0; int exception_counter=0; int no_id_counter =0; int out_of_time =0; int no_time_difference=0; int trip_with_missing_data=0;
	std::string line;
	std::getline(file_csv, line); //the headers of the columns
    
	//std::cout << line << std::endl;
	while (std::getline(file_csv, line))
	{		
		line_counter++;
		
		//std::cout << line << std::endl;
		//getchar();
		
		if (line_counter % 500000 == 0)
			std::cout << "Read " << line_counter << " lines from .csv ..." << std::endl;
		
			std::vector<std::string> tokens;
			std::string token;
			bool insideQuotes = false;

			for (char c : line) {
				if (c == '"') {
					insideQuotes = !insideQuotes; // Toggle the insideQuotes flag
				} else if (c == ',' && !insideQuotes) {
					// Found a comma outside of quotes, push the token
					if (!token.empty()) {
						if (token.front() == '"' && token.back() == '"') {
							token = token.substr(1, token.length() - 2);  // Remove opening and closing quotes
						}
						tokens.push_back(token);
					}
					token.clear(); // Clear the token for the next one
				} else {
					token += c; // Append the character to the current token
				}
			}

			// Push the last token (if any)
			if (!token.empty()) {
				if (token.front() == '"' && token.back() == '"') {
					token = token.substr(1, token.length() - 2);  // Remove opening and closing quotes
				}
				tokens.push_back(token);
			}

		
		if(tokens.size() == 13)
		{
			std::string start_buffer_str; std::string end_buffer_str;
			double start_lat; double start_lon; double end_lat; double end_lon;
			try
			{
				start_buffer_str = tokens[2];
				end_buffer_str = tokens[3];
				start_lat = std::stod(tokens[8]); start_lon = std::stod(tokens[9]);
				end_lat = std::stod(tokens[10]); end_lon = std::stod(tokens[11]);				
			} catch (const std::exception& e)
			{
				exception_counter++;
				//std::cout << line << std::endl;
				//for(int i=0; i<tokens.size(); i++)
					//std::cout << tokens[i] << " ";
				//printf("\n");
				//getchar();
				//std::cerr << "Exception encountered: " << e.what() << std::endl;
				continue;
			}
			
			std::string start_date_str, start_time_str;
			std::string end_date_str, end_time_str;

			// Find the position of the space character separating date and time
			size_t spacePos = start_buffer_str.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				start_date_str = start_buffer_str.substr(0, spacePos);
				start_time_str = start_buffer_str.substr(spacePos + 1);
			}
			spacePos = end_buffer_str.find(' ');
			if (spacePos != std::string::npos) {
				// Extract the date and time substrings
				end_date_str = end_buffer_str.substr(0, spacePos);
				end_time_str = end_buffer_str.substr(spacePos + 1);
			}
			
			int hour, minute, second;
			if (sscanf(start_time_str.c_str(), "%d:%d:%d", &hour, &minute, &second) == 3) 
			{
				if (hour < start_hour || hour > end_hour-1) 
				{
					out_of_time++; continue;
				}
			}
			
			double duration_seconds = 0.0;
			// Parse the date and time strings
			struct std::tm time1 = {};
			struct std::tm time2 = {};
			if (strptime((start_date_str + " " + start_time_str).c_str(), "%Y-%m-%d %H:%M:%S", &time1) &&
				strptime((end_date_str + " " + end_time_str).c_str(), "%Y-%m-%d %H:%M:%S", &time2)) {
				// Calculate the time difference in seconds
				duration_seconds = std::difftime(std::mktime(&time2), std::mktime(&time1));

				//std::cout << "Duration in seconds: " << duration_seconds << std::endl;
			} else {
				std::cerr << "Failed to parse time strings." << std::endl;
				no_time_difference++; std::cout << start_time_str << " " << end_time_str << std::endl; getchar();
			}
			if(duration_seconds<60.0) continue;
			//std::cout << start_time_str << std::endl;
			//getchar();
			
			TripData trip;	
			trip.start_date_str = start_date_str;
			trip.start_time_str = start_time_str;
			trip.end_date_str = end_date_str;
			trip.end_time_str = end_time_str;
			trip.duration = duration_seconds;
			trip.start_id = -1; trip.end_id = -1;
			
			for (const auto& station : data.stations)
			{
				if(std::abs(station.lat - start_lat) < 0.0001 && std::abs(station.lon - start_lon) < 0.0001)
				{
					trip.start_id = station.id; break;
				}
				if(station.name == tokens[4])
				{
					trip.start_id = station.id; break;
				}	
				
			}
			for (const auto& station : data.stations)
			{
				if(std::abs(station.lat - end_lat) < 0.0001 && std::abs(station.lon - end_lon) < 0.0001)
				{
					trip.end_id = station.id; break;
				}
				if(station.name == tokens[6])
				{
					trip.end_id = station.id; break;
				}		
			}
			if(trip.start_id==-1 || trip.end_id==-1)
			{
				no_id_counter++;
				//trip.Show();
				//std::cout << "start_lat:" << start_lat << " start_lon:" << start_lon << std::endl; 
				//std::cout << "end_lat:" << end_lat << " end_lon:" << end_lon << std::endl;
				//getchar();
				continue; //printf("Could not find station ids ... \n"); //getchar();
			}
			
			//trip.Show();
			//getchar();
			data.trips.push_back(trip);
			
			//to count the data days ...
			std::tm date = {0};
			std::istringstream dateStream(trip.start_date_str);
			dateStream >> std::get_time(&date, "%Y-%m-%d"); 
			if (firstDate) 
			{
				earliestDate = latestDate = date;
				firstDate = false;
			}else 
			{
				if (mktime(&date) < mktime(&earliestDate)) 
				{
					earliestDate = date;
				}
				if (mktime(&date) > mktime(&latestDate)) 
				{
					latestDate = date;
				}
			}
			if (data.dateToDayMap.find(trip.start_date_str) == data.dateToDayMap.end()) 
			data.dateToDayMap[trip.start_date_str] = data.dayCounter++;					
			
		}else {
			//std::cerr << ".csv line with missing data." << std::endl;
			//std::cout << "tokens size:" << (int)tokens.size() << std::endl;
			//for(const std::string& tok : tokens)
				//printf("%s \n",tok.c_str());
			//printf("\n");
			//getchar();
			trip_with_missing_data++;
		}
		
	}
	file_csv.close();
	printf("Sample trip object loaded:\n"); data.trips[0].Show();
	
	printf("The .csv has %d rows ...\n",(int)data.trips.size()+exception_counter+no_id_counter+out_of_time+no_time_difference+trip_with_missing_data);
	printf("There were %d rows with exceptions in long long ints or floats (trips that didn't have ending data)... \n",exception_counter);
	printf("There were %d rows that could not be added because station_id was not found from the .json file ... \n",no_id_counter);
	printf("There were %d rows where the time difference could not be computed from the .csv ... \n",no_time_difference);
	printf("There were %d rows with some missing data ... \n",trip_with_missing_data);
	printf("There were %d rows representing trips out_of_time [%d,%d] ... Continuing ... \n",out_of_time,start_hour,end_hour);
	
	printf("Loaded %d trips and %d stations.\n",(int)data.trips.size(),(int)data.stations.size());
	// Print elapsed time
    clock_t endTime = clock();
    double elapsedSeconds = (double)(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << std::setprecision(3) << elapsedSeconds << " seconds" << std::endl;	
	
	// Calculate the difference in days
    double diffInSeconds = difftime(mktime(&latestDate), mktime(&earliestDate));
    data.nb_days_of_data = (int)(diffInSeconds / (60 * 60 * 24));
	printf("nb_days_of_data:%d dayCounter:%d\n",data.nb_days_of_data,data.dayCounter);	
}	

