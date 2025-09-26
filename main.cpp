/**
* @brief Main program for AED software that decides to shock, or not.
* 
* Step 1. Output a start message & a done message
* Step 2. Open "ecg.dat" & Read and store the data (3 vectors)
*        & Output some values
* Step 3. Call clean_ecg (boolean func) & Bring "cleaner.h"
* Step 4. Call visualize_ecg => Output: png file
* Step 5. Compute the baseline
* Step 6. Compute the average amplitude
* Step 7. Compute the BPM
* Step 8. Compute the uniformity
* Step 9. Determine {Y, N}
* @note Written by Andrew Kim
* @note Northwestern University
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>

#include "cleaner.h"
#include "visualize.h"

using namespace std;
/** Step 5 Function (vector -> double) - compute the baseline (median)
* This function processes the input - which is the vector<double> of ECG values (y-values) - and sorts into
* an ascending order, so that we can find median.
*
* @param vector<double> datapoint
* @return double baseline
*/
double func_baseline(vector<double> datapoint) {
    double baseline=0.0; 
    sort(datapoint.begin(), datapoint.end()); // sort for finding median
    
    int n;
    
    // # of datapoint = {odd, even}
    if (datapoint.size()%2 == 1) {
    	n = datapoint.size()/2;
    	baseline = datapoint[n];
    }
    else {
    	n = datapoint.size()/2;
    	baseline = (datapoint[n-1]+datapoint[n])/2.0;
    }
    return baseline;
}

/** Step 6 Function (vector datapoint, vector rpeak, double baseline -> double) - compute the average amplitude
* This function reads the two vectors (datapoint, rpeak) and a double (baseline) to compute the average amplitude.
* The computing algorithm includes for loop through rpeak, so that we can extract rpeak = 1 cases.
* if statements are included because if there's no rpeak, than 0 is returned. 
*
* @param two vectors<double> (datapoint, rpeak) & double baseline
* @return double avg_amp (= computed average amplitude) => 0 if there are no rPeaks / avg_amp else
*/
// avg_amp = avg of abs(each Rpeak-baseline) - iff rpeak = 1
double func_avg_amp(vector<double> datapoint, vector<int> rpeak, double baseline) {
    double total=0.0, avg_amp;
    int rpeak1=0;
    
    // we have to extract rpeak = 1 cases
    for (size_t i=0; i<rpeak.size(); i++) {
        if (rpeak[i]==1) {
            total += abs(datapoint[i] - baseline);
            rpeak1++; // cnt denominator
        }
    }
    
    if (rpeak1!=0) {
	    avg_amp = total/rpeak1;
	    return avg_amp;
	  }
    else return 0;
}

/** Step 7 Function - compute the bpm
*                   (=60/average time between successive rpeak)
* This function processes vector<double> timestamp and vector<int> rpeak. By Looping through rpeak vector,
* we extract rpeak == 1 cases, and when there are two values, we extract the interval times.
* Then, we store timebetween, and cnt, and previous timestamp to continue extracting bpm.
*
* @param: vector<double> timestamp, vector<int> rpeak
* return: 0 if rpeak = 1 cases < 2 ; bpm(=60.0/avg ) else
*/

double func_bpm (vector<double> timestamp, vector<int> rpeak) {
		double timebetween=0.0, prev= -1.0, cnt = 0, avg;
		for (size_t i=0; i < rpeak.size(); i++) {
        if (rpeak[i]==1) { // Again, we only consider rpeak = 1
            if (prev!=-1.0) { // We need two values, and unless there are two values, prev > 0;
	            timebetween += timestamp[i] - prev;
	            cnt++; // needed for avg
            }
            prev = timestamp[i]; // updating previous timestamp
        }
    }
    
    if (cnt<2) return 0.0;
    else {
	    avg= timebetween/cnt;
	    double bpm = 60.0/avg;
	    return bpm;
	   }
}

/* Step 8 Function - compute the uniformity
* cnt = 0 ; Y1-Y2/R1-R2
	// increment cnt if Y1>baseline, Y2<=baseline, or vice versa
* uniformity = std of cnt
*
* This function inputs vector<double> datapoint, vector<int> rpeak, and double baseline.
* Data structure count is initialized to save the number of deflections in each rpeak=1 interval.
* This function incorporates nested for-loops: i loop to extract rpeak==1 cases, and
* j loop to loop through datapoints from the previous rpeak position until the curren rpeak position (rpeak pairs).
* j loop to loop through datapoints from the previous rpeak position until the curren rpeak position (rpeak pairs).
*
* After the loop, since rpeak < 2, uniformity = 0, check the # of rpeaks
* If rpeak >= 2, compute the standard deviation.
*
* @Param: vector<double> datapoint, vector<int> rpeak, double baseline
* @Return: 0 if rpeak < 2, else uniformity value
*/
double func_uniformity(vector<double> datapoint, vector<int> rpeak, double baseline) {
    vector<int> count;      // data structure to save count 
    int prev_peak= -1, rpeak1=0;
		
		// i loop: looping through rpeak & j loop: looping through datapoint
    for (size_t i=0; i<rpeak.size(); i++) {
        if (rpeak[i]==1) {
        	rpeak1++;
            if (prev_peak != -1) {
                int cnt=0; // cnt has to be reset when i changes
                // within previous rpeak and current rpeak 
                for (int j=prev_peak; j<(int)i; j++) {
                    double y1 = datapoint[j];
                    double y2 = datapoint[j+1];
										// conditionals to check deflections
                    if ((y1>baseline && y2<=baseline) || (y1<=baseline && y2>baseline)) {
                        cnt++;
                    }
                }
                count.push_back(cnt);
            }
            prev_peak = i;
        }
    }

    if (rpeak1<2) return 0.0;

    // Computing standard deviation (= uniformity)
    double sum=0.0, var=0.0, sd=0.0;
    for (int x : count) {
	    sum += x;
	  }
    double mean= sum/count.size();
    
    for (int x : count) {
	    var += (x-mean)*(x-mean);
    }
    var /= count.size();
    sd = sqrt(var);
    return sd;
}



// Function Section End

int main()
{
    cout << "** Starting AED Software **" << endl;
    cout << endl;
    
    // Step 2: opening "ecg.dat"
    ifstream file;
    file.open("ecg.dat");
    
    // check whether the file has been opened
    if (file.fail()) {
        cout << "ERROR: input file not found" << endl;
        return 0;
    }
    
    // if opened - create 3 vectors for timestamp, datapoint, rpeak
    vector<double> timestamp;
    vector<double> datapoint;
    vector<int> rpeak;
    
    // put data into vectors
    while (true) {
        double x, y;
        int z;
        
        // bring timestamp & datapoint (y value)
        file >> x >> y;
        if (file.fail()) break;
        timestamp.push_back(x);
        datapoint.push_back(y);
        
        // bring rpeak
        file >> z;
        if (file.fail()) break;
        rpeak.push_back(z);
    }
    
    // Step 3: Call clean_ecg function
    bool is_clean = clean_ecg(datapoint);
    
    // if statement as is_clean = {T, F}
    if (is_clean) {
        cout << "Is signal clean? YES" << endl;
    }
    else {
        cout << "Is signal clean? NO, DO NOT SHOCK" << endl;
        cout << endl;
        cout << "** Done **" << endl;
        return 0;
    }
    
    // Step 4: Visualize ECG Signal
    visualize_ecg(datapoint, timestamp); // calling the function
    
    // Step 5: Compute the baseline (func_baseline)
    double baseline = func_baseline(datapoint);
    cout << "Baseline? " << baseline << endl;
    
    // Step 6: Compute the average amplitude (func_avg_amp)
    double avg_amp = func_avg_amp(datapoint, rpeak, baseline);
    cout << "Average amplitude? " << avg_amp << endl;
    
    // Step 7: Compute the BPM (func_bpm)
    double bpm = func_bpm(timestamp, rpeak);
    cout << "BPM? " << bpm << endl;
    
    // Step 8: Compute the uniformity
    double uniformity = func_uniformity(datapoint, rpeak, baseline);
    string s;
   
    if (uniformity<1.0) {
	    s = "YES";
	}
    else s = "NO";
    cout << "Organized? " << s << " (" << uniformity << ")" << endl;
    
    // Step 9: Determine to shock or not
    if (avg_amp<0.1 || baseline>1.0 || (uniformity>=1.0 && bpm<200.0) || bpm<=150.0) {
	    s = "NO, DO NOT SHOCK";
	}
	else {
	    s = "YES, SHOCK!";
	}
	cout << "Shock patient? " << s << endl;
    
    
    
    /** Check data is printed correctly
    *cout << "Data:" << endl;
    *for (int i=0; i<5; i++) {
    *    cout << timestamp[i] << " " << datapoint[i] << " " << rpeak[i] << endl;
    *}
    *
    */
    cout << endl;
    cout << "** Done **" << endl;
    
    return 0;
}
