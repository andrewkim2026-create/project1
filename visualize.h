/*visualize.h*/

/**
 * @brief Plots the ECG signal.
 *
 * This is the header file that enables you to call visualize_ecg(),
 * provided in compiled form in file "visualize.o".
 */

#pragma once

#include <vector>

 /**
  * @brief Visualize the stored ECG signal and save a PNG plot using Gnuplot.
  *
  * This function:
  *  - Calls Gnuplot to render a PNG image of the ECG waveform.
  *
  * The output omage is written to the working directory with a unique filename:
  *   @c ecg_plot_YYYYMMDD_HHMMSS_N.png
  * where @c N is a random digit to avoid filename collisions.
  *
  * @param ecg_data Vector of raw ECG signal samples.
  * @param ecg_time_stamps Vector of time stamps, one per raw ECG data point.
  * @return void Nothing is returned, instead an output file is written.
  *
  * @note Requires Gnuplot to be installed and available on the system PATH.
  *
  * @warning Uses @c system() calls and @c rm for cleanup. Replace with safer
  *          process and file management routines if required by your environment.
  */
void visualize_ecg(std::vector<double> ecg_data, std::vector<double> ecg_time_stamps);

