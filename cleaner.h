/*cleaner.h*/

/**
 * @brief Cleans noise from the incoming ECG signal.
 *
 * This is the header file that enables you to call clean_ecg(),
 * provided in compiled form in file "cleaner.o".
 */

#pragma once

#include <vector>

 /**
  * @brief Clean an ECG signal by removing noise and artifacts.
  *
  * This function processes the input ECG signal to remove noise and artifacts,
  * returning a cleaned version of the signal. The cleaning algorithm may include
  * filtering, baseline correction, and other signal processing techniques. If the 
  * function returns false, this means the signal could not be cleaned, and cannot
  * be trusted as a valid signal.
  *
  * @param ecg_data Vector of raw ECG signal samples.
  * @return true if cleaning was successful, false otherwise.
  */
bool clean_ecg(std::vector<double>& ecg_data);
