// Copyright 2018 Nima Chamanara.  All Rights Reserved.
// Use of this source code is governed by the GNU General Public License v3.0.

#include <chrono>       // chrono::steady_clock, chrono::duration

#include "number_types.h"
#include "fdtd1d.h"


int main() {
  auto t_start = std::chrono::steady_clock::now();

  // x axis grid parameters
  fdtd1d::RealNumber x0(-10.0);
  fdtd1d::RealNumber x1(10.0);
  fdtd1d::RealNumber dx(0.01);
  // simulation time and stability factor
  fdtd1d::RealNumber t_final(22.0);
  fdtd1d::RealNumber stabilityFactor(0.99); 
  // number of threads
  int num_threads(1);
  
  fdtd1d::FDTD1D fdtd;
  fdtd.SetXAxisRangeAndGridSpacing(x0, x1, dx);
  fdtd.InitializeAndResetEMFieldArrays();
  fdtd.SetStabilityFactorAndTimeResolution(stabilityFactor);
  fdtd.SetSimulationTime(t_final);
  fdtd.SetNumberOfThreads(num_threads);
  
  //electric source j
  fdtd1d::RealNumber j_position(0.0);
  fdtd1d::RealNumber j_amplitude(1.0);
  fdtd1d::RealNumber j_t_center(1.0);
  fdtd1d::RealNumber j_t_decay(0.2);
  fdtd.InsertGaussianPointSource(j_position, j_amplitude, 
                                 j_t_center, j_t_decay);
  
  fdtd.PrintParameters();
  fdtd.SetTheWriteToFileFlag(false);
  //fdtd.SetOutputCSVFileName("efiled-utput.csv");

  fdtd.CreateThreadsAndRun();
  
  //fdtd.PrintEFieldValues();
  
  auto t_end = std::chrono::steady_clock::now();
  std::chrono::duration<double> time_span(t_end - t_start);
  
  std::cout << std::endl << "It took " << time_span.count() 
            << " seconds." << std::endl;

  return 0;  
}

