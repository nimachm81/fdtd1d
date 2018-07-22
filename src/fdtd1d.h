// Copyright 2018 Nima Chamanara.  All Rights Reserved.
// Use of this source code is governed by the GNU General Public License v3.0.


#ifndef FDTD_FDTD1D_H_
#define FDTD_FDTD1D_H_

// Simulates one dimensional Maxwell equations on a uniform grid using the 
// finite difference time domain (FDTD) method. Maxwell equations:
//
// curl E(x, t) = -mu_0*dH(x, t)/dt                    Maxwell-Faraday
// curl H(x, t) = J(x, t) + epsilon_0*dE(x, t)/dt      Maxwell-Ampere
//
// are descretized and solved for the electric (E) and magnetic (H) fields 
// using finite differences.The E and H fields are defined on a Yee staggered
// grid. For more details see:
// "Taflove, A., & Hagness, S. C. (2005). Computational electrodynamics: 
// the finite-difference time-domain method. Artech house."
//

#include <string>         // std::string
#include <cmath>          // std::lround
#include <iostream>       // std::cout
#include <fstream>        // std::ofstream
#include <vector>         // std::vector
#include <thread>         // std::thread
#include <memory>         // std::unique_ptr
#include <atomic>         // std::atomic

#include "em_source.h"
#include "number_types.h"
#include "physical_constants.h"

namespace fdtd1d {

enum class UpdateState {
  kUpdateE = 1,
  kUpdateH = 2,
  kUpdateOutputFile = 3,
};


class FDTD1D {
  public:
  FDTD1D();
  void SetXAxisRangeAndGridSpacing(const RealNumber x0, const RealNumber x1, 
                              const RealNumber dx);
  void InitializeAndResetEMFieldArrays();
  void SetStabilityFactorAndTimeResolution(const RealNumber stability_factor);
  void SetSimulationTime(const RealNumber t_final);
  void SetNumberOfThreads(const int num_threads);
  
  int get_num_threads();
  void PrintParameters();
  
  // adds a gaussian point source to the problem
  void InsertGaussianPointSource(const RealNumber position, 
                                 const RealNumber amplitude, 
                                 const RealNumber t_center, 
                                 const RealNumber t_decay);
  // at each time step the electric fields are updated using this function based
  // on the Maxwell-Ampere equation
  void UpdateElectricENodes(const int thread_index, 
                            const UpdateState& nextState);
                            
  // at each time step the magnetic fields are updated using this function based
  // on the Maxwell-Faraday equation
  void UpdateMagneticHNodes(const int thread_index, 
                            const UpdateState& nextState);
                            
  void UpdateFieldsCuncurrently(const int thread_index);
  void UpdateFieldsAndWriteToFileCuncurrently(const int thread_index);
  void CreateThreadsAndRun();
  
  // prints the values of the electric field at the end of the simulation
  void PrintEFieldValues();
  
  void SetOutputCSVFileName(const std::string& file_name);
  void SetTheWriteToFileFlag(bool write_fields_to_file);
  void WriteEfieldValuesToCSVFile(const std::string& file_name);
  
  private:
  RealNumber x0_;               // [x0_, x1_] : computational domain range
  RealNumber x1_;
  RealNumber dx_;               // grid point spacing
  RealNumber t_final_;          // simulation stops at t_final_
  RealNumber dt_;               // duration of each time step
  IntNumber num_x_;             // total number of spatial grid points
  IntNumber num_t_;             // total number of time steps
  IntNumber ind_t_;             // current time index ---> t = ind_t_*dt
  RealNumber stability_factor_; // numerical stability factor
  
  // the electric (e) and magnetic (h) field arrays
  std::unique_ptr<RealNumber[]> e_field_ = nullptr;
  std::unique_ptr<RealNumber[]> h_field_ = nullptr;
  
  // the electric current sources are hold in this container 
  std::vector<GaussianSource> point_sources_;
  
  // sets the current state of the algorithm. The threads watch the state of 
  // this flag to decide their next action.
  // update_state_ = kUpdateE ---> the electric field should be updated
  // update_state_ = kUpdateH ---> the magnetic field should be updated
  // update_state_ = kUpdateOutputFile ---> the fields should be written to 
  //   the output file
  std::atomic<UpdateState> update_state_;
  
  int num_threads_ = 1;         // number of threads
  
  // holds the beginning and end of each chunk of grid handled by a given thread 
  std::unique_ptr<IntNumber[]> thread_data_chunk_bounds_ = nullptr;

  // counts the number of thread that completed their job at a given states
  std::atomic<int> num_of_e_updated_threads_;
  std::atomic<int> num_of_h_updated_threads_;
  
  // write the output electric field to the output file after each time step
  // the saved values can then be used to visualize the fields.
  bool write_fields_to_file_ = false;
  std::string output_file_name_;
  
};

}  // namespace fdtd1d

#endif  // FDTD_FDTD1D_H_


