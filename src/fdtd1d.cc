// Copyright 2018 Nima Chamanara.  All Rights Reserved.
// Use of this source code is governed by the GNU General Public License v3.0.

#include <assert.h>

#include "fdtd1d.h"

namespace fdtd1d {

FDTD1D::FDTD1D() 
    : num_of_e_updated_threads_(0), num_of_h_updated_threads_(0),
      ind_t_(0), num_threads_(1), 
      output_file_name_("output.csv" /* default output file name */) {} 

void FDTD1D::SetXAxisRangeAndGridSpacing(const RealNumber x0, 
                                         const RealNumber x1, 
                                         const RealNumber dx) {
  x0_ = x0;
  x1_ = x1;
  num_x_ = static_cast<IntNumber>((x1 - x0) / dx);
  
  // dx_ is recalculated based on the total number of grid points num_x_ to take
  // into account the rounding error that was introduced in the floating point 
  // number to integer conversion in calculating num_x_
  dx_ = (x1 - x0) / num_x_;
}

void FDTD1D::InitializeAndResetEMFieldArrays() {
  // the H field points are staggered with respect to the E field points. Each 
  // H point is located between two E points. Therefore the number of H points
  // is smaller by 1 unit.
  e_field_.reset(new RealNumber[num_x_]);
  h_field_.reset(new RealNumber[num_x_ - 1]);
  
  for(IntNumber i=0; i<num_x_ - 1; ++i){
    e_field_[i] = 0.0;
    h_field_[i] = 0.0;
  }
  e_field_[num_x_ - 1] = 0.0;
}

void FDTD1D::SetStabilityFactorAndTimeResolution(
    const RealNumber stability_factor) {
  stability_factor_ = stability_factor;
  
  // the duration of time step is calculated using the grid spacing dx_ and the
  // specified stability factor
  dt_ = stability_factor*dx_ / PhysicalConstants::c;
}

void FDTD1D::SetSimulationTime(const RealNumber t_final) {
  t_final_ = t_final;
  num_t_ = static_cast<IntNumber>(t_final_ / dt_);
} 

void FDTD1D::SetNumberOfThreads(const int num_threads) {
  num_threads_ = num_threads;
  
  // calculate the chunk associated to each thread
  thread_data_chunk_bounds_.reset(new IntNumber[num_threads + 1]);
  thread_data_chunk_bounds_[0] = 0;
  thread_data_chunk_bounds_[num_threads] = num_x_;
  for(int i = 1; i < num_threads; ++i){
    thread_data_chunk_bounds_[i] = 
      static_cast<IntNumber>(i*std::lround(num_x_ / num_threads_));
  }
}

void FDTD1D::SetTheWriteToFileFlag(bool write_fields_to_file) {
  write_fields_to_file_ = write_fields_to_file;
}

void FDTD1D::InsertGaussianPointSource(const RealNumber position, 
                                       const RealNumber amplitude, 
                                       const RealNumber t_center, 
                                       const RealNumber t_decay) {
  GaussianSource j_gaussian(position, amplitude, t_center, t_decay);
  IntNumber ind_x = static_cast<IntNumber>((position - x0_)/dx_);
  j_gaussian.set_index_x(ind_x);
  point_sources_.emplace_back(j_gaussian);
}

void FDTD1D::UpdateElectricENodes(const int thread_index, 
                                  const UpdateState& nextState) {
  RealNumber dt_dx_eps0 = dt_/(dx_*PhysicalConstants::epsilon_0);
  
  // Maxwell-Ampere law 
  // update the inside of each chunk except its end points
  for (IntNumber i = thread_data_chunk_bounds_[thread_index] + 1;
      i < thread_data_chunk_bounds_[thread_index + 1]; ++i) {
      e_field_[i] -= (h_field_[i] - h_field_[i - 1])*dt_dx_eps0;
  }
  
  // update the end point of the chunk if necessary
  // each thread updates the beginning point of its associated chunk
  // the end point is therefore automatically updated by the next thread
  IntNumber ind_bound_0 = thread_data_chunk_bounds_[thread_index];
  if (thread_index != 0) {
    e_field_[ind_bound_0] -= 
        (h_field_[ind_bound_0] - h_field_[ind_bound_0 - 1])*dt_dx_eps0;
  }
  
  // takes into account the effect of the electric current
  RealNumber t = ind_t_*dt_;
  for (auto& j : point_sources_) {
    auto ind_j = j.get_index_x();
    auto j_value = j.GetCurrentValue(t);
    if (ind_j >= thread_data_chunk_bounds_[thread_index] && 
        ind_j < thread_data_chunk_bounds_[thread_index + 1]) {
      e_field_[ind_j] -= j_value*dt_dx_eps0;
    }
  }

  // job done --> increase num_of_e_updated_threads_
  num_of_e_updated_threads_ += 1;
  
  // if all the threads have completed their job set the state to the next state
  if (num_of_e_updated_threads_ >= num_threads_) {
    assert(num_of_e_updated_threads_ == num_threads_);
    //std::cout << num_of_e_updated_threads_ << std::endl;
    num_of_e_updated_threads_ = 0;
    update_state_ = nextState;
  }
}

void FDTD1D::UpdateMagneticHNodes(const int thread_index, 
                                  const UpdateState& nextState) {
  RealNumber dt_dx_mu0 = dt_/(dx_*PhysicalConstants::mu_0);

  // Maxwell-Faraday law 
  // the magnetic fields are all located inside the chunks
  // update the inside of each chunk except its last point inside the chunk
  for (IntNumber i = thread_data_chunk_bounds_[thread_index];
      i < thread_data_chunk_bounds_[thread_index + 1] - 1; ++i) {
    h_field_[i] -= (e_field_[i+1] - e_field_[i])*dt_dx_mu0;
  }
  
  // update the last point inside the chunk as well, except for the last thread,
  // for which it is calculated in the previous loop
  IntNumber ind_bound_1 = thread_data_chunk_bounds_[thread_index + 1] - 1;
  if (thread_index != num_threads_ - 1) {
    h_field_[ind_bound_1] -= (e_field_[ind_bound_1+1] - e_field_[ind_bound_1])*
                             dt_dx_mu0;
  }
  
  // job done --> increase num_of_h_updated_threads_
  num_of_h_updated_threads_ += 1;
  
  // if all the threads have completed their job set the state to the next state
  if (num_of_h_updated_threads_ >= num_threads_) {
    assert(num_of_h_updated_threads_ == num_threads_);
    //std::cout << num_of_h_updated_threads_ << std::endl;
    num_of_h_updated_threads_ = 0;
    update_state_ = nextState;
  }
}

// It starts with the initial state in kUpdateE, therefore the first while loop
// is executed and all the threads update the E nodes in their chunks, and the 
// next state is set to kUpdateH when the job is complete and all the threads
// have returned.
// then the second while loop is executed and all the threads update magnetic 
// fields in their corresponding chunk, and at the end the state is again set
// to kUpdateE.
// he process is repeated until the final time step is passed.
void FDTD1D::UpdateFieldsCuncurrently(const int thread_index) {
  UpdateState nextState = UpdateState::kUpdateH;
  for (IntNumber i = 0; i < num_t_; ++i) {
    if (thread_index == 0) {
      ind_t_ = i;
    }
    
    nextState = UpdateState::kUpdateH;
    while(update_state_ != UpdateState::kUpdateE) {
    }
    FDTD1D::UpdateElectricENodes(thread_index, nextState);

    nextState = UpdateState::kUpdateE;
    while(update_state_ != UpdateState::kUpdateH) {
    }
    FDTD1D::UpdateMagneticHNodes(thread_index, nextState);
  }
}

// similar to FDTD1D::UpdateFieldsCuncurrently except after the second while
// loop the fields are written to the output file.
void FDTD1D::UpdateFieldsAndWriteToFileCuncurrently(const int thread_index) {
  UpdateState nextState = UpdateState::kUpdateH;
  for (IntNumber i = 0; i < num_t_; ++i) {
    if (thread_index == 0) {
      ind_t_ = i;
    }

    nextState = UpdateState::kUpdateH;
    while(update_state_ != UpdateState::kUpdateE) {
    }
    FDTD1D::UpdateElectricENodes(thread_index, nextState);

    nextState = UpdateState::kUpdateOutputFile;
    while(update_state_ != UpdateState::kUpdateH) {
    }
    FDTD1D::UpdateMagneticHNodes(thread_index, nextState);
    
    nextState = UpdateState::kUpdateE;
    if (thread_index == 0) {
      WriteEfieldValuesToCSVFile(output_file_name_);
      update_state_ = nextState;
    }
    
  }
}

void FDTD1D::CreateThreadsAndRun() {
  std::vector<std::thread> threads;
  
  std::cout << "Initializing " << num_threads_ << " threads..." << std::endl;
  update_state_ = UpdateState::kUpdateE;
  for (int i = 0; i < num_threads_; ++i)
    if (!write_fields_to_file_) {
      threads.emplace_back(
        std::thread(&FDTD1D::UpdateFieldsCuncurrently, this, i));
    } else {
      remove(output_file_name_.c_str());  // write over the last file
      threads.emplace_back(
        std::thread(&FDTD1D::UpdateFieldsAndWriteToFileCuncurrently, this, i));
    }
  
  for (int i = 0; i < num_threads_; ++i) {
    threads[i].join();
  }
}

void FDTD1D::PrintEFieldValues() {
  std::cout << std::endl << "Electric field values: " << std::endl;
  for (int i = 0; i < num_x_; ++i) {
    std::cout << e_field_[i] << " ";
  }
  std::cout << std::endl;
}

void FDTD1D::SetOutputCSVFileName(const std::string& file_name) {
  output_file_name_ = file_name;
}

void FDTD1D::WriteEfieldValuesToCSVFile(const std::string& file_name) {
  std::ofstream ofs(file_name, std::ofstream::app);
  for (int i = 0; i < num_x_ - 1; ++i) {
    ofs << e_field_[i] << ", ";
  }
  ofs << e_field_[num_x_ - 1] << std::endl;
  ofs.close();
}

void FDTD1D::PrintParameters() {
  std::cout << "Grid: " << std::endl;
  std::cout << "x0 : " << x0_ << std::endl;
  std::cout << "x1 : " << x1_ << std::endl;
  std::cout << "dx : " << dx_ << std::endl;
  std::cout << "t1 : " << t_final_ << std::endl;
  std::cout << "dt : " << dt_ << std::endl;
  std::cout << "Nx : " << num_x_ << std::endl;
  std::cout << "Nt : " << num_t_ << std::endl;
  
  
  std::cout << "\nThread data chunk bounds: " << std::endl;
  for (int i = 0; i <= num_threads_; ++i) {
    std::cout << thread_data_chunk_bounds_[i] << " ";
  }
  std::cout << std::endl;
  
  
  std::cout << "\nPoint sources: " << std::endl;
  for (auto J : point_sources_) {
    J.PrintParameters();
  }
}

}  // namespace fdtd1d


