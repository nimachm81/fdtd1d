// Use of this source code is governed by the GNU General Public License v3.0.

#include "em_source.h"

#include <iostream>       // std::cout


namespace fdtd1d {

GaussianSource::GaussianSource(const RealNumber position, 
                               const RealNumber amplitude, 
                               const RealNumber t_center, 
                               const RealNumber t_decay) {
  position_ = position;
  amplitude_ = amplitude;
  t_center_ = t_center;
  t_decay_ = t_decay;
}

void GaussianSource::set_index_x(const IntNumber ind_x) {
  index_x_ = ind_x;
}

IntNumber GaussianSource::get_index_x() {
  return index_x_;
}

RealNumber GaussianSource::GetCurrentValue(const RealNumber t) {
  auto temp = (t - t_center_) / t_decay_;
  return amplitude_*exp(-temp*temp);
}

void GaussianSource::PrintParameters() {
  std::cout << "position : " << position_ << std::endl;
  std::cout << "amplitude : " << amplitude_ << std::endl;
  std::cout << "t_center : " << t_center_ << std::endl;
  std::cout << "t_decay : " << t_decay_ << std::endl;
  std::cout << "ind_x : " << index_x_ << std::endl;
}

}  // namespace fdtd1d

