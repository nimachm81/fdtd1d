// Use of this source code is governed by the GNU General Public License v3.0.


#ifndef FDTD_SOURCE_H_
#define FDTD_SOURCE_H_

// Defines electromagnetic sources with predefined temporal variations.

#include <cmath>    // exp

#include "number_types.h"

namespace fdtd1d {

// defines a point source that has a Gaussian temporal dependence. 
class GaussianSource {
  public:
  GaussianSource(const RealNumber position, 
                 const RealNumber amplitude, 
                 const RealNumber t_center, // see t_center_ for description
                 const RealNumber t_decay   // see t_decay_
                 );
                 
  void set_index_x(const IntNumber ind_x);
  IntNumber get_index_x();
  
  // Returns the value of the electromagnetic current at a given time
  RealNumber GetCurrentValue(const RealNumber t);
  void PrintParameters();
  
  private:
  // position of the point source
  RealNumber position_;
  
  // the amplitude of the Gaussian
  RealNumber amplitude_;
  
  // The Gaussian profile is centered in time around t_center_
  RealNumber t_center_;
  
  // The decay time of the Gaussian profile. the smaller t_decay_ the narrower
  // the generated electromagnetic pulse
  RealNumber t_decay_;
  
  // the grid index describing the position of the source on the x axis
  IntNumber index_x_;
};

}  //namespace fdtd1d

#endif  // FDTD_SOURCE_H_


