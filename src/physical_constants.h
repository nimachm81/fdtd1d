// Copyright 2018 Nima Chamanara.  All Rights Reserved.
// Use of this source code is governed by the GNU General Public License v3.0.


#ifndef FDTD_CONSTANTS_H_
#define FDTD_CONSTANTS_H_

#include <cmath>        // sqrt

#include "number_types.h"

namespace fdtd1d {

struct PhysicalConstants {
  // electric permittivity of vacuum (normalized units)
  static constexpr RealNumber epsilon_0 = static_cast<RealNumber>(1.0);
  
  // magnetic permeability of vacuum (normalized units)
  static constexpr RealNumber mu_0 = static_cast<RealNumber>(1.0);
  
  // the velocity of light in vacuum (normalized units)
  static constexpr RealNumber c = 
    static_cast<RealNumber>(1.0 / sqrt(epsilon_0*mu_0));
};

}  // namespace fdtd1d

#endif  // FDTD_CONSTANTS_H_

