// Use of this source code is governed by the GNU General Public License v3.0.


#ifndef FDTD_REALTYPE_H_
#define FDTD_REALTYPE_H_

// Defines floating point data type used in the FDTD algorithm
// and the integer types that describe indices for large arrays.
// Since the size of the arrays may pass INT_MAX, int64_t is provided as the
// int type.
// "RealNumber = double" provides higher accuracy, while "RealNumber = float"
// provides memory efficiency and it can be vectorized more efficiently.

#include <cstdint>    // std::int64_t


namespace fdtd1d {

using RealNumber = double;
using IntNumber = std::int64_t;

}  // namespace fdtd1d

#endif  // FDTD_REALTYPE_H_




