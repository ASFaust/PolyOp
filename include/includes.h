#ifndef INCLUDES_H
#define INCLUDES_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>

namespace py = pybind11;

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <tuple>
#include <set>
#include <math.h>

using namespace std;
using namespace Eigen;

typedef Vector3d vec;
typedef Matrix<double,Dynamic,Dynamic,RowMajor> Mx;

#endif