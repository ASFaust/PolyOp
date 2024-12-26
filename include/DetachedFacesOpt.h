#ifndef DETACHEDFACESOPT_H
#define DETACHEDFACESOPT_H

#include <vector>
#include <string>
#include <unordered_set>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "Polyhedron.h"

using namespace std;
using namespace Eigen;

typedef Vector3d vec;
typedef Matrix<double,Dynamic,Dynamic,RowMajor> Mx;

class DetachedFacesOpt{
    public:
        DetachedFacesOpt(Polyhedron& polyhedron_);
        Polyhedron& polyhedron;

        vector<vector<vec>> detached_faces; //the faces that are detached from the polyhedron

        void step(double lr); //take a step in the optimization
};

#endif