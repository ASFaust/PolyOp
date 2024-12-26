#ifndef CIRCLEPACKINGREPR_H
#define CIRCLEPACKINGREPR_H

// a class to represent a polyhedron as a circle packing.
// this is used to find a 3d realization of a polyhedron

#include <vector>
#include <string>
#include <unordered_set>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "Polyhedron.h"

using namespace std;

namespace py = pybind11;

typedef Eigen::Vector3d vec;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> Mx;

vec fit_circle(const std::vector<vec>& on_other, const std::vector<vec>& on_this, vec center, double lr);

class CirclePackingRepr{
    public:
        CirclePackingRepr(Polyhedron& polyhedron_,vector<int> show_closest_ = vector<int>());
        Polyhedron& polyhedron;

        vector<vec> centers; //the centers of the circles in the packing. we do not need to store the radii because ||c||^2 + ||r||^2 = 1

        void step(double lr); //take a step in the optimization

        vector<vec> get_closest_points(int vertex);

        Mx get_pos(); //get the positions of the circles

        vector<vec> closest; //the closest points on the neighboring circles

        Mx draw(int res, double thickness); //draw the circle packing

        vector<int> show_closest; //the vertices to show the closest points for

        void move_to_origin(); //move the circle packing to the origin, eliminating a möbius transformation degree of freedom
};

#endif // CIRCLEPACKINGREPR_H