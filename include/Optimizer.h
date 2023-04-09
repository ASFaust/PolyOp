#ifndef OPTIMIZER_H
#define OPTIMIZER_H

class Polyhedron;

#include <Eigen/Dense>
using namespace Eigen;

typedef Vector3d vec;
typedef Matrix<double,Dynamic,Dynamic,RowMajor> Mx;

class Optimizer
{
    public:
        Optimizer(Polyhedron& polyhedron_);
        void step();
        void set_momentum(double decay, double strength);

    private:
        Polyhedron& polyhedron;
        Mx last_pos;
        Mx acc;
        double decay, strength;
};

#endif // EDGE_H