#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <vector>
#include <string>
#include <unordered_set>
#include "Edge.h"

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

class Optimizer;
class Renderer;

//the load operation should not be class method because of the readonly nature of the v_count and f_count vars.
class Polyhedron{
    public:
        Polyhedron(vector<vector<int>> f2v);
        //the operators from conway polyhedron notation
        Polyhedron dual();
        Polyhedron trunc_selected(vector<int> vertices);
        Polyhedron trunc(int degree);
        Polyhedron kis_selected(vector<int> faces);
        Polyhedron kis(int degree);
        Polyhedron ambo();
        Polyhedron join();
        Polyhedron gyro();
        Polyhedron propeller();
        Polyhedron triangulate(vector<int> faces,string type); //either "center" or "ear"
        //scrape off combines trunc, ambo, dual
        Polyhedron scrape_off(double t);
        Polyhedron detach_faces(double,double,double,bool);
        //face and vertex selection for fine grained control
        vector<int> select_vertices(int order);
        vector<int> select_faces(int order);
        //the save function -----------------------------------------
        void save_obj(string fname,bool overwrite,bool create_folders);
        //the checking function -------------------------------------
        void check();

        //string representation
        string string_repr();

        //positional adjustment functions
        void to_origin();
        void scale(double factor);
        void unskew(); //would like to make it canonical scale, but that has weird edge cases.
        void normalize(); //to_origin + scale + unskew
        void rotate(double angle,int axis);
        double align(int steps); //aligns the polyhedron by tumbling
        void spherical();
        const vec& get_center();

        //optimization functions
        double move_to_plane();
        double move_to_face_centers();
        double move_to_neighbors();
        double move_to_canonical();
        double move_to_equal_edge_length();

        Mx get_pos();
        void set_pos(Mx& pos);

        //optimizer creation. tries to speed up the optimization functions, but is not necessary.
        Optimizer optimizer();

        //weird functions --------------------------------------------
        Mx get_planar_drawing(int v0,int axis);

        Mx get_matrix(string type);

        void realize_spectral(Mx mat, vector<int> indices);
        void hang_axis(int);

        //state variables ---------------------------------------------
        void set_f2v_sorted(bool);
        void set_log(bool);

        //rendering functions -----------------------------------------
        MatrixXd render(int res, int thickness, int aa_factor);

        const Edges& get_edges(); //Edges are unordered_set<Edges..., see Edges.h.
        const vector<vector<int>>& get_v2v();


    //private:
        //the data ----------------------------------------------------
        vector<vector<int>> f2v_;
        vector<vector<int>> v2f_;
        vector<vector<int>> v2v_;
        vector<vector<Edge>> f2e_;
        vector<vector<Edge>> v2e_;
        Edges edges_;
        Mx pos_;
        Mx f2c_, f2n_, force_;
        vec center_;
        int v_count;
        int f_count;

        //state variables ---------------------------------------------
        bool f2c_computed;
        bool f2n_computed;
        bool force_initialized;
        bool center_computed;
        bool winding_fixed;
        bool v2f_sorted;
        bool f2v_sorted;
        // ------------------------------------------------------------

        //on-demand data generation & supply functions
        const vector<vector<int>>& get_f2v();
        const vector<vector<int>>& get_sorted_f2v();
        const vector<vector<int>>& get_v2f();
        const vector<vector<int>>& get_sorted_v2f();
        const vector<vector<Edge>>& get_f2e();
        const vector<vector<Edge>>& get_v2e();
        const Mx& get_f2n();
        const Mx& get_f2c();
        Mx& get_force();

        //sorting functions -------------------------------------------
        void sort_v2f();
        void sort_f2v();
        void fix_winding();
        void fix_insideout(); //a position operation
        Polyhedron flip_winding(); //to be able to generate both chiral options
        // ------------------------------------------------------------

        //resetting functions -----------------------------------------
        void position_changed();
        void reset();
        //-------------------------------------------------------------


        //the force class is a friend of the polyhedron class
        friend class Optimizer;
        friend class Renderer;

        //this is internally used during the spectral realization
        vector<int> select_indices(vector<int> indices,VectorXd D, Mx V,int mat_size);

        //wrapping vertex index to allow -1 etc.
        int fix_vertex_index(int);

        bool log_ = false;
        //logging functions
        void log(string);

};

inline int get_v_count(vector<vector<int>>& f2v){
    int v_count = 0;
    for(auto& f : f2v){
        for(auto& v : f){
            v_count = max(v_count,v);
        }
    }
    return v_count+1;
}

#endif
