#ifndef RENDERER_H
#define RENDERER_H
#include "Polyhedron.h"

class Renderer{
    public:
        Renderer(
            Polyhedron &poly_,
            int width_,
            int height_,
            double rx_,
            double ry_,
            int aa_level_,
            int num_threads_,
            int thickness_);
        //render to a eigen double matrix
        Eigen::MatrixXd& render();
    private:
        void draw_line(double x0, double y0, double z0, double x1, double y1, double z1);
        void draw_point(int x, int y, double z, double color); //draws a point to the zbuffer with thickness :)
        void downsampling(); //fills image_downsampled, level is the downsampling factor, aka aa_level
        void calculate_zrange(); //fills z_max and z_min
        Eigen::MatrixXd image, zbuffer, image_downsampled;

        Polyhedron &poly;
        int width;
        int height;
        int true_width;
        int true_height;
        int thickness; // not scaled by aa_factor lol
        int sq_thickness;
        double rx;
        double ry;
        double z_far;
        double z_max, z_min;
        int aa_level;
        int num_threads; // currently not used
};

#endif