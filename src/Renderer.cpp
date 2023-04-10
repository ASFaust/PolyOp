
#include "Polyhedron.h"
#include "Renderer.h"
/*
class Renderer{
    public:
        Renderer(
            Polyhedron &poly_,
            int width_,
            int height_,
            double rx_,
            double ry_,
            int aa_level_,
            int num_threads_);
        //render to a eigen double matrix
        Eigen::MatrixXd& render();
    private:
        void draw_line(double x0, double y0, double z0, double x1, double y1, double z1);
        void downsampling(); //fills image_downsampled, level is the downsampling factor, aka aa_level

        Eigen::MatrixXd image, zbuffer, image_downsampled;

        Polyhedron &poly;
        int width;
        int height;
        double rx;
        double ry;
        int aa_level;
        int num_threads;
};
*/

Renderer::Renderer(
    Polyhedron &poly_,
    int width_,
    int height_,
    double rx_,
    double ry_,
    int aa_level_,
    int num_threads_,
    int thickness_):
    poly(poly_),
    width(width_ * aa_level_),
    height(height_ * aa_level_),
    rx(rx_),
    ry(ry_),
    true_width(width_),
    true_height(height_),
    aa_level(aa_level_),
    num_threads(num_threads_),
    thickness(thickness_),
    sq_thickness(thickness_ * thickness_)
{
    z_far = 1000;
    image = Eigen::MatrixXd::Zero(height, width);
    //initialize the zbuffer to be very far away
    zbuffer = Eigen::MatrixXd::Constant(height, width, z_far);
    image_downsampled = Eigen::MatrixXd::Zero(true_height, true_width);
}

Eigen::MatrixXd& Renderer::render(){
    //draw the polyhedron
    calculate_zrange();
    auto& edges = poly.get_edges();
    for(auto& e : edges){
        vec pos1 = poly.pos_.row(e.v1);
        vec pos2 = poly.pos_.row(e.v2);
        draw_line(pos1(0), pos1(1), pos1(2), pos2(0), pos2(1), pos2(2));
    }
    //downsample the image
    downsampling();
    return image_downsampled;
}

void Renderer::draw_line(double x0, double y0, double z0, double x1, double y1, double z1) {
    int ix0 = (x0 + rx) * width / (2 * rx);
    int iy0 = (y0 + ry) * height / (2 * ry);
    int ix1 = (x1 + rx) * width / (2 * rx);
    int iy1 = (y1 + ry) * height / (2 * ry);

    // Ensure the coordinates are within bounds
    if((ix0 < 0) || (ix0 >= width) || (iy0 < 0) || (iy0 >= height) ||
       (ix1 < 0) || (ix1 >= width) || (iy1 < 0) || (iy1 >= height)) {
        return;
    }

    int dx = abs(ix1 - ix0);
    int dy = abs(iy1 - iy0);
    int sx = (ix0 < ix1) ? 1 : -1;
    int sy = (iy0 < iy1) ? 1 : -1;
    int err = dx - dy;

    double total_distance = sqrt(dx * dx + dy * dy);

    while (true) {
        // Calculate the fractional progress along the line
        double progress = sqrt((ix1 - ix0) * (ix1 - ix0) + (iy1 - iy0) * (iy1 - iy0)) / total_distance;

        // Interpolate the z value and color based on the progress
        double z = z0 + (1.0 - progress) * (z1 - z0);
        //use z_max and z_min to normalize z:
        double color = 1.0 - (z - z_min) / (z_max - z_min);

        if (color < 0.1) {
            color = 0.1;
        }

        draw_point(ix0, iy0, z, color);

        if (ix0 == ix1 && iy0 == iy1) {
            break;
        }

        int e2 = 2 * err;

        if (e2 > -dy) {
            err -= dy;
            ix0 += sx;
        }

        if (e2 < dx) {
            err += dx;
            iy0 += sy;
        }
    }
}

void Renderer::calculate_zrange(){ //fills z_max and z_min by going through poly.pos_
    z_max = -1e9;
    z_min = 1e9;
    for (int i = 0; i < poly.pos_.rows(); i++){
        double z = poly.pos_(i, 2);
        if (z > z_max){
            z_max = z;
        }
        if (z < z_min){
            z_min = z;
        }
    }
}

void Renderer::draw_point(int x, int y, double z, double color){
    //dont forget to check the zbuffer at each pixel
    for (int i = -thickness; i <= thickness; i++){
        for (int j = -thickness; j <= thickness; j++){
            if (i * i + j * j <= sq_thickness){
                if (z < zbuffer(y + i, x + j)){
                    image(y + i, x + j) = color;
                    zbuffer(y + i, x + j) = z;
                }
            }
        }
    }
}

void Renderer::downsampling(){
    //downsample the image
    //the downsampling factor is aa_level
    //the image is stored in image_downsampled

    //for each pixel in image_downsampled
    for (int i = 0; i < true_height; i++){
        for (int j = 0; j < true_width; j++){
            //compute the average of the pixels in image
            //that correspond to the pixel in image_downsampled
            //and store the average in image_downsampled
            double sum = 0;
            for (int k = 0; k < aa_level; k++){
                for (int l = 0; l < aa_level; l++){
                    sum += image(i * aa_level + k, j * aa_level + l);
                }
            }
            image_downsampled(i, j) = sum / (aa_level * aa_level);
        }
    }
}
