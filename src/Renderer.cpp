
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
    if (total_distance == 0) {
        double color = 1.0 - (z0 - z_min) / (z_max - z_min);
        draw_point(ix0, iy0, z0, color);
        return;
    }

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
    //also apply dof: draw distant lines blurred
    int max_dof = 5; //max kernel size. must be odd.

    double dof = 1.0 - (z - z_min) / (z_max - z_min);
    int kernel_size = (int)(dof * max_dof);
    MatrixXd& draw_point = get_draw_point(kernel_size); //predrawn blurred point of given thickness and blurred with given kernel size
    //then apply the draw_point to the image, but do a z test for each pixel
    for (int i = 0; i < draw_point.rows(); i++){
        int x_ = x + i - draw_point.rows() / 2;
        if (x_ < 0 && x_ >= width){
            continue;
        }
        for (int j = 0; j < draw_point.cols(); j++){
            int y_ = y + j - draw_point.cols() / 2;
            if (y_ >= 0 && y_ < height){
                if (z < zbuffer(x_, y_)){
                    image(x_, y_) = draw_point(i, j);
                    zbuffer(x_, y_) = z;
                }
            }
        }
    }
}

MatrixXd& Renderer::get_draw_point(int kernel_size){
    //we have a map<int, MatrixXd> draw_points;
    //where the key is the kernel_size

    //if the kernel_size is not in the map, create a new MatrixXd
    //and store it in the map
    //then return the MatrixXd
    double eps = 1e-6;
    if (draw_points.count(kernel_size) == 0) {
        int total_matrix_size = (kernel_size * 2) + (thickness * 2 + 1); //the total diameter of the matrix that we need.
        MatrixXd point = MatrixXd::Zero(total_matrix_size, total_matrix_size);
        for (int i = -thickness; i <= thickness; i++){
            for (int j = -thickness; j <= thickness; j++){
                if (i * i + j * j <= thickness * thickness){
                    int x = i + thickness + kernel_size;
                    int y = j + thickness + kernel_size;
                    point(x,y) = 1.0;
                }
            }
        }
        if(kernel_size > 0){
            double sigma = kernel_size / 3.0;
            //then blur the point with a gaussian kernel
            MatrixXd kernel = MatrixXd::Zero(kernel_size * 2 + 1, kernel_size * 2 + 1);
            for (int i = -kernel_size; i <= kernel_size; i++){
                for (int j = -kernel_size; j <= kernel_size; j++){
                    kernel(i + kernel_size, j + kernel_size) = exp(-(i * i + j * j) / (2 * sigma * sigma));
                }
            }
            kernel /= kernel.sum();
            //then convolve the point with the kernel
            //the point is already zero-padded
            for(int i = 0; i < total_matrix_size; i++){
                for(int j = 0; j < total_matrix_size; j++){
                    double sum = 0;
                    double count = 0;
                    for(int k = -kernel_size; k <= kernel_size; k++){
                        for(int l = -kernel_size; l <= kernel_size; l++){
                            int x = i + k;
                            int y = j + l;
                            if (x >= 0 && x < total_matrix_size && y >= 0 && y < total_matrix_size){
                                sum += point(x, y) * kernel(k + kernel_size, l + kernel_size);
                                count += kernel(k + kernel_size, l + kernel_size);
                            }
                        }
                    }
                    if (count > eps){
                        sum /= count;
                        if(sum > 1.0){
                            sum = 1.0;
                        }
                        if(sum > eps){
                            point(i, j) = sum;
                        }else{
                            point(i, j) = 0;
                        }
                    }else{
                        point(i, j) = 0;
                    }
                }
            }
        }
        cout << point << endl;
        draw_points[kernel_size] = point;
    }
    return draw_points[kernel_size];
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


