#include "CirclePackingRepr.h"
#include <iostream>

using namespace std;

CirclePackingRepr::CirclePackingRepr(Polyhedron& polyhedron_,vector<int> show_closest_): polyhedron(polyhedron_) {
    //initialize the centers

    const Mx pos = polyhedron.get_pos();
    centers = vector<vec>(pos.rows());
    for(int i = 0; i < pos.rows(); i++){
        centers[i] = pos.row(i).normalized();
    }
    vec center = vec::Zero(3);
    for(auto c : centers){
        center += c;
    }
    center /= centers.size();
    for(auto& c : centers){
        c -= center;
        c = c.normalized();
        c *= 0.999; //the initial circles are very small
    }
    show_closest = show_closest_;
}

void CirclePackingRepr::move_to_origin(){
    //move the circle packing to the origin
    vec center = vec::Zero(3);
    for(auto c : centers){
        center += c;
    }
    center /= centers.size();
    for(auto& c : centers){
        c -= center;
    }
}

struct TupleHash {
    template <typename T1, typename T2>
    size_t operator()(const tuple<T1, T2>& t) const {
        size_t h1 = hash<T1>()(get<0>(t)); // Hash first element
        size_t h2 = hash<T2>()(get<1>(t)); // Hash second element
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2)); // Combine hashes
    }
};

void CirclePackingRepr::step(double lr){
    //algorithm outline:
    //for each circle, we compute the points on the neighboring circles that are closest to the center of this circle.
    //these are coplanar since they lie on the unit sphere.
    //we then use the previously engineered algorithm to find a closest-matching circle to these points.
    //this is trivially possible thanks to the fact that the points are coplanar. if they were not, we'd have a problem.
    //this yields us the new center position. since the circles are very small in the beginning, we add a damping factor to the update.

    closest.clear();
    move_to_origin();

    unordered_map<tuple<int,int>,vec,TupleHash> points;

    auto v2v = polyhedron.get_v2v();

    //first compute all the closest points on the circles
    for(int i = 0; i < centers.size(); i++) {
        vector<vec> current_closest = get_closest_points(i);
        for(int j = 0; j < current_closest.size(); j++){
            points[make_tuple(i,v2v[i][j])] = current_closest[j];
        }
    }

    vector<vec> new_centers(centers.size());

    for(int i = 0; i < centers.size(); i++){
        vector<vec> on_other, on_this;
        for(int j = 0; j < v2v[i].size(); j++){
            vec p1 = points[make_tuple(i,v2v[i][j])];
            vec p2 = points[make_tuple(v2v[i][j],i)];
            on_other.push_back(p1);
            on_this.push_back(p2);
        }
        new_centers[i] = fit_circle(on_other,on_this,centers[i],lr);
    }

    //update the centers
    for(int i = 0; i < centers.size(); i++){
        centers[i] = new_centers[i];
    }

    //no normalizing!
}

vector<vec> CirclePackingRepr::get_closest_points(int vertex){
    auto& v2v = polyhedron.get_v2v();
    vector<vec> closest_points;
    vec c1 = centers[vertex];
    bool show = false;
    for(auto& s : show_closest){
        if(s == vertex){
            show = true;
            break;
        }
    }
    if(show){
        closest.push_back(c1);
    }
    // Loop over neighbors
    for(int j = 0; j < v2v[vertex].size(); j++) {
        // c2 is the neighbor center
        vec c2 = centers[v2v[vertex][j]];

        // We want p1 on the unit sphere with angle = angle(c2, p1)
        // and lying in the plane spanned by {0, c1, c2},
        // and specifically the one that is "closest" to c1.

        // 1) e2 = c2 / ||c2||
        double normC2 = c2.norm();
        if(normC2 < 1e-12) {
            // c2 is nearly zero => skip or handle edge case
            closest_points.push_back(c1);
            cerr << "Warning: c2 is nearly zero" << endl;
            continue;
        }
        vec e2 = c2 / normC2;

        // 2) w = c1 - (c1 ⋅ e2) e2, then e2_perp = w / ||w||
        vec w = c1 - (c1.dot(e2))*e2;
        double normW = w.norm();
        if(normW < 1e-12) {
            // c1 is parallel or nearly parallel to c2 => skip or handle edge case
            closest_points.push_back(c1);
            cerr << "Warning: c1 is parallel or nearly parallel to c2" << endl;
            continue;
        }
        vec e2_perp = w / normW;

        // 3) The circle on the unit sphere for c2 has angle theta2 = arccos(normC2).
        //    So points on that circle can be written:
        //        pCandidate = normC2*e2 ± sqrt(1 - normC2^2)*e2_perp.
        double rad = std::sqrt(1.0 - normC2*normC2);

        vec pCandidate1 = normC2*e2 + rad*e2_perp;  // "plus" direction
        vec pCandidate2 = normC2*e2 - rad*e2_perp;  // "minus" direction

        // Decide which candidate is closer to c1
        double dist1 = (pCandidate1 - c1).norm();
        double dist2 = (pCandidate2 - c1).norm();

        vec p1 = (dist1 < dist2) ? pCandidate1 : pCandidate2;

        p1 = p1.normalized(); // Normalize to ensure it's on the unit sphere

        // ------------------------------------------------------------
        // p1 is now the point on the circle from c2 that is closest to c1.

        closest_points.push_back(p1);


        if(show){
            closest.push_back(p1);
        }
    }
    return closest_points;
}

Mx CirclePackingRepr::get_pos(){
    //return the centers as a matrix
    Mx pos(centers.size(),3);
    for(int i = 0; i < centers.size(); i++){
        pos.row(i) = centers[i];
    }
    return pos;
}

vec fit_circle(const std::vector<vec>& on_other, const std::vector<vec>& on_this, vec center, double lr)
{
    // Algorithm to fit a circle to a set of points on the unit sphere.
    // The points are given in two sets: on_other and on_this. the points (0,0,0), center, on_other[i] and on_this[i] are coplanar.
    // we need to locally linearize the problem, since on a sphere, we cannot add forces in the same way as on a plane.
    const double radius = sqrt(1.0 - center.squaredNorm());

    vec force = vec::Zero(3);
    double r_force = 0;
    for(int i = 0; i < on_other.size(); i++){
        force += on_other[i] - on_this[i];
        vec d_other = on_other[i] - center;
        vec d_this = on_this[i] - center;
        r_force += d_other.norm() - d_this.norm();
    }
    force /= on_other.size();
    r_force /= on_other.size();
    vec new_center = center + force * lr;
    double new_r = radius + r_force * lr;
    //it holds: ||c||^2 + ||r||^2 = 1
    //therefore ||c|| = sqrt(1 - ||r||^2)
    new_center.normalize();
    new_center *= sqrt(1.0 - new_r*new_r);

    return new_center;
}

Mx CirclePackingRepr::draw(int res, double thickness) {
    // thickness is our "epsilon" for how close to the plane we must be
    Mx img(res, res);
    img.setZero();

    for (int y = 0; y < res; y++) {
        for (int x = 0; x < res; x++) {
            // Map pixel to [-1, 1] in both x and y
            // Center the pixel in the cell with +0.5
            double nx = 2.0 * (x + 0.5) / res - 1.0;  // in [-1,1]
            double ny = 2.0 * (y + 0.5) / res - 1.0;  // in [-1,1]

            double r2 = nx*nx + ny*ny;
            // If r2 > 1, then we are "off" the sphere in this orthographic projection
            if (r2 <= 1.0) {
                // z is the "upper" hemisphere. If you want the lower, use -sqrt(1 - r2).
                double nz = std::sqrt(1.0 - r2);

                // 3D point on the sphere
                vec p(nx, ny, nz);

                // Check each circle in your packing
                // If the pixel is within "thickness" of that circle, color it
                for (size_t i = 0; i < centers.size(); i++) {
                    const vec &c = centers[i];
                    double planeDist = p.dot(c) - c.squaredNorm();

                    if (std::fabs(planeDist) < thickness) {
                        // Mark this pixel as on the circle
                        // e.g., set to 1.0 (or some color value)
                        img(y, x) = 1.0;
                        break;
                    }
                }
            }
        }
    }
    //also draw the closest points. but only the ones that are visible. we do that by only drawing the 50% closest points.
    for(auto& p : closest){
        double nx = p[0];
        double ny = p[1];
        double x = (nx + 1.0) * res / 2.0;
        double y = (ny + 1.0) * res / 2.0;
        img(y, x) = 0.5;
    }

    return img;
}
