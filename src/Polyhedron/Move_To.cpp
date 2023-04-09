#include "Polyhedron.h"

double Polyhedron::move_to_plane(){
    auto& f2c = get_f2c();
    auto& f2n = get_f2n();
    auto& v2f = get_v2f();
    auto& f2v = get_f2v();
    auto& force = get_force(); // just an empty matrix

    for(int f = 0; f < f_count; f++){
        auto& center = f2c.row(f);
        auto& normal = f2n.row(f);
        for(int v : f2v[f]){
            force.row(v) -= normal * normal.dot(pos_.row(v) - center) / double(v2f[v].size());
        }
    }
    pos_ += force;
    position_changed();
    return force.rowwise().norm().mean() / v_count;
}

double Polyhedron::move_to_face_centers(){
    auto& f2c = get_f2c();
    auto& v2f = get_v2f();

    vec target = vec::Zero();
    double ret = 0;

    for(int v = 0; v < v_count; v++){
        target = vec::Zero();
        for(int f : v2f[v]){
            target += f2c.row(f);
            //force.row(v) += (f2c.row(f) - pos.row(v)) / double(v2f[v].size()) * strength;
        }
        target /= double(v2f[v].size());
        ret += (target.transpose() - pos_.row(v)).norm();
        pos_.row(v) = target;
    }
    position_changed();
    return ret / v_count;
}

double Polyhedron::move_to_neighbors(){
    auto& v2v = get_v2v();
    auto& v2f = get_v2f();
    //auto& force = get_force(); // just an empty matrix
    Vector3d target = Vector3d::Zero();
    double ret = 0;
    for(int v = 0; v < v_count; v++){
        target = Vector3d::Zero();
        for(int n : v2v[v]){
            target += pos_.row(n);
        }
        target /= double(v2f[v].size());
        ret += (target.transpose() - pos_.row(v)).norm();
        pos_.row(v) = target;
    }
    position_changed();
    return ret / v_count;
}

double Polyhedron::move_to_canonical(){  //they operate in-place
    auto& edges = get_edges();
    auto& v2f = get_v2f();
    vec center = get_center();
    auto& force = get_force(); // just an empty matrix

    for(auto& edge : edges){
        vec pos_v1 = pos_.row(edge.v1);
        vec pos_v2 = pos_.row(edge.v2);

        vec edge_center = (pos_v1 + pos_v2) / 2.0;
        vec dir_to_center = edge_center - center;
        vec edge_dir = (pos_v2 - pos_v1).normalized();
        double t = dir_to_center.dot(edge_dir);

        vec closest_point;

        if(t <= 0) {
            closest_point = pos_v1;
        } else if(t >= (pos_v2 - pos_v1).norm()) {
            closest_point = pos_v2;
        } else {
            closest_point = pos_v1 + t * edge_dir;
        }

        vec ideal_closest_point = closest_point.normalized();
        vec vertex_force = ideal_closest_point - closest_point;

        force.row(edge.v1) += vertex_force / v2f[edge.v1].size();
        force.row(edge.v2) += vertex_force / v2f[edge.v2].size();
    }
    pos_ += force;
    position_changed();
    return force.rowwise().norm().mean() / v_count;
}

double Polyhedron::move_to_equal_edge_length(){
    auto& edges = get_edges();
    auto& v2f = get_v2f();
    auto& force = get_force(); // just an empty matrix

    double avg_edge_length = 0;

    for(auto& edge : edges){
        avg_edge_length += (pos_.row(edge.v1) - pos_.row(edge.v2)).norm();
    }

    avg_edge_length /= edges.size();

    for(auto& edge : edges){
        vec pos_v1 = pos_.row(edge.v1);
        vec pos_v2 = pos_.row(edge.v2);
        vec delta = pos_v2 - pos_v1;
        vec mid = (pos_v1 + pos_v2) / 2.0;
        vec ideal_v1 = mid - delta.normalized() * avg_edge_length / 2.0;
        vec ideal_v2 = mid + delta.normalized() * avg_edge_length / 2.0;
        force.row(edge.v1) += (ideal_v1 - pos_v1) / v2f[edge.v1].size();
        force.row(edge.v2) += (ideal_v2 - pos_v2) / v2f[edge.v2].size();
    }

    pos_ += force;
    position_changed();
    return force.rowwise().norm().mean() / v_count;
}

Mx& Polyhedron::get_force(){
    if(!force_initialized){
        force_ = MatrixXd::Zero(v_count, 3);
        force_initialized = true;
    }
    force_.setZero();
    return force_;
}
