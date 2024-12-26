#include "DetachedFacesOpt.h"

DetachedFacesOpt::DetachedFacesOpt(Polyhedron& polyhedron_) : polyhedron(polyhedron_){
    //initialize detached_faces
    detached_faces = vector<vector<vec>>();

    //initialize normals

    auto pos = polyhedron.get_pos(); //matrix of shape (v_count,3)
    //normalize the positions to be on a sphere
    for(int i = 0; i < polyhedron.v_count; i++){
        pos.row(i) /= pos.row(i).norm();
    }
    polyhedron.set_pos(pos);

    const auto& f2v = polyhedron.get_f2v();
    const auto& v2f = polyhedron.get_v2f();
    const auto& f2n = polyhedron.get_f2n();
    const auto& f2c = polyhedron.get_f2c();

    for(int f = 0; f < polyhedron.f_count; f++){
        vector<vec> detached_face;
        vec normal = f2n.row(f);
        vec center = f2c.row(f);
        for(int v : f2v[f]){
            vec cpos = pos.row(v);
            vec projected_pos = cpos - (cpos - center).dot(normal) * normal;
            detached_face.push_back(projected_pos);
        }
        detached_faces.push_back(detached_face);
    }
}

void DetachedFacesOpt::step(double lr){
    Mx mean_vertex_positions = Mx::Zero(polyhedron.v_count,3);
    const auto& f2v = polyhedron.get_f2v();
    const auto& v2f = polyhedron.get_v2f();

    for(int f = 0; f < polyhedron.f_count; f++){
        for(int i = 0; i < detached_faces[f].size(); i++){
            int v = f2v[f][i];
            mean_vertex_positions.row(v) += detached_faces[f][i];
        }
    }
    for(int v = 0; v < polyhedron.v_count; v++){
        mean_vertex_positions.row(v) /= (double)v2f[v].size();
    }

    //after doing that, we can recompute the normals by using polyhedron.set_pos
    polyhedron.set_pos(mean_vertex_positions);
    const auto& f2n = polyhedron.get_f2n();

    const auto& f2c = polyhedron.get_f2c();
    double sum_error = 0;

    for(int i = 0; i < polyhedron.f_count; i++){
        const auto& face = f2v[i];
        auto& detached_face = detached_faces[i];
        vec normal = f2n.row(i);
        vec center = f2c.row(i);
        for(int v = 0; v < face.size(); v++){
            vec cpos = mean_vertex_positions.row(face[v]);
            double corr = (cpos - center).dot(normal);
            sum_error += corr < 0 ? -corr : corr;
            vec projected_pos = cpos - corr * normal;
            detached_face[v] = projected_pos;
        }
    }
    cout << "sum_error: " << sum_error << endl;
}