#include "Polyhedron.h"
#include "utils.h"
#include <queue>

Mx Polyhedron::get_pos(){
    return pos_;
}

void Polyhedron::set_pos(Mx& pos){
    /*
    sets the position of the vertices. The matrix must have at least as many rows as vertices and exactly 3 columns.
    */
    if(pos.rows() < v_count){
        cout << "wrong number of vertices in position matrix" << endl;
        return;
    }
    if(pos.cols() != 3){
        cout << "wrong number of columns in position matrix" << endl;
        return;
    }
    pos_ = pos.block(0,0,v_count,3);
    position_changed();
}

void Polyhedron::rotate(double angle, int axis){
    //inline Eigen::Matrix3d rotationMatrix(int axis, double eps) {

    pos_ *= rotationMatrix(axis, angle);
    position_changed();
}

void Polyhedron::position_changed(){
    /*
    this function needs to be called after pos has changed,
    so that the face normals and centers get recomputed.
    since get_pos returns a reference, this function needs to be called manually after changing the position.
    */
    f2c_computed = false;
    f2n_computed = false;
    center_computed = false;
}

const Mx& Polyhedron::get_f2c(){
    /*
    this function computes the midpoints of the faces by averaging the
    vertices of the faces.
    */
    if(!f2c_computed){
        auto& f2v = get_f2v();
        for(int f = 0; f < f_count; f++){
            //face center calculation
            vec center(0,0,0);
            for(int v : f2v[f]){
                center += pos_.row(v);
            }
            center /= double(f2v[f].size());
            f2c_.row(f) = center;
        }
        f2c_computed = true;
    }
    return f2c_;
}

const Mx& Polyhedron::get_f2n(){
    /*
    this function computes the normals of the faces of the mesh.
    it does so by solving the total least squares problem of fitting a plane to the vertices of the face,
    such that the orthogonal distances of the vertices to the plane are minimized.
    */
    if(!f2n_computed){
        //compute normals of faces based on f2v pos
        auto& f2c = get_f2c();
        auto& f2v = get_f2v();
        for(int f = 0; f < f_count; f++){
            //make the matrix containing the positions
            Mx position_matrix(f2v[f].size(), 3);
            int i = 0;
            for(int v : f2v[f]){
                position_matrix.row(i++) = pos_.row(v) - f2c.row(f);
            }

            //compute the SVD of the position matrix
            Eigen::JacobiSVD<Mx> svd(position_matrix, Eigen::ComputeFullV);

            //the last column of V contains the normal vector
            f2n_.row(f) = svd.matrixV().col(2).transpose().normalized();
        }
        f2n_computed = true;
    }
    return f2n_;
}

const vec& Polyhedron::get_center(){
    if(!center_computed){
        center_ = pos_.colwise().mean();
    }
    return center_;
}

void Polyhedron::to_origin(){
    /*
        moves the average position to the origin
    */
    pos_.rowwise() -= get_center().transpose();
    position_changed();
}

void Polyhedron::scale(double factor){
    to_origin();
    //make maxCoeff distance to origin 1
    //originally was mean distance to origin, but for rendering purposes, this is better.
    pos_ *= factor / pos_.rowwise().norm().maxCoeff();
    position_changed();
}

void Polyhedron::unskew(){
    /*
    uses gram_schmidt to unskew the 3 dimensions. moves to origin but preserves scale
    so that it can be used in conjunction with canonical force application.
    because canonical force changes the scale.
    since we use differing definitions of "scale" for canonical() and scale().
    scale() uses average vertex distance to center.
    canonical() uses minimal edge distance to center.
    */

    to_origin();

    const double before_scale = pos_.rowwise().norm().maxCoeff();

    pos_ /= before_scale;

    pos_.col(0).normalize();

    const double sc1 = pos_.col(0).dot(pos_.col(1));
    const double sc2 = pos_.col(0).dot(pos_.col(2));

    pos_.col(1) -= pos_.col(0) * sc1;
    pos_.col(1).normalize();
    const double sc3 = pos_.col(1).dot(pos_.col(2));

    pos_.col(2) -= pos_.col(0) * sc2 + pos_.col(1) * sc3;
    pos_.col(2).normalize();

    to_origin(); //although this shouldnt be necessary, mabye it is beneficial.

    pos_ *= before_scale / pos_.rowwise().norm().maxCoeff(); //rescaling to original size

    position_changed();
}

void Polyhedron::normalize(){
    /*
    this function normalizes the positions of the vertices.
    it does so by subtracting the mean of each dimension,
    decorrelating the 3 dimensions, so that skewed shapes are fixed,
    and then normalizing the mean norm of the vertices to 1.
    for decorrelation, it uses gram-schmidt orthogonalization.
    */

    scale(1.0);

    pos_.col(0).normalize();

    const double sc1 = pos_.col(0).dot(pos_.col(1));
    const double sc2 = pos_.col(0).dot(pos_.col(2));

    pos_.col(1) -= pos_.col(0) * sc1;
    pos_.col(1).normalize();
    const double sc3 = pos_.col(1).dot(pos_.col(2));

    pos_.col(2) -= pos_.col(0) * sc2 + pos_.col(1) * sc3;
    pos_.col(2).normalize();

    scale(1.0);

    position_changed();
}

vector<int> Polyhedron::select_indices(vector<int> indices, VectorXd D, Mx V, int mat_size){
    //select the indices of the biggest eigenvalues
    vector<int> ret;
    if(indices.size() != 3){
        //cout << "using standard indices for spectral realization" << endl;
        ret.push_back(mat_size-2);
        ret.push_back(mat_size-3);
        ret.push_back(mat_size-4);
        //TODO: select indices based on D
    }else{
        //wrap indices
        for(int i = 0; i < 3; i++){
            int index = indices[i];
            while(index < 0){
                index += mat_size;
            }
            index %= mat_size;
            ret.push_back(index);
        }
    }
    return ret;
}

void Polyhedron::realize_spectral(Mx mat, vector<int> indices){
    if((mat.rows() < v_count) || (mat.cols() < v_count)){
        cout << "No valid matrix provided for spectral realization. Using face mean matrix." << endl;
        //default case: adjacency matrix
        mat = get_matrix("face_mean");
    }
    Eigen::SelfAdjointEigenSolver<Mx> es(mat);
    Mx V = es.eigenvectors();
    VectorXd D = es.eigenvalues();

    indices = select_indices(indices,D,V,mat.cols());

    pos_.col(0) = V.col(indices[0]).head(v_count);
    pos_.col(1) = V.col(indices[1]).head(v_count);
    pos_.col(2) = V.col(indices[2]).head(v_count);

    position_changed();
}

void Polyhedron::hang_axis(int axis){
    /*
    this function does the following:
    1. find the vertex with the smallest distance to the origin
    2. find the face with the smallest distance to the origin
    3. if the face is closer, triangulate it and call hang_axis on the triangulation
    4. if the vertex is closer, set the vertex to 0 and set all other vertices to their distance to the vertex
        along the specified axis
    this function fixes the hyperbolic component of the singular value decomposition if it occurs.
    it doesnt always produce convex meshes, but it does produce meshes with a good realization.
    */
    double min_vertex_distance = INFINITY;
    double min_face_distance = INFINITY;
    int v0 = 0;
    int f0 = 0;


    //set to zero
    pos_.col(axis) = Mx::Zero(v_count,1);
    position_changed();
    auto& centers = get_f2c();
    for(int v = 0; v < v_count; v++){
        //get closest point to origin
        double dist = pos_.row(v).norm();
        if(dist <= min_vertex_distance){
            min_vertex_distance = dist;
            v0 = v;
        }
    }
    for(int f = 0; f < f_count; f++){
        //get closest point to origin
        double dist = centers.row(f).norm();
        if(dist < min_face_distance){
            min_face_distance = dist;
            f0 = f;
        }
    }
    if(min_face_distance < min_vertex_distance){
    //we hang the fac
        auto t = triangulate({f0},"center");
        t.hang_axis(axis);
        auto tpos = t.get_pos();
        pos_.col(axis) = tpos.col(axis).head(v_count);
        return;
    }

    pos_(v0,axis) = 0;
    auto& v2v = get_v2v();

    vector<int> distances(v_count, -1); // store the distance from v0 to each vertex
    distances[v0] = 0; // distance from v0 to v0 is 0
    queue<int> q;
    q.push(v0);

    while (!q.empty()) {
        int v = q.front();
        q.pop();

        for (int n : v2v[v]) {
            if (distances[n] == -1) {
                distances[n] = distances[v] + 1;
                q.push(n);
                pos_(n, axis) = distances[n];
            }
        }
    }
    position_changed();
}

Mx Polyhedron::get_planar_drawing(int v0,int axis){
 /*
 finds a planar drawing by first spectral realization, then the same as hang
 but not assigning the third axis but instead the radius from the origin for each vertex
 the return type is positions for each vertex.
 */
    v0 = fix_vertex_index(v0); //wrap around, be able to use -1 and so on.
    pos_(v0,axis) = 0;
    auto& v2v = get_v2v();
    vector<int> distances(v_count, -1); // store the distance from v0 to each vertex
    distances[v0] = 0; // distance from v0 to v0 is 0
    queue<int> q;
    q.push(v0);
    double max_dist = 0;
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        vector<int> cv = v2v[v];
        for (const int n : cv) {
            if (distances[n] == -1) {
                distances[n] = distances[v] + 1;
                q.push(n);
                pos_(n,axis) = 0;
                if(pos_.row(n).norm() > 0){
                    pos_.row(n).normalize();
                }
                pos_.row(n) *= distances[n];
                if(distances[n] > max_dist){
                    max_dist = distances[n];
                }
            }
        }
    }

    //normalize the distances to be maximally 1
    pos_ /= max_dist;
    //try later: scale(1.0); //should work also.
    position_changed();
    return pos_;
}

void Polyhedron::spherical(){
    /*
    set all vertices to be on unit sphere
    */
    pos_.rowwise().normalize();
    position_changed();
}

void Polyhedron::fix_insideout(){
    /*
    this is an addon to winding fixing. while winding fixing calculates the correct winding, it
    doesn't look at the vertex positions. so the whole mesh could be inside out.
    this function checks if the mesh is inside out and flips it if necessary, so that the
    f2v winding adheres to obj file format conventions.
    it does so by first normalizing the mesh, then checking if it gets smaller or bigger when moving
    along cross product computed normals. it only looks at the first 3 vertices of each face to compute
    the normals, it should be alright.
    */
    fix_winding();
    //now the average distance to 0 is 1.
    //we can use this to check if we are inside or outside
    //we calculate the normal by cross product of the first 3 vertices
    double avg_dist_1 = 0;
    double avg_dist_2 = 0;
    vec center = get_center();
    auto& f2v = get_sorted_f2v(); //const reference!
    for(int f = 0; f < f_count; f++){
        const int v1 = f2v[f][0];
        const int v2 = f2v[f][1];
        const int v3 = f2v[f][2];
        Vector3d p1 = pos_.row(v1);
        Vector3d p2 = pos_.row(v2);
        Vector3d p3 = pos_.row(v3);
        const auto n = (p2-p1).cross(p3-p1).normalized();
        avg_dist_1 += (p1 + n*0.1 - center).norm();
        avg_dist_2 += (p1 - n*0.1 - center).norm();
        //avg_dist_1 += (p1 + n*0.1).norm();
        //avg_dist_2 += (p1 - n*0.1).norm();
    }
    ////cout << "avg dist 1 " << avg_dist_1 << " avg dist 2 " << avg_dist_2 << endl;
    if(avg_dist_1 < avg_dist_2){
        //cout << "flipping all faces" << endl;
        for(int f = 0; f < f_count; f++){
            std::reverse(f2v_[f].begin(), f2v_[f].end());
        }
    }
}

double Polyhedron::align(int steps) {
    to_origin();

    //we align by rotating +-eps and look at the projected vertex distances.
    //we can decay eps

    double decay = 0.9;
    double min_eps = 0.00001;

    double eps = 3.1415926535 * 1.9; //radians

    double delta = 0;

    for (int i = 0; i < steps; i++){
        for(int axis = 0; axis < 3; axis++){
            //compute rotation matrices +-eps
            Matrix3d rot1 = rotationMatrix(axis, eps);
            Matrix3d rot2 = rotationMatrix(axis, -eps);
            //rotate
            MatrixXd pos1 = pos_ * rot1;
            MatrixXd pos2 = pos_ * rot2;
            //compute distances. leave out the current axis
            double dist1 = 0;
            double dist2 = 0;
            for(int v = 0; v < v_count; v++){
                for(int j = 0; j < 3; j++){
                    if(j != axis){
                        dist1 += abs(pos1(v,j));
                        dist2 += abs(pos2(v,j));
                    }
                }
            }
            delta = abs(dist1 - dist2);
            if(dist1 > dist2){
                pos_ = pos1;
            } else {
                pos_ = pos2;
            }
        }
        eps *= decay;
        if(eps < min_eps){
            eps = min_eps;
        }
    }
    return delta;
}


/*
void Polyhedron::crystallize(double min_angle){
    //look at the face normals and push them apart if they are too close
    auto& f2n = get_f2n();
    auto& f2c = get_f2c();
    auto& v2f = get_v2f();
    auto& pos = get_pos();
    Mx normal_deltas(f_count,3);
    normal_deltas.setZero();
    auto& edge_map = get_edge_map();

    for(const auto& [vertices, faces] : edge_map) {
        const int f1 = faces[0];
        const int f2 = faces[1];
        const auto& n1 = f2n.row(f1);
        const auto& n2 = f2n.row(f2);
        if(acos(n1.dot(n2)) <= min_angle){
            //the faces are too coplanar.
            //we push them apart by pushing apart the normals.
            const auto& c1 = f2c.row(f1);
            const auto& c2 = f2c.row(f2);
            const Vector3d delta = (c2 - c1).normalized();
            normal_deltas.row(f1) -= delta * 0.01;
            normal_deltas.row(f2) += delta * 0.01;
        }
    }{
        const int f1 = faces[0];
        const int f2 = faces[1];
        const auto& n1 = f2n.row(f1);
        const auto& n2 = f2n.row(f2);
        if(acos(n1.dot(n2)) <= min_angle){
            //the faces are too coplanar.
            //we push them apart by pushing apart the normals.
            const auto& c1 = f2c.row(f1);
            const auto& c2 = f2c.row(f2);
            const Vector3d delta = (c2 - c1).normalized();
            normal_deltas.row(f1) -= delta * 0.01;
            normal_deltas.row(f2) += delta * 0.01;
        }
    }
    Mx pos_deltas(v_count,3);
    pos_deltas.setZero();

    for(int f = 0; f < f_count; f++){
        auto normal = (f2n.row(f) + normal_deltas.row(f)).normalized();
        auto center = f2c.row(f);
        for(int v : f2v[f]){
            pos_deltas.row(v) += normal * normal.dot(pos.row(v) - center) / double(v2f[v].size());
        }
    }
    //crystallizer is broken right now. might be wrong direction to move the vertices in.
    pos += pos_deltas;
    position_changed();
}
*/
