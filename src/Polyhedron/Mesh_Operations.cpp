#include "Polyhedron.h"
#include "ChiralMap.h"
#include "utils.h"
#include <stdexcept>
#include <cassert>

Polyhedron Polyhedron::kis_selected(vector<int> faces){
    vector<vector<int>> r_f2v;
    auto& f2v = get_sorted_f2v();
    int r_v_count = v_count;
    if(faces.size() == 0){
    //if no faces are given, add all faces
        for(int i = 0; i < f_count; i++){
            faces.push_back(i);
        }
    }
    
    for(int f : faces){
        const int c_f_count = f2v[f].size();
        int v0 = r_v_count++;
        int v1 = f2v[f][0];
        for(int i = 0; i < c_f_count; i++){
            int v2 = f2v[f][(i + 1) % c_f_count];
            r_f2v.push_back({v0,v1,v2});
            v1 = v2;
        }
    }

    //we need to add the faces that haven't been altered.
    for(int f = 0; f < f_count; f++){
        bool found = false;
        for(int f2 : faces){
            if(f == f2){
                found = true;
                break;
            }
        }
        if(!found){
            r_f2v.push_back(f2v[f]);
        }
         //if(find(faces.begin(),faces.end(),f) == faces.end()){
        //}
    }

    auto ret = Polyhedron(r_f2v);
    ret.set_f2v_sorted(true);
    //place the vertices at the barycenter if the positions are initialized
    normalize();
    auto& f2c = get_f2c();
    Mx new_pos = Mx::Zero(r_v_count,3);
    //each row is a 3d position
    //and we have r_v_count rows
    //f2c is a matrix of size 3 x f_count and each row is a face baricenter
    new_pos.block(0,0,v_count,3) = pos_;
    int i = v_count;
    for(int f : faces){
        new_pos.row(i++) = f2c.row(f);
    }
    ret.set_pos(new_pos);
    return ret;
}

Polyhedron Polyhedron::dual(){
    auto ret = Polyhedron(get_sorted_v2f());
    ret.set_f2v_sorted(true);
    auto f2c = get_f2c();
    ret.set_pos(f2c);
    ret.normalize();
    return ret;
}

Polyhedron Polyhedron::kis(int degree){
    auto faces = select_faces(degree);
    return kis_selected(faces);
}

Polyhedron Polyhedron::trunc_selected(vector<int> vertices){
    return this->dual().kis_selected(vertices).dual();
}

Polyhedron Polyhedron::trunc(int degree){
    auto vertices = select_vertices(degree);
    return trunc_selected(vertices);
}

Polyhedron Polyhedron::join(){
    vector<vector<int>> r_f2v;
    const Edges& edges = get_edges();

    for(auto& edge : edges){
        r_f2v.push_back({edge.v1,edge.f1 + v_count,edge.v2,edge.f2 + v_count});
    }

    auto ret = Polyhedron(r_f2v);

    normalize();

    auto& f2c = get_f2c();
    Mx new_pos = Mx::Zero(v_count + f_count,3);
    new_pos.block(0,0,v_count,3) = pos_;
    new_pos.block(v_count,0,f_count,3) = f2c;
    ret.set_pos(new_pos);
    ret.set_f2v_sorted(true); //Todo: when join is broken, comment this line & try again.
    return ret;
}

Polyhedron Polyhedron::ambo(){
    auto ret = this->join().dual(); //or dual().join() ?
    return ret;
}

Polyhedron Polyhedron::gyro(){
    fix_winding();
    //gyro adds a vertex to each face, and 2 vertices to each edge.
    //then it adds edges between the center vertex and the other new vertices.
    //but not to every new vertex. it is a chiral operator. per face, it only adds one edge per original edge.
    //and not 2 edges per original edge, even though there are 2 new vertices per original edge.
    //and it is always the first vertex of the pair of new vertices. the second one is never connected, it is
    //connected on the other face, because the winding makes it the first vertex for this edge of the other face.
    //we need to compute the new faces also.
    //we use typedef unordered_map<tuple<int,int,int>,vector<int>,KeyHash3,KeyEqual3> EdgeMap3;
    //first we create a map from edges to new vertices:

    int new_v_it = v_count;

    auto& edges = get_edges();

    auto& f2v = get_sorted_f2v();

    int new_v_count = v_count + f_count + edges.size() * 2;

    normalize();

    Mx new_pos = Mx::Zero(new_v_count,3);
    new_pos.block(0,0,v_count,3) = pos_;

    ChiralMap new_vertices;

    for(auto& edge : edges) {
        const int v1 = edge.v1;
        const int v2 = edge.v2;

        const int nv1 = new_v_it++;
        const int nv2 = new_v_it++;

        new_vertices.emplace(ChiralKey(v1,v2,v1),nv1);
        new_vertices.emplace(ChiralKey(v1,v2,v2),nv2);

        new_pos.row(nv1) = (pos_.row(v1) * 2 + pos_.row(v2)) / 3;
        new_pos.row(nv2) = (pos_.row(v1) + pos_.row(v2) * 2) / 3;
    }

    //now we create the new faces:
    vector<vector<int>> r_f2v;
    for(int f = 0; f < f_count; f++){
        const int c_f_count = f2v[f].size();
        int v0 = new_v_it++; //the center vertex. they make up the last f_count vertices of the new Polyhedron.
        for(int i = 0; i < c_f_count; i++){
            int v1 = f2v[f][i];
            int v2 = f2v[f][(i + 1) % c_f_count];
            int v3 = f2v[f][(i + 2) % c_f_count];

            //might be an issue with the retrieval?
            int nv1 = new_vertices[ChiralKey(v1,v2,v1)];
            int nv2 = new_vertices[ChiralKey(v1,v2,v2)];
            int nv3 = new_vertices[ChiralKey(v2,v3,v2)];

            r_f2v.push_back({v0,nv1,nv2,v2,nv3});
        }
    }
    auto ret = Polyhedron(r_f2v);
    ret.set_f2v_sorted(true); //this could be the culprit - the new sorting algorithm might not be working properly.

    //add center positions last
    auto& f2c = get_f2c();

    for(int f = 0; f < f_count; f++){
        int v0 = new_v_count - f_count + f;
        new_pos.row(v0) = f2c.row(f);
    }

    ret.set_pos(new_pos);
    ret.normalize();
    return ret;
}

Polyhedron Polyhedron::propeller(){
    fix_winding();
    //gyro adds a vertex to each face, and 2 vertices to each edge.
    //then it adds edges between the center vertex and the other new vertices.
    //but not to every new vertex. it is a chiral operator. per face, it only adds one edge per original edge.
    //and not 2 edges per original edge, even though there are 2 new vertices per original edge.
    //and it is always the first vertex of the pair of new vertices. the second one is never connected, it is
    //connected on the other face, because the winding makes it the first vertex for this edge of the other face.
    //we need to compute the new faces also.
    //we use typedef unordered_map<tuple<int,int,int>,vector<int>,KeyHash3,KeyEqual3> EdgeMap3;
    //first we create a map from edges to new vertices:
    int new_v_it = v_count;

    auto& edges = get_edges();
    auto& f2v = get_sorted_f2v();
    int new_v_count = v_count + edges.size() * 2;

    normalize();
    Mx new_pos = Mx::Zero(new_v_count,3);
    new_pos.block(0,0,v_count,3) = pos_;

    ChiralMap new_vertices;

    for(auto& edge : edges) {
        const int v1 = edge.v1;
        const int v2 = edge.v2;

        const int nv1 = new_v_it++;
        const int nv2 = new_v_it++;

        new_vertices.emplace(ChiralKey(v1,v2,v1),nv1);
        new_vertices.emplace(ChiralKey(v1,v2,v2),nv2);

        new_pos.row(nv1) = (pos_.row(v1) * 2 + pos_.row(v2)) / 3;
        new_pos.row(nv2) = (pos_.row(v1) + pos_.row(v2) * 2) / 3;
    }

    //now we create the new faces:
    vector<vector<int>> r_f2v;
    for(int f = 0; f < f_count; f++){
        vector<int> middle_face;
        const int c_f_count = f2v[f].size();
        for(int i = 0; i < c_f_count; i++){
            int v1 = f2v[f][i];
            int v2 = f2v[f][(i + 1) % c_f_count];
            int v3 = f2v[f][(i + 2) % c_f_count];
            int nv1 = new_vertices[ChiralKey(v1,v2,v1)];
            int nv2 = new_vertices[ChiralKey(v1,v2,v2)];
            int nv3 = new_vertices[ChiralKey(v2,v3,v2)];
            r_f2v.push_back({nv1,nv2,v2,nv3});
            middle_face.push_back(nv3);
        }
        r_f2v.push_back(middle_face);
    }
    auto ret = Polyhedron(r_f2v);
    ret.set_pos(new_pos);
    ret.normalize();
    ret.set_f2v_sorted(true);
    return ret;
}

Polyhedron Polyhedron::triangulate(vector<int> face_ids, string mode){
    auto& f2v = get_sorted_f2v();
    if(face_ids.size() == 0){
        for(int i = 0; i < f_count; i++){
            if(f2v[i].size() > 3){
                face_ids.push_back(i);
            }
        }
    }
    if(mode == "center"){
        return kis_selected(face_ids);
    }
    else if(mode == "ear"){ //i think broken.
        //cout << "ear clipping " << endl;
        //cout << "v_count: " << v_count << endl;
        vector<vector<int>> r_f2v;
        for(int f : face_ids){
            int c_f_count = f2v[f].size();
            //we need to check if the face is a triangle
            //cout << "edge count: " << c_f_count << endl;
            if(c_f_count == 3){
                continue;
            }
            //ear clipping. create new faces from v0
            int v0 = f2v[f][0];
            for(int i = 1; i < c_f_count - 1; i++){
                int v1 = f2v[f][i];
                int v2 = f2v[f][i+1];
                r_f2v.push_back({v0,v1,v2});
                //cout << "new face: " << v0 << " " << v1 << " " << v2 << endl;
            }
        }
        //add triangles to the new Polyhedron
        for(int f = 0; f < f_count; f++){
            if(f2v[f].size() == 3){
                r_f2v.push_back(f2v[f]);
            }
        }
        //cout << "original face count: " << f_count << endl;
        //cout << "new face count: " << r_f2v.size() << endl;

        //cout << "done" << endl;
        auto ret = Polyhedron(r_f2v);
        //cout << "ret v_count: " << ret.v_count << endl;
        ret.set_pos(pos_);
        //cout << "returning Polyhedron" << endl;
        ret.set_f2v_sorted(true);
        return ret;
    }
    return *this;
}


Polyhedron Polyhedron::scrape_off(double t){
    if(t < 0 || t > 1){
        cout << "t must be between 0 and 1" << endl;
        return *this;
    }

    vector<vector<int>> r_f2v; //the new faces

    if(t <= 0.5){ // => truncate. the positions are then found by interpolation using t.
        //every original face is copied.
        //every vertex gets replaced by v2f[v].size() vertices
        //the new vertices form a new face that is added to r_f2v
        //each original face that is part of the v2f gets 2 vertices added: [n, n+1 % v2f[v].size()]
        //maybe we need v2e? then we could
        auto& v2e = get_v2e();
        log("scrape_off: v2e got");
        //cout << "v2e got" << endl;
        //the first f_count entries are the original faces
        r_f2v.resize(f_count, vector<int>());
        int new_v_count = 0;
        for(int v = 0; v < v_count; v++){
            new_v_count += v2e[v].size();
        }
        log("scrape_off: new_v_count: " + to_string(new_v_count));

        Mx new_pos = Mx::Zero(new_v_count, 3);
        int v_it = -1;
        for(int v = 0; v < v_count; v++){
            vector<int> new_face;
            for(auto& e : v2e[v]){
                v_it++;
                r_f2v[e.f1].push_back(v_it);
                r_f2v[e.f2].push_back(v_it);
                new_face.push_back(v_it);
                new_pos.row(v_it) = (1.0-t)*pos_.row(v) + t*pos_.row(e.v2);
            }
            r_f2v.push_back(new_face);
        }
        log("new faces added");
        auto ret = Polyhedron(r_f2v);
        ret.set_pos(new_pos);
        return ret;
    }else{
        //almost the same. but we create 2 new vertices per edge.
        auto& edges = get_edges();
        log("scrape_off: edges got");
        const Mx& f2c = get_f2c();

        //the first f_count entries are still the original faces
        //after that, the next v_count entries are the scraping faces
        r_f2v.resize(f_count + v_count, vector<int>());
        int new_v_count = 0;
        for(auto& e : edges){
            new_v_count += 2;
        }
        Mx new_pos = Mx::Zero(new_v_count, 3);
        int v_it = 0; //here we start at 0
        for(auto& e : edges){
            int v1 = v_it++; //v_it is incremented after the assignment, which is what we want
            int v2 = v_it++;

            r_f2v[e.f1].push_back(v1);
            r_f2v[f_count + e.v1].push_back(v1);
            r_f2v[f_count + e.v2].push_back(v1);

            r_f2v[e.f2].push_back(v2);
            r_f2v[f_count + e.v1].push_back(v2);
            r_f2v[f_count + e.v2].push_back(v2);

            vec edge_center = (pos_.row(e.v1) + pos_.row(e.v2))/2.0;

            //pos is interpolated between edge_center and the center of the face.
            //should be at edge center for t = 0.5
            //and at face center for t = 1.0
            double t2 = 2.0*(t-0.5);

            new_pos.row(v1) = (1.0 - t2) * edge_center;// + t2 * f2c.row(e.f1);
            new_pos.row(v1) += t2 * f2c.row(e.f1);
            new_pos.row(v2) = (1.0 - t2) * edge_center;// + t2 * f2c.row(e.f2);
            new_pos.row(v2) += t2 * f2c.row(e.f2);

        }
        log("new faces added");
        auto ret = Polyhedron(r_f2v);
        ret.set_pos(new_pos);
        return ret;
    }
}

//just detach the faces and replace them with prisms or pyramids
//two factors: how much each face shrinks towards its own center
//and how much each face shrinks towards the center of the polyhedron
//scaling and moving are good names for these factors
//moving is the factor for the center of the polyhedron
//scaling is the factor for the center of the face
//oh and also the thickness of the resulting prisms
//maybe make it possible to select antiprisms instead of prisms and also pyramids
Polyhedron Polyhedron::detach_faces(double scaling, double moving, double thickness, bool pyramidal){
    //return *this;
    fix_winding();
    auto& f2n = get_f2n();
    auto& f2v = get_f2v();
    auto& f2c = get_f2c();

    //we first construct the vertex positions for the midplane of each polyhedron
    //they are already scaled and moved
    vector<Mx> midplanes; //each entry is a (f2v[f].size(), 3) matrix
    vector<vec> thickness_dirs; //each entry is a vector of length 3
    //vector<Mx> midplanes_normals; //each entry is a (f2v[f].size(), 3) matrix, we dont need that, thats just f2n
    //vector<Mx> midplanes_centers; //each entry is a (1, 3) matrix, we dont need that, thats just f2c
    for(int f = 0; f < f_count; f++){
        Mx midplane = Mx::Zero(f2v[f].size(), 3);
        vec normal = f2n.row(f).normalized();
        vec face_center = f2c.row(f);
        for(int i = 0; i < f2v[f].size(); i++){
            vec pos = pos_.row(f2v[f][i]);
            //project pos onto the plane defined by face_center and normal
            vec pos_proj = pos - (pos - face_center).dot(normal) * normal;
            //now we have the position of the vertex projected onto the plane. now we scale it towards the face_center of the face
            //and also move it towards the face_center of the polyhedron
            vec force_to_scale = face_center - pos_proj;
            vec force_to_move = normal;
            midplane.row(i) = pos + scaling * force_to_scale + moving * force_to_move;
        }
        thickness_dirs.push_back(thickness * normal * 0.5);
        midplanes.push_back(midplane);
    }

    int new_v_count = 0;
    vector<vector<int>> ret_f2v;
    vector<Mx> new_pos_array; //easier to first save them in an array and then copy them to a matrix

    for(int f = 0; f < f_count; f++){
        //per face, we make 2 + f2v[f].size() new faces. the first 2 are the top and bottom faces
        //or just a pyramid if pyramidal is true

        //rewrite the position finding code. it is just messed up.
        //the vertices seem to be alright. the positions are not.

        if(pyramidal){
            int v0 = -1;
            new_v_count++;
            Mx new_pos = Mx::Zero(1, 3);
            new_pos.row(0) = f2c.row(f);
            new_pos.row(0) += thickness_dirs[f];

            new_pos_array.push_back(new_pos);

            Mx new_pos2 = midplanes[f];
            new_pos2 -= thickness_dirs[f];

            new_pos_array.push_back(new_pos2);

            vector<int> bottom_face;
            const int start = new_v_count; //start - 1 is the bottom face.
            const int face_size = f2v[f].size();
            int v1 = 0;
            for(int i = 0; i < face_size; i++){
                bottom_face.push_back(v1 + start);
                int v2 = (i + 1) % face_size; //okay for the last one we land back on new_v_count. :)
                ret_f2v.push_back(vector<int>{v0 + start, v1 + start, v2 + start});
                v1 = v2;
            }
            new_v_count += face_size;
            ret_f2v.push_back(bottom_face);
            //ah we go through by face. so we indeed screw up the indices in the matrices.
        }else{
            //we have a bottom face and a top face
            //and instead of triangular side faces we have quadrilateral side faces
            //so no fixed v0
            Mx new_pos = midplanes[f];
            new_pos += thickness_dirs[f];
            new_pos_array.push_back(new_pos);

            Mx new_pos2 = midplanes[f];
            new_pos2 -= thickness_dirs[f];
            new_pos_array.push_back(new_pos2);

            const int start = new_v_count;
            vector<int> bottom_face;
            vector<int> top_face;
            const int face_size = f2v[f].size();
            assert(face_size == midplanes[f].rows());
            for(int i = 0; i < face_size; i++){
                //we add start later. first think locally. we have vertices from 0 to (2*f2v[f].size()) - 1.
                //we want to add start to them later.
                int v0 = i;
                int v1 = (i + 1) % face_size;
                int v2 = face_size + v0; //here, we increase by face_size. which is the number of vertices in the face
                int v3 = face_size + v1; //it should be the same as the number of rows in midplanes[f].
                ret_f2v.push_back(vector<int>{v0 + start, v1 + start, v3 + start, v2 + start});
                bottom_face.push_back(v0 + start);
                top_face.push_back(v2 + start);
            }
            new_v_count += 2 * face_size;
            ret_f2v.push_back(bottom_face);
            ret_f2v.push_back(top_face);
        }
    }
    //now we copy the new positions to a matrix and also check wether new_v_count is correct
    Mx ret_pos = Mx::Zero(new_v_count, 3);
    int i = 0;
    //we have to copy the new positions to the matrix
    for(auto& pos : new_pos_array){
        ret_pos.block(i, 0, pos.rows(), pos.cols()) = pos;
        i += pos.rows();
    }
    assert(i == new_v_count);

    //now construct the new mesh
    Polyhedron ret(ret_f2v);
    ret.set_pos(ret_pos);
    return ret;
}

