#include "Polyhedron.h"
#include "Edge.h"
#include "Optimizer.h"
#include "Renderer.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <time.h>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>    // std::sort, reverse
#include <stdexcept>
#include <map>
#include <tuple>
#include <filesystem>

using namespace std;

Polyhedron::Polyhedron(vector<vector<int>> f2v){
    // log("creating polyhedron");
    f_count = f2v.size();
    // log("f_count = " + to_string(f_count));
    v_count = get_v_count(f2v);
    // log("v_count = " + to_string(v_count));
    f2v_ = f2v;
    reset();
    // log("polyhedron creation complete");
}

void Polyhedron::reset(){
    f2v_sorted = false;
    v2f_sorted = false;
    winding_fixed = false;
    f2c_computed = false;
    f2n_computed = false;
    force_initialized = false;
    pos_ = Mx::Zero(v_count,3);
    f2n_ = Mx::Zero(f_count,3);
    f2c_ = Mx::Zero(f_count,3);
}

// void Polyhedron::log(string s){
//     if(log_){
//         cout << s << endl;
//     }
// }

void Polyhedron::check(){
    // log("checking polyhedron");
    if(v_count < 4){
        //cout << "v_count < 4" << endl;
        //cout << string_repr() << endl;
        throw invalid_argument("Polyhedron check failed, v_count < 4");
    }
    // log("v_count >= 4");
    if(f_count < 4){
        //cout << "f_count < 4" << endl;
        //cout << string_repr() << endl;
        throw invalid_argument("Polyhedron check failed, f_count < 4");
    }
    // log("f_count >= 4");
    // log("getting f2v");
    //first some f2v checks ----------------------------------------------------
    auto& f2v = get_f2v();
    if(f2v.size() != f_count){
        //cout << "f2v.size() != f_count" << endl;
        //cout << string_repr() << endl;
        throw invalid_argument("Polyhedron check failed, f2v.size() != f_count");
    }
    //log("f2v.size() == f_count");
    //check vertex usage first
    //log("checking vertex usage");
    vector<bool> used(v_count, false);
    for(int f = 0; f < f_count; f++){
        for(int v : f2v[f]){
            used[v] = true;
        }
    }
    for(int v = 0; v < v_count; v++){
        if(!used[v]){
            //cout << "vertex " << v << " is not used" << endl;
            //cout << string_repr() << endl;
            throw invalid_argument("Polyhedron check failed, a vertex is unused");
        }
    }
    //log("vertex usage ok");
    //log("checking face vertex count and duplicates");
     for(int f = 0; f < f_count; f++){
        if(f2v[f].size() < 3){
            //cout << "face " << f << " has less than 3 vertices" << endl;
            //cout << string_repr() << endl;
            throw invalid_argument("Polyhedron check failed, a face has less than 3 vertices");
        }
    }
    // log("face vertex count ok");
    // log("checking face vertex duplicates and out of range");
    for(int f = 0; f < f_count; f++){
        for(int v : f2v[f]){
            if((v >= v_count) || (v < 0)){
                // cout << "face " << f << " has vertex index out of range" << endl;
                // cout << string_repr() << endl;
                throw invalid_argument("Polyhedron check failed, a face has vertex out of range");
            }
            for(int j = v + 1; j < f2v[f].size(); j++){
                if(f2v[f][v] == f2v[f][j]){
                    // cout << "face " << f << " has duplicate vertices" << endl;
                    // cout << string_repr() << endl;
                    throw invalid_argument("Polyhedron check failed, a face has duplicate vertices");
                }
            }
        }
     }
    // log("face vertex duplicates and out of range ok");
    //then some v2f checks ----------------------------------------------------

    // log("getting v2f");
    auto& v2f = get_v2f();
    // log("got v2f");
    if(v2f.size() != v_count){
        // cout << "v2f.size() != v_count" << endl;
        // cout << string_repr() << endl;
        throw invalid_argument("Polyhedron check failed, v2f.size() != v_count");
    }
    // log("v2f.size() == v_count");
    //then some v2v checks ----------------------------------------------------
    // log("getting v2v");
    auto& v2v = get_v2v();
    // log("got v2v");
    // log("checking v2v size and 3-connectedness");
    if(v2v.size() != v_count){
        // cout << "v2v.size() != v_count" << endl;
        // cout << string_repr() << endl;
        throw invalid_argument("Polyhedron check failed, v2v.size() != v_count");
    }

    for(int v = 0; v < v_count; v++){
        if(v2v[v].size() < 3){
            cout << "vertex " << v << " is not 3-connected" << endl;
            cout << string_repr() << endl;
            throw invalid_argument("Polyhedron check failed, a vertex is not 3-connected");
        }
    }
    // log("v2v size and 3-connectedness ok");

    //check v2v euler formula
    int v2v_count = 0;
    for(int v = 0; v < v_count; v++){
        v2v_count += v2v[v].size();
    }
    if(v_count + f_count - v2v_count/2 != 2){
        // cout << "v2v Euler's formula failed" << endl;
        // cout << v_count << " " << f_count << " " << v2v_count << "/2" << endl;
        // cout << string_repr() << endl;
        throw invalid_argument("Polyhedron check failed, v2v Euler's formula failed");
    }
    // log("v2v euler formula ok");

    //then edge checks --------------------------------------------------------

    // log("getting edges");
    auto& edges = get_edges(); //std::set<Edge>
    // log("got edges");
    //check eulers formula
    if(v_count + f_count - edges.size() != 2){
        // cout << "Euler's formula failed" << endl;
        // cout << v_count << " " << f_count << " " << edges.size() << endl;
        // cout << string_repr() << endl;
        throw invalid_argument("Polyhedron check failed, Euler's formula failed");
    }
    // log("euler's formula ok");
    // log("check() finished.");
}

// string Polyhedron::string_repr(){
//     stringstream ss;
//     auto& f2v = get_f2v();
//     ss << "v_count: " << v_count << endl;
//     ss << "f_count: " << f_count << endl;
//     ss << "f2v: " << endl;
//     for(int f = 0; f < f_count; f++){
//         ss << f << ": ";
//         for(int v : f2v[f]){
//             ss << v << " ";
//         }
//         ss << endl;
//     }
//     return ss.str();
// }

int Polyhedron::fix_vertex_index(int ix){
    while(ix < 0){
        ix += v_count;
    }
    ix %= v_count;
    return ix;
}

void Polyhedron::set_log(bool log){
    this->log_ = log;
}

void Polyhedron::save_obj(string fname,bool overwrite = false, bool create_folders = true){
    //need to fix normalizing behaviour.
    fix_winding();
    //fix_insideout();

    if (!overwrite && filesystem::exists(fname)) {
        throw invalid_argument("File already exists");
    }

    // Check if folder exists, and create it if necessary
    filesystem::path pathObj(fname);
    filesystem::path folder = pathObj.parent_path();
    if (!filesystem::exists(folder)) {
        if (create_folders) {
            if (!filesystem::create_directories(folder)) {
                throw invalid_argument("Failed to create folder");
            }
        } else {
            throw invalid_argument("Folder does not exist");
        }
    }

    // Open file for writing
    ofstream fs(fname);
    if (!fs) {
        throw invalid_argument("Failed to open file for writing");
    }

    fs.precision(10);

    auto& f2v = get_sorted_f2v();

    for(int i = 0; i < v_count; i++){
        fs << "v " << pos_(i,0) << " " << pos_(i,1) << " " << pos_(i,2) << endl;
    }
    for(int f = 0; f < f_count; f++){
        fs << "f ";
        for(auto vertex : f2v[f]){
            fs << vertex + 1 << " ";
        }
        fs << endl;
    }
    fs.close();
}

Mx Polyhedron::get_matrix(string type){
    Mx ret = Mx::Zero(v_count, v_count);
    auto& v2v = get_v2v();
    if(type == "adjacency"){
        for(int v = 0; v < v_count; v++){
            for(int v2 : v2v[v]){
                ret(v,v2) = 1;
                ret(v2,v) = 1;
            }
        }
        return ret;
    }
    if(type == "laplacian"){
        for(int v = 0; v < v_count; v++){
            for(int v2 : v2v[v]){
                ret(v,v2) = 1;
                ret(v2,v) = 1;
            }
        }
        for(int v = 0; v < v_count; v++){
            ret(v,v) = -v2v[v].size();
        }
        return ret;
    }
    if(type == "face_mean"){
        auto& v2f = get_v2f();
        auto& f2v = get_f2v();
        for(int v = 0; v < v_count; v++){
            for(int f : v2f[v]){
                double fac = 1.0 / (v2f[v].size() * f2v[f].size());
                for(int v2 : f2v[f]){
                    ret(v,v2) += fac;
                }
            }
        }
        return ret;
    }
    if(type == "vertex_mean"){
        for(int v = 0; v < v_count; v++){
            for(int v2 : v2v[v]){
                ret(v,v2) += 1.0 / v2v[v].size();
            }
        }
        return ret;
    }
    // cout << "unknown matrix type: " << type << endl;
    // cout << "available types: adjacency, laplacian, face_mean, vertex_mean" << endl;
    // cout << "returning empty matrix" << endl;
    invalid_argument("unkown matrix type; pick one of adjacency, laplacian, face_mean, vertex_mean");
    // return ret;
}

vector<int> Polyhedron::select_vertices(int order){
    vector<int> ret;
    auto& v2f = get_v2f();
    if(order == -1){
        for(int v = 0; v < v2f.size(); v++){
            ret.push_back(v);
        }
    }else{
        for(int v = 0; v < v2f.size(); v++){
            if(v2f[v].size() == order){
                ret.push_back(v);
            }
        }
    }
    return ret;
}

vector<int> Polyhedron::select_faces(int order){
    vector<int> ret;
    auto& f2v = get_f2v();
    if(order == -1){
        for(int f = 0; f < f2v.size(); f++){
            ret.push_back(f);
        }
    }else{
        for(int f = 0; f < f2v.size(); f++){
            if(f2v[f].size() == order){
                ret.push_back(f);
            }
        }
    }
    return ret;
}

const Edges& Polyhedron::get_edges() {
    if(edges_.size() == 0){
        auto& f2v = get_sorted_f2v();
        auto& v2f = get_v2f(); //can't be sorted - would result in a call loop
        for (int i = 0; i < f2v.size(); i++) {
            const std::vector<int>& face = f2v[i];
            int num_verts = face.size();
            for (int j = 0; j < num_verts; j++) {
                int v1 = face[j];
                int v2 = face[(j+1) % num_verts];
                // This edge is an interior edge
                int f1 = i;
                int f2 = -1;
                for (int k = 0; k < v2f[v2].size(); k++) {
                    int other_face = v2f[v2][k];
                    if (other_face != f1 && std::find(f2v[other_face].begin(), f2v[other_face].end(), v1) != f2v[other_face].end()) {
                        f2 = other_face;
                        break;
                    }
                }
                if (f2 != -1) {
                    edges_.emplace(v1, v2, f1, f2);
                }
            }
        }
    }
    return edges_;
}

const vector<vector<int>>& Polyhedron::get_v2f(){
    if(v2f_.size() == 0){
        for(int v = 0; v < v_count; v++){
            v2f_.push_back({});
        }
        //fill v2f
        auto& f2v = get_f2v();
        for(int f1 = 0; f1 < f_count; f1++){
            for(int v : f2v[f1]){
                v2f_[v].push_back(f1);
            }
        }
        v2f_sorted = false;
    }
    return v2f_;
}

const vector<vector<int>>& Polyhedron::get_sorted_v2f(){
    sort_v2f();
    return get_v2f();
}

const vector<vector<int>>& Polyhedron::get_f2v(){
    return f2v_;
}

const vector<vector<int>>& Polyhedron::get_sorted_f2v(){
    sort_f2v();
    return f2v_;
}

const vector<vector<Edge>>& Polyhedron::get_f2e(){
    if(f2e_.size() == 0){
        auto& edges = get_edges();
        auto& f2v = get_f2v();
        f2e_.reserve(f_count);
        for(int i = 0; i < f_count; i++){
            f2e_.push_back({});
            f2e_.back().reserve(f2v[i].size());
        }
        for(auto& edge : edges) {
            f2e_[edge.f1].push_back(Edge(edge.v1,edge.v2,edge.f2,edge.f2));
            f2e_[edge.f2].push_back(Edge(edge.v1,edge.v2,edge.f1,edge.f1));
        }
    }
    return f2e_;
}

void Polyhedron::sort_f2v() {
    if(f2v_sorted){
        return;
    }
    auto& f2v = get_f2v();
    auto& v2f = get_v2f();
    for (int f = 0; f < f2v.size(); f++) {
        vector<int> sorted_face;
        int current_vertex = f2v[f][0];
        sorted_face.push_back(current_vertex);

        while (sorted_face.size() < f2v[f].size()) {
            int next_vertex = -1;
            for (int neighbor_face_index : v2f[current_vertex]) {
                if (neighbor_face_index == f) continue;
                const vector<int>& neighbor_face = f2v[neighbor_face_index];
                for (int v : neighbor_face) {
                    if (v != current_vertex && find(sorted_face.begin(), sorted_face.end(), v) == sorted_face.end() &&
                        find(f2v[f].begin(), f2v[f].end(), v) != f2v[f].end()) {
                        next_vertex = v;
                        break;
                    }
                }

                if (next_vertex != -1) break;
            }
            sorted_face.push_back(next_vertex);
            current_vertex = next_vertex;
        }
        f2v_[f] = sorted_face;
    }
    f2v_sorted = true;
}

void Polyhedron::sort_v2f(){
    if(v2f_sorted && (v2f_.size() > 0)){
        return;
    }
    auto& v2f = get_v2f();
    auto& f2e = get_f2e();
    vector<bool> visited;
    visited.resize(f_count,false);
    for(int v = 0; v < v_count; v++){
        vector<int> sorted;
        auto& vf = v2f[v];
        const int c_v_count = vf.size();
        int f1 = vf[0];
        sorted.reserve(c_v_count);
        sorted.push_back(f1);
        for(int f2 : vf){
            visited[f2] = false;
        }
        int f2;
        while(sorted.size() < c_v_count){
            visited[f1] = true;
            for(auto& e : f2e[f1]){ // im stumped. the
                f2 = e.f2;
                if((!visited[f2]) && ((e.v1 == v) || (e.v2 == v))){
                    f1 = f2;
                    sorted.push_back(f1);
                    break;
                }
            }
        }
        v2f_[v] = sorted;
    }
    v2f_sorted = true;
}

const vector<vector<int>>& Polyhedron::get_v2v(){
    if(v2v_.size() == 0){
        for(int i = 0; i < v_count; i++){
            v2v_.push_back({});
        }
        auto& f2v = get_sorted_f2v();
        auto& v2f = get_sorted_v2f();

        for(int v = 0; v < v_count; v++){
            auto& vf = v2f[v];
            const int c_v_count = vf.size();
            for(int i = 0; i < c_v_count; i++){
                int f1 = vf[i];
                int f2 = vf[(i+1)%c_v_count];
                for(int v1 : f2v[f1]){
                    for(int v2 : f2v[f2]){
                        if((v1 == v2) && (v1 != v)){
                            v2v_[v].push_back(v1);
                        }
                    }
                }
            }
        }
    }
    return v2v_;
}

const vector<vector<Edge>>& Polyhedron::get_v2e(){
    if(v2e_.size() == 0){
        auto& edges = get_edges();
        auto& v2f = get_v2f();
        v2e_.reserve(v_count);
        for(int i = 0; i < v_count; i++){
            v2e_.push_back({});
            v2e_.back().reserve(v2f[i].size());
        }
        for(auto& edge : edges) {
            v2e_[edge.v1].push_back(Edge(edge.v2,edge.v2,edge.f1,edge.f2));
            v2e_[edge.v2].push_back(Edge(edge.v1,edge.v1,edge.f1,edge.f2));
        }
    }
    return v2e_;
}

//now that we fixed the winding we can do gyro!!:-)
void Polyhedron::fix_winding(){
    //fixes the winding order of the faces. this is necessary for some operations and also for
    //obj export. but obj export also depends on the vertex positions, that gets fixed in a separate function
    //this function only makes sure that the edges never occur in the same order in the two faces that share them
    //this ensures that the winding order is consistent.
    if(winding_fixed){
        return;
    }
    //we need to select a random face and set that as reference
    //then we flip the winding for each face
    //at the end we try something to see wether we did it inside or outside
    vector<bool> winding_flag;
    vector<int> winding_index;
    winding_flag.resize(f_count,false);
    int current_face = 0;
    int last_face = 0;
    int v1,v2;
    //first find a new face that shares an edge
    //cout << "getting f2e " << endl;
    auto& f2e = get_f2e();
    auto& f2v = get_sorted_f2v();
    //cout << "done" << endl;
    while(true){
        winding_index.push_back(current_face);
        winding_flag[current_face] = true;
        bool found_face = false;
        //first search for a face that needs to be checked that is connected to the already checked faces
        //to do this we iterate over the corrected faces (winding_index)
        //much better to iterate in reverse. since new faces are more likely to have unfixed neighbours.
        for(int i = winding_index.size() - 1; i >= 0; i--){
            //we then check every edge
            for(auto& e : f2e[winding_index[i]]){
                //f2e[face_id] maps to a list of entries of [other_face_id,v1,v2]
                if(!winding_flag[e.f2]){
                    //face found
                    found_face = true;
                    current_face = e.f2;
                    last_face = winding_index[i];
                    v1 = e.v1;
                    v2 = e.v2;
                    //cout << "found face " << current_face << " last face " << last_face << endl;
                    break;
                }
            }
            if(found_face){
                break;
            }
        }
        if(!found_face){
            break;
        }

        const int f1_size = f2v[last_face].size();
        const int f2_size = f2v[current_face].size();

        //here we perform a swap of v1 and v2 so that they are in the right order for the winding check
        for(int j = 0; j < f1_size; j++){
            if((v1 == f2v[last_face][(j+1) % f1_size]) && (v2 == f2v[last_face][j])){
                swap(v1,v2);
                //cout << "swapped v1 and v2" << v1 << " " << v2 << endl;
                break;
            }
        }

        //the actual winding check for the current face
        for(int j = 0; j < f2_size; j++){
            int v3 = f2v[current_face][j];
            int v4 = f2v[current_face][(j+1) % f2_size];
            if((v1 == v3) && (v2 == v4)){
                //cout << "flipping face " << current_face << endl;
                reverse(f2v_[current_face].begin(), f2v_[current_face].end());
                break;
            }
        }
    }
    winding_fixed = true;
}

Polyhedron Polyhedron::flip_winding(){
    fix_winding();
    for(auto& f : f2v_){
        reverse(f.begin(), f.end());
    }
    return *this;
}

Optimizer Polyhedron::optimizer(){
    return Optimizer(*this);
}

void Polyhedron::set_f2v_sorted(bool b){
    f2v_sorted = b;
}

MatrixXd Polyhedron::render(int res, int thickness, int aa_factor){
    /*
            Renderer(
            Polyhedron &poly_,
            int width_,
            int height_,
            double rx_,
            double ry_,
            int aa_level_,
            int num_threads_,
            int thickness_);
    */
    Renderer renderer(*this, res, res, 1.2, 1.2, aa_factor, 1, thickness);
    return renderer.render();
}
