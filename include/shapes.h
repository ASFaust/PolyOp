#ifndef SHAPES_H
#define SHAPES_H

#include "Polyhedron.h"

inline Polyhedron cone(int sides){
    if(sides < 3){
        throw std::invalid_argument("sides must be more than 2!");
    }
    vector<vector<int>> f2v;
    vector<int> bottom_face;
    for(int i = 0; i < sides; i++){
        f2v.push_back({i,(i-1 + sides) % sides,sides});
        bottom_face.push_back(i);
    }
    f2v.push_back(bottom_face);
    auto ret = Polyhedron(f2v);
    ret.realize_spectral(ret.get_matrix("face_mean"),{-2,-3,-4});
    //ret.hang_axis(2); //noo dont hang axis haha
    ret.normalize();
    //cones need to be realized with the hanging vertex method
    return ret;
}

inline Polyhedron antiprism(int sides){
    if(sides < 3){
        throw std::invalid_argument("sides must be more than 2!");
    }
    vector<vector<int>> f2v;
    vector<int> f1,f2;
    for(int i = 0; i < sides; i++){
        int j = (i+1) % sides;
        f2v.push_back({i,j+sides,j});
        f2v.push_back({i,j+sides,i+sides});
        f1.push_back(i);
        f2.push_back(i+sides);
    }
    f2v.push_back(f1);
    f2v.push_back(f2);
    auto ret = Polyhedron(f2v);
    auto tr = ret.triangulate({},"center");
    //wait isnt face mean the same as vertex mean for triangulated polyhedron?
    tr.realize_spectral(tr.get_matrix("face_mean"),{-2,-3,-4});
    auto pos = tr.get_pos(); //this is a copying behaviour
    ret.set_pos(pos);
    ret.normalize();
    //antiprisms need to be realized with the triangulated form. let's see.
    return ret;
}

inline Polyhedron prism(int sides){
    auto ret = cone(sides).trunc(sides);
    ret.realize_spectral(ret.get_matrix("face_mean"),{-2,-3,-4});
    ret.normalize();
    //prisms can be realized with the standard matrix
    return ret;
}

//all platonic solids can be realized with the standard matrix
inline Polyhedron hexahedron(){
    return prism(4);
}

inline Polyhedron octahedron(){
    auto ret = hexahedron().dual();
    ret.realize_spectral(ret.get_matrix("face_mean"),{-2,-3,-4});
    ret.normalize();
    return ret;
}

inline Polyhedron tetrahedron(){
    return cone(3);
}

inline Polyhedron icosahedron(){
    auto ret = antiprism(5).kis(5);
    ret.realize_spectral(ret.get_matrix("face_mean"),{-2,-3,-4});
    ret.normalize();
    return ret;
}

inline Polyhedron dodecahedron(){
    auto ret = icosahedron().dual();
    ret.realize_spectral(ret.get_matrix("face_mean"),{-2,-3,-4});
    ret.normalize();
    return ret;
}

#endif
