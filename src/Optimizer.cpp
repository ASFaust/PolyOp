#include "Optimizer.h"
#include "Polyhedron.h"
#include "Edge.h"
#include "includes.h"

Optimizer::Optimizer(Polyhedron& p_) : polyhedron(p_) {
    last_pos = polyhedron.get_pos();
    acc = polyhedron.get_pos();
    acc.setZero();
    decay = 0.9;
    strength = 0.1;
}

//set decay and strength
void Optimizer::set_momentum(double d, double s){
    decay = d;
    if(decay >= 1){
        cerr << "Decay must be < 1. Setting to 0.99999" << endl;
        decay = 0.99999;
    }else if(decay <= 0){
        cerr << "Decay must be positive. Setting to 0." << endl;
        decay = 0;
    }
    strength = s;
    if(strength <= 0){
        cerr << "Strength must be greater than 0. Setting to 0.001." << endl;
        strength = 0.001;
    }
    if(strength > 1){
        cerr << "Strength must be <= 1. Setting to 1." << endl;
        strength = 1;
    }
}

void Optimizer::step(){
    auto force = polyhedron.pos_ - last_pos;
    acc = decay * acc + strength * force;
    polyhedron.pos_ = last_pos + acc;
    last_pos = polyhedron.pos_;
    //ADAM variant:
}
