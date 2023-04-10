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
    if (d <= 0 || d >= 1){
	std::invalid_argument("decay must be greater than 0 and less than 1");
    }
    decay = d;
    strength = s;
}

void Optimizer::step(){
    auto force = polyhedron.pos_ - last_pos;
    acc = decay * acc + strength * force;
    polyhedron.pos_ = last_pos + acc;
    last_pos = polyhedron.pos_;
    //ADAM variant:
}
