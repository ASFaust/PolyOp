#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
//for throwing invalid argument:
#include <stdexcept>

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[";
    bool first = true;
    for (const auto& elem : v) {
        if (!first) {
            os << ", ";
        }
        os << elem;
        first = false;
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>>& v) {
    os << "[";
    bool first = true;
    for (const auto& inner_vec : v) {
        if (!first) {
            os << "\n\t";
        }
        os << inner_vec;
        first = false;
    }
    os << "]";
    return os;
}


// Axis-aligned rotation matrix function
inline Eigen::Matrix3d rotationMatrix(int axis, double eps) {
    Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();
    double s = std::sin(eps);
    double c = std::cos(eps);

    switch (axis) {
        case 0: // Rotate around x-axis
            rot << 1, 0, 0,
                   0, c, -s,
                   0, s, c;
            break;
        case 1: // Rotate around y-axis
            rot << c, 0, s,
                   0, 1, 0,
                   -s, 0, c;
            break;
        case 2: // Rotate around z-axis
            rot << c, -s, 0,
                   s, c, 0,
                   0, 0, 1;
            break;
        default:
            throw std::invalid_argument("Invalid axis value. Must be 0, 1, or 2.");
            break;
    }

    return rot;
}

#endif // UTILS_H
