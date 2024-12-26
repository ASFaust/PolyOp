#include "Polyhedron.h"
#include "Optimizer.h"
#include "CirclePackingRepr.h"

#include "shapes.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

typedef Matrix<double,Dynamic,Dynamic,RowMajor> Mx;
vector<int> non;

PYBIND11_MODULE(polyop, m) {
    py::class_<Polyhedron>(m, "Polyhedron")
        .def("trunc",&Polyhedron::trunc_selected,py::arg("vertices") = non)
        .def("trunc",&Polyhedron::trunc,py::arg("degree") = -1)
        .def("ambo",&Polyhedron::ambo)
        .def("dual",&Polyhedron::dual)
        .def("kis",&Polyhedron::kis_selected,py::arg("faces") = non)
        .def("kis",&Polyhedron::kis,py::arg("degree") = -1)
        .def("join",&Polyhedron::join)
        .def("gyro",&Polyhedron::gyro)
        .def("propeller",&Polyhedron::propeller)
        .def("triangulate",&Polyhedron::triangulate,py::arg("faces") = non,py::arg("mode") = "center")
        .def("scrape_off",&Polyhedron::scrape_off,py::arg("t"))
        .def("detach_faces",&Polyhedron::detach_faces)

        .def("save_obj",&Polyhedron::save_obj,py::arg("filename"),py::arg("overwrite") = false,py::arg("create_folders") = true)
        .def("check",&Polyhedron::check)
	.def("__str__",&Polyhedron::string_repr)

        .def("optimizer",&Polyhedron::optimizer)
        .def("render",&Polyhedron::render)

        .def("to_origin",&Polyhedron::to_origin)
        .def("scale",&Polyhedron::scale,py::arg("factor") = 1.0)
        .def("unskew",&Polyhedron::unskew)
        .def("normalize",&Polyhedron::normalize) //to_origin,unskew,scale
        .def("align",&Polyhedron::align,py::arg("steps") = 100)
        .def("rotate",&Polyhedron::rotate,py::arg("angle"),py::arg("axis"))

        .def("move_to_plane",&Polyhedron::move_to_plane)
        .def("move_to_face_centers",&Polyhedron::move_to_face_centers)
        .def("move_to_neighbors",&Polyhedron::move_to_neighbors)
        .def("move_to_canonical",&Polyhedron::move_to_canonical)
        .def("move_to_equal_edge_length",&Polyhedron::move_to_equal_edge_length)

        .def("realize_spectral",&Polyhedron::realize_spectral,py::arg("matrix") = Mx::Zero(0,0),py::arg("indices") = non)
        .def("get_matrix",&Polyhedron::get_matrix)
        .def("hang_axis",&Polyhedron::hang_axis,py::arg("axis") = 2)
        .def("set_log",&Polyhedron::set_log)
        .def_readonly("v_count",&Polyhedron::v_count)
        .def("set_pos",&Polyhedron::set_pos);


        //.def("spherical",&Polyhedron::spherical)
        //.def("flatten_faces",&Polyhedron::flatten_faces)
        //.def("get_pos",&Polyhedron::get_pos)
        //.def_readonly("faces",&Polyhedron::f2v)
        //.def("flip_winding",&Polyhedron::flip_winding)
        //.def("face_mean",&Polyhedron::face_mean)
        //.def("vertex_mean",&Polyhedron::vertex_mean)
        //.def("get_vertices",&Polyhedron::get_vertices)
        //.def("get_edges",&Polyhedron::get_edges)

        //.def("get_planar_drawing",&Polyhedron::get_planar_drawing)
        //.def_readwrite("pos",&Polyhedron::_pos)
        //.def("scrape_off",&Polyhedron::scrape_off);
    py::class_<Optimizer>(m, "Optimizer")
        .def("set_momentum",&Optimizer::set_momentum)
        .def("step",&Optimizer::step);

    py::class_<CirclePackingRepr>(m, "CirclePackingRepr")
        .def(py::init<Polyhedron&, vector<int>>())
        .def("step",&CirclePackingRepr::step)
        .def("get_pos",&CirclePackingRepr::get_pos)
        .def("draw",&CirclePackingRepr::draw,py::arg("res") = 100,py::arg("thickness") = 0.01);

    m.def("prism",&prism,py::arg("sides"));
    m.def("antiprism",&antiprism,py::arg("sides"));
    m.def("cone",&cone,py::arg("sides"));
    m.def("cube",&hexahedron);
    m.def("hexahedron",&hexahedron);
    m.def("tetrahedron",&tetrahedron);
    m.def("octahedron",&octahedron);
    m.def("dodecahedron",&dodecahedron);
    m.def("icosahedron",&icosahedron);
};
