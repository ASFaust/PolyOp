#examples that do not require matplotlib

import polyop as pp
import numpy as np

dodecahedron = pp.dodecahedron()
icosahedron = pp.icosahedron()
octahedron = pp.octahedron()
tetrahedron = pp.tetrahedron()
cube = pp.cube()
cones = []
for i in range(4, 8):
    cones.append(pp.cone(i))
prisms = []
for i in range(5, 8):
    prisms.append(pp.prism(i))
antiprisms = []
for i in range(3, 8):
    antiprisms.append(pp.antiprism(i))

dodecahedron.save_obj("results/dodecahedron.obj")
icosahedron.save_obj("results/icosahedron.obj")
octahedron.save_obj("results/octahedron.obj")
tetrahedron.save_obj("results/tetrahedron.obj")
cube.triangulate(mode="center").save_obj("results/cube.obj")

for i in range(len(cones)):
    cones[i].triangulate(mode="center").save_obj("results/cone{}.obj".format(i + 4))

for i in range(len(prisms)):
    prisms[i].triangulate(mode="center").save_obj("results/prism{}.obj".format(i + 5))

for i in range(len(antiprisms)):
    antiprisms[i].triangulate(mode="center").save_obj("results/antiprism{}.obj".format(i + 3))


