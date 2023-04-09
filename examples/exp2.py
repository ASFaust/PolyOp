import polyop

#this creates a planar drawing of a truncated cube
#triangulating the first face. the new vertex is added to the end of the vertex list
a = polyop.cube().trunc().triangulate([0],mode="center")

#spectral realization using standard arguments
a.realize_spectral()

#but iterating is necessary for some reason
#it seems to converge better with face_mean and also normalize is needed 
#because the gram-schmidt orthogonalization aligns the symmetry axes with the 
#coordinate system axes due to the non-equal treatment of the coordinate system axes
opt = a.optimizer()
for i in range(1000):
    a.move_to_plane()
    a.normalize()
    opt.step()

#getting a planar drawing centered on the last vertex, which can be obtained via get_vertices()[-1]
#and setting dimension 0 to zero
#which dimension needs to be set to 0 depends on the polyhedron
a.get_planar_drawing(a.get_vertices()[-1], 0)

#copying over the result to a fresh truncated cube because we want to get rid of the center vertex
#it was only used to get the planar drawing centered on a face and also to create a symmetry axis which
#the gram schmidt orthogonalization aligns to the coordinate system.
#if the symmetry axis wasn't aligned with a coordinate system axis, we wouldn't be able to just flatten that dimension
#and still get a symmetric planar drawing
b = polyop.cube().trunc()
b.set_pos(a.get_pos())
b.save_obj("results/planar_tricube.obj")

