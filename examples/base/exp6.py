#examples that do not require matplotlib

import polyop as pp
import numpy as np

#Polyhedron Polyhedron::detach_faces(double scaling, double moving, double thickness, bool pyramidal){
p = pp.cube()
p.align()

p = p.detach_faces(0.1, 0.1, 0.1, False) #detaching faces is broken right now.

#to fix:
#1. normals are inconsistently pointing inwards or outwards, leading to skewed results when using
#them as thickness and moving vectors
#problem is that orienting the normals on a previously detached polyhedron will be impossible,
#because how should we know which normals are pointing inwards and which are pointing outwards
#when there are multiple meshes? would need to determine for each mesh, but the meshes are only independent
#implicitly. so we would need to check for "islands" first and treat each "island" separately. 

#2. one position is not set correctly, i don't know which one and i don't know why.
#i dont really understand what is happening, and it probably would be best to rewrite the whole thing



p.save_obj("results/exp6/result.obj",overwrite=True)
p.check()
