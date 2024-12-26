#examples that do not require matplotlib

import polyop as pp
import numpy as np
import cv2

p = pp.tetrahedron()
#here, the conway operators are applied in the order they are listed
p = p.ambo().gyro().trunc(degree = 3).dual().ambo().trunc().propeller()
#truncate is only applied to vertices with degree 3

#then the positions are optimized
#first the vertices are moved towards the average of the centers of the faces they are part of
#this has been shown to produce nice-looking results
for i in range(200):
    p.move_to_face_centers()
    p.normalize()

#but it produces faces that are not planar
#so the vertices are then moved to the planes of the faces they are part of:
for i in range(200):
    p.move_to_plane()
    p.normalize()

#the normalization command moves the polyhedron to the origin
#and scales it so that the average distance between the vertices and the origin is 1
#it also uses gram-schmidt orthogonalization to eliminate any skewness
#it is necessary to call this after every step that moves the vertices because the
#vertices are scaled in each move_to command

#align the symmetry axes with the coordinate system axes

#save the result
p.save_obj("results/exp1/output.obj",overwrite=True)

p.align(steps = 100)
# create a numpy array

# show
while True:
    p.rotate(0.02,2)
    p.rotate(0.01,1)
    p.move_to_face_centers()
    p.normalize()
    img_array = p.render(800, 4, 2).copy()
    img_array /= img_array.max()
    img_array = ((1.0 - img_array) * 255).astype(np.uint8)
    cv2.imshow("image", img_array)
    cv2.waitKey(1)
