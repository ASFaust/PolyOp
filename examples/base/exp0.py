import polyop as pp
import matplotlib.pyplot as plt
import numpy as np

p = pp.icosahedron().ambo().dual().trunc().propeller().gyro()

for i in range(1000):
    p.move_to_face_centers()
    p.normalize()
p.compute_automorphisms()


p.align()
img_array = p.render(800, 4, 2).copy() # render the polyhedron to a numpy array
img_array /= img_array.max() # normalize the array to the range [0,1]
img_array = 1.0 - img_array # invert the colors so that the lines are white and the background is black
plt.imshow(img_array, cmap="gray") # show the image
plt.axis("off") # hide the axes
plt.savefig("results/exp0/image.png", bbox_inches="tight", pad_inches=0) # save the image

#save the result
#p.save_obj("results/exp0/cube.obj",overwrite=True)
