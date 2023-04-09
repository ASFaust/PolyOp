import polyop as pp
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.tri as mtri
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection

a = pp.prism(6).trunc().gyro().ambo().gyro().dual()
for i in range(1000):
    a.face_mean()
#for i in range(100):
#    x = a.flatten_faces()

a.save_obj("test.obj")
#triangulating gives a bigger chance of success
