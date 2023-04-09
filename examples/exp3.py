import polyop as pp

a = pp.tetrahedron().gyro().gyro()

a.check() #performs a series of tests on the integrity of the polyhedron

for i in range(1000):
    a.move_to_face_centers()
    a.normalize()

for i in range(1000):
    a.move_to_plane()
    a.normalize()

for i in range(100):
    a.move_to_canonical()

print(a)

a.save_obj("results/gyrated_tetrahedron.obj")
