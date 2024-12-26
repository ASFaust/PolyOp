import polyop as pp

a = pp.cube() #.gyro()

b = a.gyro()
b.check() #performs a series of tests on the integrity of the polyhedron

for i in range(1000):
    b.move_to_face_centers()
    b.normalize()

for i in range(1000):
    b.move_to_plane()
    b.normalize()


for i in range(101):
    c = a.scrape_off(i/100.0)
    c.save_obj("results/exp4/anim_{:03d}.obj".format(i))

for i in range(101):
    c = b.scrape_off(i/100.0)
    c.save_obj("results/exp4/gyro_anim_{:03d}.obj".format(i))

#print(a)
#a.save_obj("results/gyrated_tetrahedron.obj")
