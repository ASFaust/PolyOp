import polyop as pp

p = pp.cube()

p.align()

#save the result
p.save_obj("results/exp0/cube.obj",overwrite=True)
