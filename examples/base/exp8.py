import polyop as pp

icosahedron = pp.icosahedron()

#align the icosahedron with the x-axis
icosahedron.align(steps=1000)
icosahedron.normalize()

icosahedron.save_obj("results/aligned_icos.obj", overwrite = True)

prism = pp.prism(6) #create a prism with 6 sides
prism.align(steps=1000)
prism.normalize()
prism.save_obj("results/aligned_prism.obj", overwrite = True)

