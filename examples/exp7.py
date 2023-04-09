#examples that do not require matplotlib

import polyop as pp
import numpy as np
import imageio
import cv2

thickness = 8

p = pp.dodecahedron()

resolution = 512
aa_factor = 4

p.align(steps=100)
#here, the conway operators are applied in the order they are listed

def sequence(i):
    p2 = p.scrape_off(i/100.0)
    p2.normalize()
    img_array = p2.render(resolution, thickness, aa_factor).copy()
    img_array /= img_array.max()
    img_array = ((1.0 - img_array) * 255).astype(np.uint8)
    cv2.imshow("image", img_array)
    cv2.waitKey(1)
    return img_array

frames = []

for i in range(20):
    frames.append(sequence(i))

p = p.scrape_off(20/100.0)

for i in range(101):
    frames.append(sequence(i))

p = p.dual()

#for i in range(33):
#    sequence(i)
#
#p = p.scrape_off(0.33)

p.normalize()

p2 = p.propeller()
img_array2 = p2.render(resolution, thickness, aa_factor).copy()
img_array2 /= img_array2.max()
img_array2 = ((1.0 - img_array2) * 255).astype(np.uint8)
img_array = p.render(resolution, thickness, aa_factor).copy()
img_array /= img_array.max()
img_array = ((1.0 - img_array) * 255).astype(np.uint8)

#make a slow fade in of the propeller

for i in range(75):
    img_array3 = img_array * (1.0 - i/74.0) + img_array2 * (i/74.0)
    frames.append(img_array3)
    cv2.imshow("image", img_array3.astype(np.uint8))
    cv2.waitKey(1)

p = p2

opt = p.optimizer()

#decay, strength
opt.set_momentum(0.9,0.01)

for i in range(75):
    p.move_to_face_centers()
    p.normalize()
    opt.step()

    img_array = p.render(resolution, thickness, aa_factor).copy()
    img_array /= img_array.max()
    img_array = ((1.0 - img_array) * 255).astype(np.uint8)
    cv2.imshow("image", img_array)
    cv2.waitKey(1)
    frames.append(img_array)

#opt.set_momentum(0.99, 0.01)

for i in range(101):
    frames.append(sequence(i))

p = p.dual()
opt = p.optimizer()

#decay, strength
opt.set_momentum(0.9,0.1)

for i in range(75):
    p.move_to_equal_edge_length()
    p.normalize()
    opt.step()

    img_array = p.render(resolution, thickness, aa_factor).copy()
    img_array /= img_array.max()
    img_array = ((1.0 - img_array) * 255).astype(np.uint8)
    cv2.imshow("image", img_array)
    cv2.waitKey(1)
    frames.append(img_array)

delta = 0
#now rotate it a bit
for i in range(300):
    if(i < 100):
        delta += 0.01
    if(i > 200):
        delta -= 0.01
    p.rotate(delta * 0.01, 1)
    p.rotate(delta * delta * 0.02, 2)
    p.normalize()
    #opt.step()

    img_array = p.render(resolution, thickness, aa_factor).copy()
    img_array /= img_array.max()
    img_array = ((1.0 - img_array) * 255).astype(np.uint8)
    cv2.imshow("image", img_array)
    cv2.waitKey(1)
    frames.append(img_array)

print(f"len(frames) = {len(frames)}")

imageio.mimsave('results/start.gif', frames, fps = 24, palettesize=256)
