from polyop import *
import cv2
import numpy as np


def test():
    a = cone(7).propeller().dual().gyro().dual()

    a.align(steps=1000)

    b = DetachedFacesOpt(a)

    for i in range(1000):
        img = a.render(512,8,4).copy()
        img /= img.max()
        img = ((1.0 - img) * 255).astype(np.uint8)

        print(".",end="",flush=True)
        cv2.imshow("image", img)
        cv2.waitKey(1)

        b.step(0.1)

test()