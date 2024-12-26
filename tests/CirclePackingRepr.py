from polyop import *
import cv2
import numpy as np


def test():
    a = antiprism(5)

    a.align(steps=1000)

    img = a.render(512, 8, 4).copy()
    img /= img.max()
    img = ((1.0 - img) * 255).astype(np.uint8)
    cv2.imshow("image", img)
    cv2.waitKey(1000)

    print(a)

    b = CirclePackingRepr(a,list(range(100)))

    for i in range(1000):
        img = b.draw(512,0.001).copy()
        img /= img.max()
        img = ((1.0 - img) * 255).astype(np.uint8)
        #make the image into a color image by copying the grayscale image into the three color channels
        img = cv2.merge([img,img,img])
        #we have special values where img == 127, which we want to convert to red pixels

        positions = np.where(img[:,:,0] == 127)

        if(positions[0].shape[0] > 0):

            #then use cv2 to draw circles at the positions
            for i in range(len(positions[0])):
                cv2.circle(img, (positions[1][i], positions[0][i]), 5, (0, 0, 255), -1)


        cv2.imshow("image", img)
        cv2.waitKey(1)
        b.step(0.1)

test()