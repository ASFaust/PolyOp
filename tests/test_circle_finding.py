import numpy as np
import cv2 as cv
#this script is a test for the circle finding algorithm
#for n given 3d points, it finds the circle that best fits them
def run():
    n_points = 4
    points = (np.random.randn(n_points, 2))
    center = np.random.randn(2)
    radius = np.random.rand() * 2 + 2
    #put all the points on the circle
    #for i in range(n_points):
    #    points[i] = center + radius * (points[i] - center) / np.linalg.norm(points[i] - center)

    #find the circle
    #the idea is to use a force-based algorithm: we dont estimate the radius, just the center.
    #the force is positive for points inside the circle and negative for points outside the circle

    #the force is proportional to the distance of the point to the circle

    center = np.mean(points, axis=0) #initial guess for # the center
    radius = np.mean(np.linalg.norm(points - center, axis=1))

    canvas_res = 500
    canvas_size = [-4, 4]

    last_forces = np.zeros(2)

    for i in range(100):
        forces = np.zeros(2)
        for point in points:
            distance = np.linalg.norm(point - center)
            forces += (point - center) * (distance - radius) / distance

        #we add a dynamic damping to the forces to avoid oscillations. the damping is proportional to the scalar product of the forces
        damp_factor = np.dot(forces, last_forces) / (np.linalg.norm(forces) * np.linalg.norm(last_forces) + 1e-6) > 0

        last_forces *= damp_factor
        last_forces += forces / n_points
        center += last_forces
        print(f"{i:03d} magnitude of forces: {np.linalg.norm(last_forces)}")
        #recompute the radius
        radius = np.mean(np.linalg.norm(points - center, axis=1))

        #draw the points and the circle on a canvas, using opencv
        canvas = np.ones((canvas_res, canvas_res, 3)) * 255
        for point in points:
            x = int((point[0] - canvas_size[0]) / (canvas_size[1] - canvas_size[0]) * canvas_res)
            y = int((point[1] - canvas_size[0]) / (canvas_size[1] - canvas_size[0]) * canvas_res)
            cv.circle(canvas, (x, y), 5, (0, 0, 255), -1)
        x = int((center[0] - canvas_size[0]) / (canvas_size[1] - canvas_size[0]) * canvas_res)
        y = int((center[1] - canvas_size[0]) / (canvas_size[1] - canvas_size[0]) * canvas_res)
        r = int(radius / (canvas_size[1] - canvas_size[0]) * canvas_res)
        cv.circle(canvas, (x, y), r, (0, 255, 0), 2)
        cv.imshow("canvas", canvas.astype(np.uint8))

        cv.waitKey(10)



while True:
    run()
