import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Number of points
n_points = 5

# Generate random 3D points
points = np.random.randn(n_points, 3)

# Define a random circle
center = np.random.randn(3)
radius = np.random.rand() * 2 + 2

# Place all points approximately on the circle surface for testing
for i in range(n_points):
    direction = np.random.randn(3)
    direction /= np.linalg.norm(direction)
    points[i] = center + radius * direction

# Initial guess for the center and radius
center = np.mean(points, axis=0)
radius = np.mean(np.linalg.norm(points - center, axis=1))

# Iterative force-based adjustment
last_forces = np.zeros(3)

for i in range(1000):
    forces = np.zeros(3)

    for point in points:
        distance = np.linalg.norm(point - center)
        forces += (point - center) * (distance - radius) / distance

    # Dynamic damping to avoid oscillations
    damp_factor = np.dot(forces, last_forces) / (np.linalg.norm(forces) * np.linalg.norm(last_forces) + 1e-6) > 0

    last_forces *= damp_factor
    last_forces += forces / n_points
    center += last_forces

    print(f"{i:03d} magnitude of forces: {np.linalg.norm(last_forces)}")

    # Recompute the radius
    radius = np.mean(np.linalg.norm(points - center, axis=1))

    # Stop iteration if forces are negligible
    if np.linalg.norm(last_forces) < 1e-6:
        break

# Visualization using Matplotlib
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Draw the points
ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='r', label='Points')

# Draw the estimated circle (as a sphere for visualization)
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
z = center[2] + radius * np.outer(np.ones_like(u), np.cos(v))
ax.plot_surface(x, y, z, color='b', alpha=0.3, label='Fitted Circle')

# Set labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.legend()
plt.show()
