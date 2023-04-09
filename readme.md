# `polyop`: Conway Polyhedron Operators

<video width="512" height="512" controls>
  <source src="examples/results/start.webm" type="video/webm">
  Your browser does not support the video tag.
</video>

This is a Python library written in C++ for working with [Conway polyhedron operators](https://en.wikipedia.org/wiki/Conway_polyhedron_notation).

## Installation
To install the library, run the following commands in the root directory of this repository:
```
pip install numpy setuptools pybind11 imageio opencv-python
pip install .
```
## Features 

### Seed polyhedra
 * Platonic Solids (Tetrahedron, Cube, Octahedron, Dodecahedron, Icosahedron)
 * N-sided Prisms, Antiprisms and Pyramids

### Conway Polyhedron Operators
* Truncation
* Ambo
* Dual
* Gyration
* Propeller
* Compound operators of the above

### Shape optimization for the resulting polyhedra
#### Iterative refinement forces
Iterative refinement forces are a set of forces that can be applied to a polyhedron to alter its shape.
They have in common that their respective energy is not minimized by a single application of the force.
So the forces need to be applied iteratively, and the polyhedron is moved to minimize the total energy of the respective forces.
  * Move vertices to average of neighbor vertex positions
  * Move vertices to average of the centers of the faces they are part of
  * Move vertices to equalize edge length over the whole polyhedron
  * Move vertices to make edges tangent to unit sphere
  * Move vertices to lie in the total least squares planes of the faces they are part of

#### Spectral Realization

Spectral realization provides a non-iterative way of optimizing the shape of a polyhedron. 
It works by finding the eigenvectors of a matrix that is related to the polyhedron, and then assigning certain eigenvectors to the vertex positions.

##### Available Matrices:
via the `Polyhedron.get_matrix(string name)` function:
* Adjacency matrix
* Face-mean-matrix (works best)
* Vertex-mean-matrix
* (Laplacian, but it's broken rn)

Face-mean matrix is the default & equivalent to the minimum energy of moving the vertices to the average of their face centers.
The iterative force method is more robust and faster in most cases since it is sparsely computed, so it is especially faster in more complex polyhedra.

#### Other Shape Altering Functions
* Unskewing using Gram-Schmidt-Orthogonalization
* Scaling and rotating the polyhedron
* Automatic rotational alignment of the polyhedron to the coordinate axes

### Visualization
  * Small render engine (WIP) for rendering the polyhedron edges in Orthographic projection (See movie above)
  * Saving the polyhedron as an obj file


## Usage
See the [examples](/examples/) folder for examples of how to use the library. 

### Quickstart
The code below instantiates a cube, then truncates it, dualizes it,
and saves it as an obj file.
```python
from polyop import cube
p = cube()
p = p.trunc().dual()
p.save_obj('notacube.obj')
```

### Documentation


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.