import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay  # For tetrahedral mesh generation

# Physical Constants
mu_0 = 4 * np.pi * 1e-7  # Permeability of free space


def biot_savart_3d_fem(nodes, elements, J, evaluation_points, quadrature_points=3):
    """
    Calculates the magnetic field from a 3D cylindrical current distribution using the Biot-Savart Law.

    Args:
        nodes (np.array): Array of node coordinates (N x 3 for 3D).
        elements (np.array): Array of element connectivity (M x num_nodes_per_element). Tetrahedral or Hexahedral mesh
        J (function): Current density function J(x, y, z) vector with (3) components.
        evaluation_points (np.array): Array of points where to evaluate the magnetic field (P x 3).
        quadrature_points (int): Number of quadrature points for integration (default: 3).

    Returns:
        B (np.array): Magnetic field at the evaluation points (P x 3).
    """
    B = np.zeros(evaluation_points.shape)
    gauss_points, gauss_weights = np.polynomial.legendre.leggauss(
        quadrature_points
    )  # Gauss-Legendre quadrature
    for element in elements:
        # loop through elements, using the J function and its properties to calculate results

        element_nodes = nodes[element]  # Nodes for the element
        # Numerical Integration - Add support to other types of geometry - This example supports simple tetrahedron mesh

        for i, point in enumerate(evaluation_points):
            # Evaluate center of elements
            center_point = np.mean(element_nodes, axis=0)
            # Calculate dl, calculate r, calculate the Biot Savart function, implement to Gauss integration and apply
            # Calculate r
            r = point - center_point
            # Apply the Biot Savart law
            J_value = J(
                center_point[0], center_point[1], center_point[2]
            )  # Current Density
            dB = (mu_0 / (4 * np.pi)) * (
                np.cross(J_value, r) / (np.linalg.norm(r) ** 3)
            )
            B[i] += dB

    return B


def cylindrical_current(x, y, z):
    """
    Example current density function for a cylindrical conductor.
    """
    radius = 0.1  # Radius of the cylinder
    current_density = 100  # Set J value
    if np.sqrt(x**2 + y**2) <= radius:
        return np.array([0, 0, current_density])  # Current in the z direction
    else:
        return np.array([0, 0, 0])  # Outside the wire


def create_cylinder_mesh(radius, height, n_radial, n_circumferential, n_axial):
    """
    Creates a structured mesh of a cylinder using NumPy and then converts it to tetrahedral elements using Delaunay.

    Args:
        radius (float): Radius of the cylinder.
        height (float): Height of the cylinder.
        n_radial (int): Number of elements in the radial direction.
        n_circumferential (int): Number of elements in the circumferential direction.
        n_axial (int): Number of elements in the axial direction.

    Returns:
        nodes (np.array): Array of node coordinates (N x 3).
        elements (np.array): Array of element connectivity (M x 4) for tetrahedra.
    """

    # Create Structured Mesh
    r = np.linspace(0, radius, n_radial + 1)  # Radius Vector
    theta = np.linspace(0, 2 * np.pi, n_circumferential + 1)  # Circumferential points
    z = np.linspace(0, height, n_axial + 1)  # Z points
    R, Theta, Z = np.meshgrid(r, theta, z)  # 3D Grid points

    # Cartesian Coordinates from Cylindrical
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)

    # Merge XYZ into nodes
    nodes = np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T

    # Delaunay Tessellation
    tri = Delaunay(nodes[:, :2])  # XY
    simplices = tri.simplices

    # To account for 3D, project Z to triangles
    elements = tri.simplices  # This needs to be calculated.

    # Remove duplicate elements (This can be problematic)
    elements = np.unique(elements, axis=0)

    return nodes, elements


# Example Usage - Run Calculations and Test it
radius = 0.1  # Radius of the wire
height = 1  # Length
nodes, elements = create_cylinder_mesh(
    radius, height, n_radial=4, n_circumferential=8, n_axial=10
)

# Print elements and number of points to analyse possible geometry
print("Number of Nodes", len(nodes))
print("Number of Elements", len(elements))

# Define the current
evaluation_points = np.array(
    [[0, 0, 0], [0.09, 0.09, 0.09], [0.001, 0.01, 0.2], [1, 1, 1]]
)
J = cylindrical_current
# B = biot_savart_3d_fem(nodes, elements, J, evaluation_points) #Run Biot Savart

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

# Create 3d quiver - It will not work directly. It needs to account for 3D. It needs to calculate and display something.
# X = nodes[:,0]
# Y = nodes[:,1]
# Z = nodes[:,2]
# U = B[:,0]
# V = B[:,1]
# W = B[:,2]

# ax.quiver(X, Y, Z, U, V, W, length=0.05, normalize=True)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.savefig("Bios3D.png")
