import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting

# Physical Constants
mu_0 = 4 * np.pi * 1e-7  # Permeability of free space


def biot_savart_fem(nodes, elements, current, evaluation_points):
    """
    Calculates the magnetic field using the Biot-Savart Law by dividing the source into a set of finite elements,
    integrating the Biot-Savart Law at each finite element, and performing superposition on the B field.

    Args:
        nodes (np.array): Array of node coordinates (N x 3 for 3D).
        elements (np.array): Array of element connectivity (M x 2). Each element is a line segment.
        current (float): Current flowing in the wire (Amperes).
        evaluation_points (np.array): Array of points where to evaluate the magnetic field (P x 3).

    Returns:
        B (np.array): Magnetic field at the evaluation points (P x 3).
    """

    B = np.zeros(evaluation_points.shape) # Initialize a B field with the correct shape

    # Loop through each of the line elements
    for element in elements:
        # Element Nodes
        node_1 = nodes[element[0]]
        node_2 = nodes[element[1]]

        # Vector that goes from node 1 to node 2
        dl = node_2 - node_1 #Direction and Magnitude

        # Element Length
        element_length = np.linalg.norm(dl)

        # Midpoint of the wire element (Good Approximation)
        midpoint = (node_1 + node_2)/2

        # Loop through all evaluation points for the B field
        for i, point in enumerate(evaluation_points):

            # r vector and magnitude
            r = point - midpoint
            r_mag = np.linalg.norm(r)

            #Biot-Savart Law
            dB = (mu_0 / (4 * np.pi)) * (current * np.cross(dl, r) / (r_mag**3)) #Superposition
            B[i] += dB

    return B


# ---------------------- Example Usage ----------------------
# Define the geometry of the wire (example: a straight wire along the z-axis)
domain_start = 0.0
domain_end = 1.0
n_segments = 10 # number of segments
n_nodes = n_segments + 1
nodes = np.zeros((n_nodes, 3))  # Three Columns for x, y, z
nodes[:, 2] = np.linspace(0.*(domain_start+domain_end), 1*(domain_start+domain_end), n_nodes) #Equally Spaced Points - 1 meter wire

# Element Connectivity (line segments) - 1 less element than the number of nodes.
elements = np.zeros((n_segments, 2), dtype=int)
for i in range(n_segments):
    elements[i, 0] = i
    elements[i, 1] = i + 1


# Define the current
current = 1.0  # Ampere

# Define the points where to evaluate the magnetic field
n_eval_points = 50
evaluation_points = np.zeros((n_eval_points, 3))
evaluation_points[:, 0] = np.linspace(-0.5, 0.5, n_eval_points) # Equally spaced x values
evaluation_points[:, 1] = 0.1  # y location of the Evaluation Points
evaluation_points[:, 2] = 0.5 # z location of the Evaluation Points

# Calculate the magnetic field
B = biot_savart_fem(nodes, elements, current, evaluation_points)

# ---------------------- Plotting ----------------------
# Plot the magnetic field (By component)
plt.plot(evaluation_points[:, 0], B[:, 1], marker='o', label='By') #Plot By
plt.xlabel("x (m)")
plt.ylabel("By (T)")
plt.title("Magnetic Field from Straight Wire")
plt.grid(True)
plt.legend()
plt.savefig("Bios.png")

#3D Quiver plot example with more points and with current on the x direction
#Setup
n_segments = 20 # number of segments
nodes = np.zeros((n_segments + 1, 3))  # Three Columns for x, y, z
nodes[:, 0] = np.linspace(domain_start, domain_end, n_segments + 1) #Equally Spaced Points - 1 meter wire
print(nodes)
elements = np.zeros((n_segments, 2), dtype=int)
for i in range(n_segments):
    elements[i, 0] = i
    elements[i, 1] = i + 1
current = 1 #Amps
# Create a grid of evaluation points for the 3D Plot.
num_points = 5
x = np.linspace(-0.5, 0.5, num_points)
y = np.linspace(-0.5, 0.5, num_points)
z = np.linspace(-0.5, 0.5, num_points)
X, Y, Z = np.meshgrid(x, y, z)

# Flatten the Grid
evaluation_points = np.vstack((X.flatten(), Y.flatten(), Z.flatten())).T

#Calculate the B field
B = biot_savart_fem(nodes, elements, current, evaluation_points)

#Plot - The results are not great because there is a very coarse grid, a more refined model can be implemented to increase the results.
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#Normalize B - Without normalization the quivers will be scaled improperly
norm = np.linalg.norm(B, axis=1, keepdims=True)
B_normalized = B / norm

# ax.plot(elements[:,0], elements[:,1], ls="-", color="r")
ax.quiver(evaluation_points[:, 0], evaluation_points[:, 1], evaluation_points[:, 2],
              B_normalized[:, 0], B_normalized[:, 1], B_normalized[:, 2], length=0.2)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.title('B Field from Straight Wire')
plt.savefig("Bios1.png")


## Good integral

# def biot_savart_fem(nodes, elements, current, evaluation_points, quadrature_points=3):
#     """
#     Biot-Savart Law with Gaussian Quadrature

#     """
#     # Gaussian Quadrature
#     gauss_points, gauss_weights = np.polynomial.legendre.leggauss(quadrature_points)
#     #...
#     #Loop though point

#             #Numerical integration of the element
#             dB = np.zeros(point.shape)
#             for j in range(quadrature_points):

#                 # Map the Gaussian point to the element
#                 xi = 0.5 * ((node_2 - node_1) * gauss_points[j]) + 0.5 * (node_1 + node_2) # point in the integration
#                 dl = node_2 - node_1
#                 r = point - xi
#                 r_mag = np.linalg.norm(r)

#                 dB += (mu_0 / (4 * np.pi)) * (current * np.cross(dl, r) / (r_mag**3)) * gauss_weights[j]


#             B[i] += dB #Add superposition
#     return B

# import numpy as np

# def create_straight_wire_geometry(start_point, end_point, n_segments):
#     """
#     Creates nodes and elements for a straight wire in 3D space.

#     Args:
#         start_point (np.array): 3D coordinates of the starting point of the wire.
#         end_point (np.array): 3D coordinates of the ending point of the wire.
#         n_segments (int): Number of line segments (elements) to divide the wire into.

#     Returns:
#         nodes (np.array): Array of node coordinates (N x 3).
#         elements (np.array): Array of element connectivity (M x 2).
#     """

#     nodes = np.zeros((n_segments + 1, 3))
#     elements = np.zeros((n_segments, 2), dtype=int)

#     #Calculate node positions
#     for i in range(n_segments + 1):
#         t = i / n_segments
#         nodes[i, :] = start_point + t * (end_point - start_point) # Linear interpolation

#     #Define element connectivity
#     for i in range(n_segments):
#         elements[i, 0] = i
#         elements[i, 1] = i + 1

#     return nodes, elements

# # Example Usage - Define a straight wire
# start_point = np.array([0, 0, 0]) #Starts at origin
# end_point = np.array([1, 0, 0]) #Ends at (1, 0, 0)
# n_segments = 10

# nodes, elements = create_straight_wire_geometry(start_point, end_point, n_segments)

# print("Nodes:\n", nodes)
# print("Elements:\n", elements)