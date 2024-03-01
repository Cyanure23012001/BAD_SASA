import numpy as np
import numba
# Constants
ELEMENT_SYMBOLS = np.array(["H", "HE", "C", "N", "O", "F", "NA", "MG", "P", "S", 
                            "CL", "K", "CA", "NI", "CU", "ZN", "SE", "BR", "CD", 
                            "I", "HG"])
ELEMENT_RADII = np.array([1.200, 1.400, 1.700, 1.550, 1.520, 1.470, 2.270, 1.730, 
                          1.800, 1.800, 1.750, 2.750, 2.310, 1.630, 1.400, 1.390, 
                          1.900, 1.850, 1.580, 1.980, 1.550])
SOLVENT_RADIUS = 1.4

def parse_pdb_file(pdb_file):
    """
    Parses a PDB file to extract atomic coordinates and element symbols.

    Parameters:
    pdb_file (str): Path to the PDB file.

    Returns:
    tuple: A tuple containing two numpy arrays, one for numeric data 
           (atomic coordinates) and one for string data (element symbols).
    """
    numeric_data = []
    string_data = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x_coordinate = float(line[30:38].strip())
                y_coordinate = float(line[38:46].strip())
                z_coordinate = float(line[46:54].strip())
                numeric_data.append([x_coordinate, y_coordinate, z_coordinate])
                
                atom_name = line[12:16].strip()
                element_symbol = line[76:78].strip()
                string_data.append([atom_name, element_symbol])

    return np.array(numeric_data), np.array(string_data)

@numba.jit(nopython=True)
def get_radius_vectorized(elements, symbols, radii):
    """
    Finds the radii for a given array of element symbols.

    Parameters:
    elements (numpy array): Array of element symbols.
    symbols (numpy array): Array of available element symbols.
    radii (numpy array): Array of corresponding radii for the elements.

    Returns:
    numpy array: Array of radii for the given element symbols.
    """
    index = np.searchsorted(symbols, elements)
    valid_indices = index < len(radii)
    radii_found = np.empty(len(elements), dtype=np.float64)  # Specify dtype directly
    radii_found.fill(np.nan)  
    for i in range(len(elements)):
        radii_found[i] = radii[index[i]]
        
    return radii_found

@numba.jit(nopython=True)
def calc_points_vectorized(center, radius, n):
    """
    Calculates points on the surface of a sphere given a center and radius.

    Parameters:
    center (list): The center of the sphere (x, y, z coordinates).
    radius (float): The radius of the sphere.
    n (int): Number of points to calculate on the sphere's surface.

    Returns:
    numpy array: An array of points on the surface of the sphere.
    """
    points = np.zeros((n, 3))
    N = n
    theta = np.arccos(-1 + 2 * np.linspace(0, N - 1, N) / float(N - 1))
    theta[0] += 1e-6  # Small offset to avoid division by zero
    theta[-1] -= 1e-6

    phi_increment = 3.6 / np.sqrt(N * (1 - np.square(np.cos(theta))))
    phi = phi = np.cumsum(np.concatenate((np.array([0]), phi_increment[:-1]))) % (2 * np.pi)


    x = center[0] + radius * np.sin(phi) * np.sin(theta)
    y = center[1] + radius * np.cos(phi) * np.sin(theta)
    z = center[2] - radius * np.cos(theta)

    return np.column_stack((x, y, z))

@numba.jit(nopython=True)
def is_accessible_vectorized(point, current_atom_index, other_atoms, atomic_radii):
    """
    Determines if a point on a sphere is accessible or not.

    Parameters:
    point (numpy array): The point to check.
    current_atom_index (int): Index of the current atom.
    other_atoms (numpy array): Array of coordinates of other atoms.
    atomic_radii (numpy array): Array of radii of the atoms.

    Returns:
    bool: True if the point is accessible, False otherwise.
    """
    mask = np.arange(other_atoms.shape[0]) != current_atom_index
    filtered_atoms = other_atoms[mask]
    filtered_radii = atomic_radii[mask]

    distances = np.sqrt(np.sum(np.square(filtered_atoms - point), axis=1))
    return np.all(distances >= filtered_radii + SOLVENT_RADIUS)

@numba.jit(nopython=True)
def calculate_sasa_vectorized(numeric_atoms, string_data, symbols, radii, num_points=100):
    sasa = 0.0
    atomic_radii = get_radius_vectorized(string_data[:, 1], symbols, radii)
    adjusted_radii = atomic_radii + SOLVENT_RADIUS
    area_per_point = 4 * np.pi * np.square(adjusted_radii) / num_points

    for i in range(numeric_atoms.shape[0]):
        surface_points = calc_points_vectorized(numeric_atoms[i], adjusted_radii[i], num_points)
        # Convert the boolean list comprehension to a NumPy array before summing
        accessible_points = np.array([is_accessible_vectorized(point, i, numeric_atoms, atomic_radii) for point in surface_points]).sum()
        sasa += accessible_points * area_per_point[i]

    return sasa

numeric_atoms, string_data = parse_pdb_file("small_h.pdb")
sasa = calculate_sasa_vectorized(numeric_atoms, string_data, ELEMENT_SYMBOLS, ELEMENT_RADII)
print(f"SASA: {sasa} Å²")
