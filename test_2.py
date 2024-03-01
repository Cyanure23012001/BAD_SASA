import math
import numpy as np
import pandas as pd

ATOMIC_RADII = {
    "H": 1.200,
    "HE": 1.400,
    "C": 1.700,
    "N": 1.550,
    "O": 1.520,
    "F": 1.470,
    "NA": 2.270,
    "MG": 1.730,
    "P": 1.800,
    "S": 1.800,
    "CL": 1.750,
    "K": 2.750,
    "CA": 2.310,
    "NI": 1.630,
    "CU": 1.400,
    "ZN": 1.390,
    "SE": 1.900,
    "BR": 1.850,
    "CD": 1.580,
    "I": 1.980,
    "HG": 1.550,
}
SOLVENT_RADIUS = 1.4

def parse_pdb_file(pdb_file):
    data = []
    residues = {}
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x_coordinate = float(line[30:38].strip())
                y_coordinate = float(line[38:46].strip())
                z_coordinate = float(line[46:54].strip())
                atom_name = line[12:16].strip()
                element_symbol = line[76:78].strip()
                residue_name = line[17:20].strip()
                residue_seq_num = int(line[22:26].strip())

                data.append([x_coordinate, y_coordinate, z_coordinate, atom_name, element_symbol, residue_name, residue_seq_num])

    return np.array(data, dtype=object)


def calc_points(center, radius, n):
    points = []
    N = n
    for k in range(1, N + 1):
        h = -1 + 2 * (k - 1) / float(N - 1)
        theta = math.acos(h)
        if k == 1 or k == N:
            phi = 0
        else:
            phi += 3.6 / math.sqrt(N * (1 - h * h))

        x = center[0] + radius * math.sin(phi) * math.sin(theta)
        y = center[1] + radius * math.cos(phi) * math.sin(theta)
        z = center[2] - radius * math.cos(theta)

        points.append([x, y, z])
        phi %= 2*math.pi

    return points

def is_accessible(point, current_atom, other_atoms):
    for atom in other_atoms:
        if np.array_equal(atom, current_atom):
            continue
        dx = point[0] - atom[0]
        dy = point[1] - atom[1]
        dz = point[2] - atom[2]
        distance = math.sqrt(dx**2 + dy**2 + dz**2)
        if distance < ATOMIC_RADII[atom[-3]] + SOLVENT_RADIUS:
            return False
    return True

def calculate_sasa(atoms, num_points=100):
    sasa = 0
    for atom in atoms:
        accessible_points = 0
        radius = ATOMIC_RADII[atom[-3]] + SOLVENT_RADIUS
        surface_points = calc_points(atom[:3], radius, num_points)

        for point in surface_points:
            if is_accessible(point, atom, np.delete(atoms, np.where(np.all(atoms == atom, axis=1)), axis=0)):
                accessible_points += 1

        area_per_point = 4 * math.pi * radius**2 / num_points
        sasa += accessible_points * area_per_point

    return sasa

def write_sasa_to_pdb(input_pdb, output_pdb, sasa_values):
    min_sasa = min(sasa_values.values())
    max_sasa = max(sasa_values.values())
    range_sasa = max_sasa - min_sasa
    normalized_sasa_values = {key: (value - min_sasa) / range_sasa if range_sasa else 0 for key, value in sasa_values.items()}

    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                residue_name = line[17:20].strip()
                residue_seq_num = int(line[22:26].strip())
                sasa = normalized_sasa_values.get((residue_name, residue_seq_num), 0)
                new_line = line[:60] + f"{sasa:6.2f}" + line[66:]
                outfile.write(new_line)
            else:
                outfile.write(line)

def calculate_sasa_for_residue(atoms, residue_name, residue_seq_num, num_points=100):
    filtered_atoms = [atom for atom in atoms if atom[5] == residue_name and atom[6] == residue_seq_num]
    
    if not filtered_atoms:
        print(f"No atoms found for residue {residue_name} {residue_seq_num}.")
        return 0

    return calculate_sasa(np.array(filtered_atoms), num_points)

def calculate_and_write_sasa_for_all_residues(input_pdb, output_pdb, atoms, num_points=100):
    # Extract unique residues
    unique_residues = set((atom[5], atom[6]) for atom in atoms)

    # Calculate SASA for each residue
    sasa_values = {}
    for residue_name, residue_seq_num in unique_residues:
        sasa = calculate_sasa_for_residue(atoms, residue_name, residue_seq_num, num_points)
        sasa_values[(residue_name, residue_seq_num)] = sasa

    # Write SASA values to the PDB file
    write_sasa_to_pdb(input_pdb, output_pdb, sasa_values)



def plot_atom_points_plotly(atoms, num_points=100):
    import plotly.graph_objects as go
    import pandas as pd
    # Define a color map for the elements
    element_colors = {
        'H': 'blue',
        'C': 'black',
        'N': 'red',
        'O': 'green',
        'P': 'orange',
        'S': 'yellow',
        # Add more elements and their colors as needed
    }

    all_points = []
    for atom in atoms:
        radius = ATOMIC_RADII[atom[-3]] + SOLVENT_RADIUS
        surface_points = calc_points(atom[:3], radius, num_points)

        for point in surface_points:
            all_points.append([point[0], point[1], point[2], element_colors.get(atom[-3], 'gray')])

    points_df = pd.DataFrame(all_points, columns=['x', 'y', 'z', 'color'])

    fig = go.Figure(data=[go.Scatter3d(
        x=points_df['x'],
        y=points_df['y'],
        z=points_df['z'],
        mode='markers',
        marker=dict(
            size=2,
            color=points_df['color'],  # Set color using the mapped values
            opacity=0.8
        )
    )])

    fig.update_layout(
        title='3D Atom Points Visualization',
        scene=dict(
            xaxis_title='X Axis',
            yaxis_title='Y Axis',
            zaxis_title='Z Axis'
        )
    )
    fig.show()

input_pdb_file = "small_h.pdb"  # Replace with your input PDB file name
output_pdb_file = "output_with_sasa.pdb"  # Output file name

# Parse the PDB file to get atoms array
atoms = parse_pdb_file(input_pdb_file)

# Calculate and write SASA for all residues
calculate_and_write_sasa_for_all_residues(input_pdb_file, output_pdb_file, atoms)



residue_name = "LYS"  # Example residue name
residue_seq_num = 48   # Example residue sequence number
sasa_residue = calculate_sasa_for_residue(atoms, residue_name, residue_seq_num)
print(f"SASA for residue {residue_name} {residue_seq_num}: {sasa_residue} Å²")

sasa = calculate_sasa(atoms)
print(f"SASA: {sasa} Å²")

output_pdb_file = "output_with_sasa.pdb"  # Output file name
