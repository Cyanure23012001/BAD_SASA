from sasa_calculator import parse_pdb_file, plot_atom_points_plotly

# Load atoms from a PDB file
atoms = parse_pdb_file("small_h.pdb")

# Plot atom points
plot_atom_points_plotly(atoms)