from sasa_calculator import parse_pdb_file, calculate_sasa_for_residue

# Load atoms from a PDB file
atoms = parse_pdb_file("small_h.pdb")

# Calculate SASA for a specific residue
residue_name = "LYS"  # Replace with the desired residue name
residue_seq_num = 48  # Replace with the desired residue sequence number
sasa_residue = calculate_sasa_for_residue(atoms, residue_name, residue_seq_num)
print(f"SASA for residue {residue_name} {residue_seq_num}: {sasa_residue} Å²")
