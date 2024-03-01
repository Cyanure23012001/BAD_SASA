import argparse
from sasa_calculator import (parse_pdb_file, calculate_sasa, 
                             calculate_sasa_for_residue, calculate_and_write_sasa_for_all_residues, 
                             plot_atom_points_plotly)

def main():
    parser = argparse.ArgumentParser(description="SASA Calculator for PDB files.")
    parser.add_argument("pdb_file", help="Path to the PDB file.")
    parser.add_argument("--residue", help="Calculate SASA for a specific residue. Format: Name,SeqNum")
    parser.add_argument("--all", action="store_true", help="Calculate and write SASA for all residues.")
    parser.add_argument("--plot", action="store_true", help="Plot atom points using Plotly.")
    parser.add_argument("--output", help="Output PDB file path for SASA values. Required with --all")

    args = parser.parse_args()

    atoms = parse_pdb_file(args.pdb_file)

    if args.residue:
        residue_name, residue_seq_num = args.residue.split(',')
        sasa_residue = calculate_sasa_for_residue(atoms, residue_name, int(residue_seq_num))
        print(f"SASA for residue {residue_name} {residue_seq_num}: {sasa_residue} Å²")
    elif args.all:
        if not args.output:
            parser.error("--all requires --output.")
        calculate_and_write_sasa_for_all_residues(args.pdb_file, args.output, atoms)
    elif args.plot:
        plot_atom_points_plotly(atoms)
    else:
        # Default behavior: calculate SASA for the entire protein
        total_sasa = calculate_sasa(atoms)
        print(f"Total SASA for the protein: {total_sasa} Å²")

if __name__ == "__main__":
    main()
