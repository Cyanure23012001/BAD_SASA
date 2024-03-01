# SASA Calculator Library and CLI

The `sasa_calculator`(.py) is a Python library designed for calculating the Solvent Accessible Surface Area (SASA) of proteins from PDB (Protein Data Bank) files. It provides a set of functions to parse PDB files, calculate SASA for individual residues or entire proteins, and visualize atom points using 3D plotting. The library also includes a command-line interface (`sasa_cli`)(.py) for easy use directly from the terminal.

## Library Functions

### `parse_pdb_file(pdb_file)`
- **Description**: Parses a PDB file and extracts atomic and residue information.
- **Parameters**:
  - `pdb_file` (str): Path to the PDB file.
- **Returns**: A NumPy array containing atomic data.

### `calculate_sasa(atoms, num_points=100)`
- **Description**: Calculates the Solvent Accessible Surface Area (SASA) for a set of atoms.
- **Parameters**:
  - `atoms` (np.ndarray): Array of atoms.
  - `num_points` (int, optional): Number of points to use in the calculation.
- **Returns**: The calculated SASA (float).

### `calculate_sasa_for_residue(atoms, residue_name, residue_seq_num, num_points=100)`
- **Description**: Calculates the SASA for a specific residue.
- **Parameters**:
  - `atoms` (np.ndarray): Array of atoms.
  - `residue_name` (str): Name of the residue.
  - `residue_seq_num` (int): Sequence number of the residue.
  - `num_points` (int, optional): Number of points to use in the calculation.
- **Returns**: SASA of the specified residue (float).

### `calculate_and_write_sasa_for_all_residues(input_pdb, output_pdb, atoms, num_points=100)`
- **Description**: Calculates SASA for all residues in a protein and writes the results to a new PDB file.
- **Parameters**:
  - `input_pdb` (str): Path to the input PDB file.
  - `output_pdb` (str): Path to the output PDB file with SASA values.
  - `atoms` (np.ndarray): Array of atoms from the input PDB file.
  - `num_points` (int, optional): Number of points to use in the calculation.

### `plot_atom_points_plotly(atoms, num_points=100)`
- **Description**: Creates a 3D plot of atoms using Plotly.
- **Parameters**:
  - `atoms` (np.ndarray): Array of atoms.
  - `num_points` (int, optional): Number of points to use for each atom's surface.

## Command-Line Interface (CLI) - `sasa_cli`

### Usage

```bash
python sasa_cli.py <pdb_file> [options]
```
### Options
- `--residue <Name,SeqNum>`: Calculate SASA for a specific residue.
- `--all`: Calculate and write SASA for all residues (requires `--output`).
- `--plot`: Plot atom points using Plotly.
- `--output <file>`: Specify the output PDB file for SASA values (required with `--all`).

### Examples

- **Default (Calculate SASA for Entire Protein)**:
  ```bash
  python sasa_cli.py your_pdb_file.pdb
  ```
- **Calculate SASA for a Specific Residue**:
  ```bash
  python sasa_cli.py your_pdb_file.pdb --residue LYS,48
  ```
- **Calculate SASA for all Residues and save it**:
  ```bash
  python sasa_cli.py your_pdb_file.pdb --all --output output_with_sasa.pdb
  ```
- **Plot Atom Points**: 
    ```bash
    python python sasa_cli.py your_pdb_file.pdb --plot
  ```
