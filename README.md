# COCαDA - Large Scale Protein Interatomic Contact Optimization by Cα Distance Matrices

## Description

COCαDA (Contact Optimization by alpha-Carbon Distance Analysis) optimizes the calculation of atomic interactions in proteins, by using a set of fine-tuned Cα distances between every pair of aminoacid residues.
The code includes a customized parser for both PDB and CIF files, containing functionalities for handling large files, filtering out specific residues and interactions, and calculating geometric properties such as centroid and normal vectors for aromatic residues.

Additionaly, as comparison for demonstrating the efficiency of our method, a Biopython-dependent script is also included (in the folder Biopython), containing the same restrictions and contact definitions.

The contact types available for calculation are:
  - Hydrophobic
  - Hydrogen Bond
  - Attractive
  - Repulsive
  - Disulfide Bond
  - Salt Bridge
  - Aromatic Stacking

## Features

- **PDB and CIF Parsing:** Supports parsing both PDB and CIF files to extract atomic and residue information.
- **Residue Filtering:** Ignores low-quality atoms, water molecules, and problematic/erroneous information.
- **Interaction Analysis:** Identifies contacts between specified types of atoms based on predefined conditions.
- **Aromatic Residue Analysis:** Computes centroids and normal vectors for aromatic residues, to determine aromatic stacking contacts.

## Installation

### Prerequisites

- For both methods:
  - Python==3.x
  - NumPy==2.0.1
  - psutil==6.0.0

- Exclusively for BioPython:
  - Biopython==1.84

### Setup

1. Clone the repository:
   ```sh
   git clone https://github.com/rplemos/Contacts.git
   ```

2. Navigate into the project directory:
   ```sh
   cd Contacts
   ```

3. Set up a virtual environment (optional but recommended):
    ```sh
    python -m venv venv
    source venv/bin/activate
    ```

4. Install dependencies:
    ```sh
    pip install -r requirements.txt
    ```
## Usage
### To run the Flexible Distances script:

1. Ensure you are in the project directory and the virtual environment is activated (if used).

2. Run the script with the path to your folder or file:
    ```sh
    python3  main.py <-f> path_to_files/ [-m] [-c] [-s] [-o] [-h]
    ```
**Parameters:**
 - <-f> <--files>: List of files in pdb/cif format (at least one required).
 - [-m] [--mode]: Select "SingleCore" or "MultiCore" mode.
 - [-c] [--ncores]: Number of cores to use (only needed on Multi mode). Default runs with all available cores.
 - [-s] [--selcore]: Select specific core to run (int value).
 - [-o] [--output]: Outputs the detailed results to files in ./outputs.
 - [-h] [--help]: Shows usage and instructions.


## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
Still not published.

## Contact
For any questions or issues, please contact:

Rafael Pereira Lemos - PhD Student in Bioinformatics @ Federal University of Minas Gerais

Email: rafaellemos42@gmail.com

GitHub: https://github.com/rplemos

## Contributions and Acknowledgements
 - Prof. Raquel Cardoso de Melo Minardi, UFMG;
 - Prof. Sabrina de Azevedo Silveira, UFV;
 - Dr. Diego César Batista Mariano, UFMG;
 - All the 'Laboratory of Bioinformatics and Systems' team.
