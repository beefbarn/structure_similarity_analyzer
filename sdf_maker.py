from rdkit import Chem
from rdkit.Chem import AllChem

# List of sample molecules (name, SMILES)
smiles_list = [
    ("Ethanol", "CCO"),
    ("Benzene", "C1=CC=CC=C1"),
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("Glucose", "C(C1C(C(C(C(O1)O)O)O)O)O"),
    ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("Morphine", "CC1CCC2C3C1CCN2CC4=C3C=CC(=O)O4"),
    ("Nicotine", "CN1CCCC1c2cccnc2"),
    ("L-DOPA", "N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O"),
]

# Create an SDF writer
sdf_filename = "sample_molecules.sdf"
writer = Chem.SDWriter(sdf_filename)

# Convert SMILES to RDKit molecules and write to SDF
for name, smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol.SetProp("_Name", name)  # Add molecule name
        AllChem.Compute2DCoords(mol)  # Generate 2D coordinates
        writer.write(mol)

writer.close()
print(f"SDF file saved as {sdf_filename}")