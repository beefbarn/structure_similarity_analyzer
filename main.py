from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, DataStructs, Draw

# Load the SDF file
sdf_file = "sample_molecules.sdf"  # Change if needed
suppl = Chem.SDMolSupplier(sdf_file)

# Convert molecules in SDF to RDKit format
database = []
mol_dict = {}  # Store molecules by name for later retrieval
for mol in suppl:
    if mol is not None:
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"
        fingerprint = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024).GetFingerprint(mol)
        database.append((name, fingerprint))
        mol_dict[name] = mol  # Save molecule for later use

# Define target molecule
target_smiles = "C1=CC=CC=C1"  # Example: Ethanol CCO
target_mol = Chem.MolFromSmiles(target_smiles)
target_fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024).GetFingerprint(target_mol)

# Compare similarity
results = []
for name, db_fp in database:
    similarity = DataStructs.TanimotoSimilarity(target_fp, db_fp)
    results.append((name, similarity))

# Sort by highest similarity
results.sort(key=lambda x: x[1], reverse=True)

# Print results
print(f"Target Molecule: {target_smiles}\n")
print("Similar Molecules in Database:")
for name, sim in results:
    print(f"{name}: Similarity {sim:.3f}")

# Generate image for the top match
if results:
    top_match = results[0][0]  # Get name of most similar molecule
    mol = mol_dict.get(top_match)  # Retrieve molecule from dictionary
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        img.show()  # Display the image
        img.save(f"{top_match}.png")  # Save the image
        print(f"Saved structure image as {top_match}.png")
