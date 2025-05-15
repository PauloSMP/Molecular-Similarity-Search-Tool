# Final Project CS50 Python – Molecular Similarity Search Tool, using PubRest API and RDKit Library
# Project Title: Molecular Similarity Search Tool
# My name is: Paulo Sérgio Marques Pereira
# My GitHub username is: PauloSMP
# My edx username is: paulo_m_pereira95
# My Country is: Portugal
# My City is: Coimbra
# Date(dd/mm/yyyy): 12/05/2025

# Importing general-use necessary libraries
import requests
import csv
import os
import sys
# import the RDKit library and the necessary modules
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, Draw


# Main Function
def main():
    
    #List of CIDs to build the database. 
    cids = [
        2405,
        2249,
        5992,
        3869,
        4171,
        4946,
        12456,
        1978,
        2369,
        59768,
        101338,
        37464,
        4828,
        71301,
        1983,
        2244,
        5288826,
        5284371,
        5284603,
        3778,
        3474,
        33613,
        6249,
        52946032,
        2764,
        4583,
        54675776,
        446598,
        441401,
    ]

    #Get the .sdf file from PubChem using the PubRest API.
    get_sdf_file(cids)

    #Get target molecule information from the user.
    smile_target_molecule, cid_target_molecule, name_target_molecule = (get_target_molecule())

    #Read the target molecule from the SMILE string input by the user. Generate a mol RDKit object.
    molecule = read_molecule(smile_target_molecule)

    #Read molecule database from the .sdf file. Generate a list of mol RDKit objects.
    molecules_db = read_molecules_db("PubChem_Dataset.sdf")

    #Compute morgan fingerprints for the target molecule.
    fingerprint_mol = compute_fingerprint(molecule)

    #Compute morgan fingerprints for the molecules in the database.
    fingerprints_mols_db = compute_fingerprints_db(molecules_db)

    #Calculate the Tanimoto similarity index between the target molecule and the mols in database.
    similarity = calculate_similarity(fingerprint_mol, fingerprints_mols_db)

    #Get the top n hits from the similarity index.
    n = int(input("Enter the number of top similar molecules to retrieve: ").strip())
    top_mols, top_scores = get_top_hits(similarity, molecules_db, n)

    #Draw the query molecule + top n molecules
    draw_molecules(molecule, top_mols, top_scores, cid_target_molecule, name_target_molecule)
    
    #Save the top n hits in a .csv file
    save_csv(top_mols, top_scores)


def get_sdf_file(cids):
    """
    Take the list of CIDs to build the database. 
    
    Get the .sdf file from PubChem using the PubRest API.
    """
    
    if os.path.exists("PubChem_Dataset.sdf"):
        print("File already exists. Skipping download.")
        return
    path = "PubChem_Dataset.sdf"
    with open(path, "ab") as f:
        for cid in cids:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF"
            response = requests.get(url)

            if response.status_code == 200:
                f.write(response.content)
                print(f"Downloaded CID: {cid} .sdf file")
            else:
                print(
                    f"Failed to download CID: {cid} .sdf file (status code: {response.status_code})"
                )


def get_target_molecule():
    
    """
    Get target molecule information from the user.
    Validate the user input.
    The user must enter the SMILE string, PubChem CID and name of the target molecule.
    
    Example:
    SMILE: CC(C)NCC(COC1=CC=C(C=C1)COCCOC(C)C)O
    SMILE --> beta-blocker known as Bisoprolol.
    
    Example:
    cid_target_molecule for bisoprolol = 2405 
    
    # Example (Common name / Iupac name):
    # Common name = Bisoprolol
    # IUPAC name  = 1-(propan-2-ylamino)-3-[4-(2-propan-2-yloxyethoxymethyl)phenoxy]propan-2-ol
    
    """

    smile_target_molecule = input("Enter SMILES string for your target molecule: ").strip()
    try:
        ms = Chem.MolFromSmiles(smile_target_molecule)
        if ms is None:
            raise ValueError(f"Invalid SMILES string: {smile_target_molecule}")
    except Exception as e:
        sys.exit(f"Error parsing SMILES string: {e}")

    cid_target_molecule = input("Enter PubChem CID for your target molecule: ").strip()
    if cid_target_molecule.isdigit():
        cid_target_molecule = int(cid_target_molecule)
    else:
        sys.exit("Invalid CID. Please enter a valid PubChem CID - (int).")

    name_target_molecule = input("Enter your target molecule name: ").strip()
    if not name_target_molecule:
        sys.exit("Empty mol name var. Please enter a name for your target molecule.")

    return smile_target_molecule, cid_target_molecule, name_target_molecule




def read_molecule(smile):
    
    """
    Read the target molecule from the SMILE string input by the user.
    
    Debugging:
    # print (len(molecules_db)) -->  = nº of molecules read/CIDs
    
    print (molecule)
    Result: <rdkit.Chem.rdchem.Mol object at 0x000002C55DC3EDC0>)
    
    """
    ms = Chem.MolFromSmiles(smile)
    return ms

    
    
def read_molecules_db(file_path):
    """
    Read molecule database from the .sdf file.
    
    Debugging:
    print (molecules_db[0]) 
    Result: <rdkit.Chem.rdchem.Mol object at 0x000002104EBFEE30>
    """
    suppl = Chem.SDMolSupplier(file_path)
    mol = [mol for mol in suppl if mol is not None]
    return mol


def compute_fingerprint(mol):
    """
    Compute fingerprints for the target molecule.
    Debugging:
    print(fingerprint_mol)
    result: <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x0000021C57D8CE40>
    """
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = fpgen.GetFingerprint(mol)
    return fp

def compute_fingerprints_db(mols_db):
    """
    Compute fingerprints for the molecules in the database.
    Debugging: Print the first fingerprint to check if it was computed correctly
    print(fingerprints_mols_db[0])
    result: <rdkit.DataStructs.cDataStructs.ExplicitBitVect object at 0x000001B22AB5FA70>
    """
    fps_list = []
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    for mol in mols_db:
        fps = fpgen.GetFingerprint(mol)
        fps_list.append(fps)
    return fps_list


def calculate_similarity(mol1, mol_list):
    """
    Calculate the Tanimoto similarity index between the target molecule and the database.
    Debugging:
    print(similarity)
    result: value of similarity between the target molecule and the database
    """
    similarity = []
    for mol in mol_list:
        sim = DataStructs.TanimotoSimilarity(mol1, mol)
        similarity.append(sim)
    return similarity


def get_top_hits(similarity, molecules_db, n):
    """
    Get the top n hits from the similarity index.
    Retun n top_mols and top_scores.
    Debugging:
    Print to check if it was computed correctly
    print(f"Top {n} molecules: {top_mols}")
    print(f"Top {n} scores: {top_scores}")
    """
    # Get indices for the top n hits
    top_hits = sorted(
        range(len(similarity)), key=lambda i: similarity[i], reverse=True
    )[:n]
    # Get the top molecules that correspond to top n hits
    top_mols = [molecules_db[i] for i in top_hits]
    # Get the top scores that correspond to top n hits
    top_scores = [similarity[i] for i in top_hits]
    # Return the top n hits, top molecules and top scores

    return top_mols, top_scores


def draw_molecules(mol, top_mols, top_scores, cid_target_molecule, name_target_molecule):
    """
    Draw the query molecule + top n molecules
    Debugging:
    check the image file "Top_n_hits.png" to see if it was generated correctly
    """
    # Add the query molecule with the top n hits
    molecules = [mol] + top_mols

    # Generate legends for the query and hits
    query_legend = f"Query_SMILE: {Chem.MolToSmiles(mol)}\nPubChemCid: {cid_target_molecule}\nName: {name_target_molecule}"
    hit_legends = [
        f'Hit_{i+1}_SMILE: {Chem.MolToSmiles(m)}\nSim: {top_scores[i]:.2f}\nPubChemCid:{m.GetProp("PUBCHEM_COMPOUND_CID")}\nIUPAC_name: {m.GetProp("PUBCHEM_IUPAC_NAME")}'
        for i, m in enumerate(top_mols)
    ]
    legends = [query_legend] + hit_legends

    # Draw the molecules (query + hits)
    img = Draw.MolsToGridImage(
        molecules, molsPerRow=3, subImgSize=(800, 800), legends=legends
    )
    img.save("Top_n_hits.png")


def save_csv(top_mols, top_scores):
    """
    Generate a CSV file with the top n hits.
    Debugging:
    check the CSV file "Top_n_hits.csv" to see if it was generated correctly
    """
    with open("Top_n_hits.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "PubChemCID", "Similarity", "IUPAC_NAME", "M-Formula"])
        for i, mol in enumerate(top_mols):
            pubchem_cid = mol.GetProp("PUBCHEM_COMPOUND_CID")
            iupac_name = mol.GetProp("PUBCHEM_IUPAC_NAME")
            m_formula = mol.GetProp("PUBCHEM_MOLECULAR_FORMULA")
            writer.writerow([i + 1, pubchem_cid, top_scores[i], iupac_name, m_formula])


if __name__ == "__main__":
    main()
