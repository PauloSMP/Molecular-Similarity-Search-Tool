# 🧪 Molecular Similarity Search Tool 👨‍🔬
    
#### 🎥 Video Demo:  <https://www.youtube.com/watch?v=ZhnIEjypgi0>
    
#### 🖹 Description: 
Welcome to my final project of CS50’s Introduction to Programming with Python. This project leverages the [PubChem PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) and the [RDKit](https://www.rdkit.org/docs/index.html) chemoinformatic library, to build a **Molecular Similarity Search Tool**.

#### 📚 Before Start:
You can find all the documentation regarding the [PubChem - PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest), including the [Tutorial](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial#section=How-PUG-REST-Works).

Regarding the [RDKit](https://www.rdkit.org/docs/index.html) library, the need documentation is all available, including:

- [An overview of the RDKit](https://www.rdkit.org/docs/Overview.html)
- [Installation](https://www.rdkit.org/docs/Install.html)
- [Getting Started with the RDKit in Python](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [The RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html)
- [Python API Reference](https://www.rdkit.org/docs/api-docs.html)

### 📚 [Molecular Similarity](https://www.drugdesign.org/chapters/molecular-similarity/)

**Molecular similarity** is a foundational concept in drug discovery and drug repurposing. It is based on the assumption that structurally similar molecules tend to exhibit similar physicochemical and biological properties. This principle is widely leveraged for identifying novel drug candidates and optimizing existing therapeutic agents.

This program implements molecular similarity screening using computational tools to identify structurally analogous compounds within a curated dataset. The process is centered on fingerprint-based similarity computation and Tanimoto scoring, applied to a set of compounds retrieved from PubChem.

---

## 🧬 Dataset and Target Molecule

The tool uses the [PubChem PUG REST API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest) to download a set of **.SDF** files corresponding to a curated list of **PubChem Compound Identifiers (CIDs)**, comprising compounds from three pharmacological classes:

- **Beta-blockers**
- **Analgesics**
- **Antibiotics**

The downloaded compounds are stored in a single `.sdf` file named **`PubChem_Dataset.sdf`**. 

This dataset serves as the screening library for identifying molecules similar to a predefined **target compound**, **`Bisoprolol`**, a beta-blocker.

The target compound is given by the user after being prompted to input a **SMILE** string* for the target molecule.

---

## 🔬 Similarity Computation Workflow

### 1. **Read Molecules**

The SMILE string of the **target molecule** is converted into an RDKit `Mol` object.

Concurrently, the tool reads the dataset of screening molecules from the `.sdf` file and converts each entry into `Mol` objects using RDKit.

- **SMILES Input**: The user must provide a valid **SMILES string**, obtainable from public chemical databases such as [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
  
- **SMILES → Mol Object**: Conversion is performed using  
  [`rdkit.Chem.MolFromSmiles`](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolFromSmiles)

- **SDF Dataset → Mol Objects**: Molecules in the `.sdf` file are parsed using  
  [`rdkit.Chem.SDMolSupplier`](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.SDMolSupplier)

**Molecular Properties of Interest:**
When reading a molecule from the `PubChem_Dataset.sdf`, each molecule is parsed into a RDKit object that  contains several annonated properties providedby PubChem.

To explore these properties you mayy use the following Python snippet:

```
#Example: Print all available properties for the first molecule
for mol in molecules_db[0:1]:
    print(list(mol.GetPropNames()))
```
This tool specifically utilizes the following properties:
- PUBCHEM_COMPOUND_CID: The unique identifier for the compound in the PubChem database.
- PUBCHEM_IUPAC_NAME: The systematic IUPAC name of the compound.
- PUBCHEM_MOLECULAR_FORMULA: The molecular formula of the compound.

These properties are leveraged in generating visual outputs (annotated molecule grids) and saved result files (CSV outputs).

**ℹ️ Note**: The user is also prompted to enter:
 - The **PubChem CID** of the target molecule  
 - The **name of the target molecule**  
 These fields are used for annotating the output `.png` image file.



### 2. **Calculate Molecule Fingerprints**

It is possible to encode molecular structures as binary vectors, capturing circular substructures around each atom. In this program we use radius-2 **Morgan Fingerprints**, using the module [`rdkit.Chem.rdFingerprintGenerator`](https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html#rdkit.Chem.rdFingerprintGenerator.GetMorganGenerator) and the **`GetFingerprint`** module.
- RETURNS: a **ExplicitBitVect containing fingerprint**.

### 3. **Tanimoto Coefficient**

The Tanimoto coefficient quantifies the similarity between two binary fingerprint vectors. It ranges from 0 to 1, where:
- **[0]** indicates no similarity.
- **[1]** indicates identical fingerprints.

The similarity between the target molecule and dataset compounds is calculated using:

- **Tanimoto Coefficient** with the module [`DataStructs.TanimotoSimilarity`](https://www.rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html#rdkit.DataStructs.cDataStructs.TanimotoSimilarity).


### 4. **Image and CSV Output Generation**

After computing similarity scores between the target molecule and all dataset entries, the top N most similar molecules are selected and stored in:

- **top_mols**: list of RDKit molecule objects

- **top_scores**: list of corresponding Tanimoto similarity values

Using the module [`rdkit.Chem.Draw.MolsToGridImage`](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html#rdkit.Chem.Draw.MolsToGridImage) to generate a grid image for each molecule in **top_mols**, to allow visual analysis and storage of the query results.

📊 `.png` Output:

**Top_n_hits.png**, includes:
- **Molecule Scaffold**
- **Similarity scores**
- **PubChem CID**
- **IUPAC Name**

📑 `.csv` Output:

**Top_n_hits.csv**, includes structured metadata such as:

- **Rank**
- **PubChem CID**
- **Similarity Score**
- **IUPAC Name**


