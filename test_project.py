# This file is intended to test project.py with pytest.

from project import (
    read_molecule,
    read_molecules_db,
    compute_fingerprint,
    compute_fingerprints_db,
    calculate_similarity,
)

# import the RDKIT library and the necessary modules
from rdkit import Chem, DataStructs
import pytest


def test_read_molecule_pass():

    mol = read_molecule("CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC")
    assert isinstance(mol, Chem.rdchem.Mol)


def test_read_molecule_error():

    mol = read_molecule("NOT_A_SMILE")
    # Will raise a TypeError because the SMILE string is invalid
    # EX:
    # [23:00:04] SMILES Parse Error: syntax error while parsing: NOT_A_SMILE
    # [23:00:04] SMILES Parse Error: check for mistakes around position 3:
    # [23:00:04] NOT_A_SMILE
    # [23:00:04] ~~^
    # [23:00:04] SMILES Parse Error: Failed parsing SMILES 'NOT_A_SMILE' for input: 'NOT_A_SMILE'
    with pytest.raises(TypeError):
        read_molecule(mol)


def test_read_molecules_db_pass():
    path = "PubChem_Dataset.sdf"
    mols_db = read_molecules_db(path)
    # checkif the output is a list
    assert isinstance(mols_db, list)
    # check if the list is not empty
    assert len(mols_db) > 0
    # Check if all elements in the list are rdkit.Chem.rdchem.Mol objects
    for mol in mols_db:
        assert isinstance(mol, Chem.rdchem.Mol)


def test_read_molecules_db_error():
    path = "invalid_path.sdf"
    with pytest.raises(OSError):
        read_molecules_db(path)


def test_compute_fingerprint_pass():

    mol = read_molecule("CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC")
    test_var = compute_fingerprint(mol)

    assert isinstance(test_var, DataStructs.ExplicitBitVect)


def test_compute_fingerprints_db_pass():
    mols_db = read_molecules_db("PubChem_Dataset.sdf")
    fps_list = compute_fingerprints_db(mols_db)

    # Check if the output is a list
    assert isinstance(fps_list, list)
    # Check if the list is not empty
    assert len(fps_list) > 0
    # Check if all elements in the list are rdkit.DataStructs.ExplicitBitVect objects
    for fingerprint in fps_list:
        assert isinstance(fingerprint, DataStructs.ExplicitBitVect)


def test_calculate_similarity_pass():

    mol = read_molecule("CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC")
    mol1 = compute_fingerprint(mol)

    mols_db = read_molecules_db("PubChem_Dataset.sdf")
    fps_list = compute_fingerprints_db(mols_db)

    similarity = calculate_similarity(mol1, fps_list)

    # Check if the output is a list
    assert isinstance(similarity, list)
    # Check if len of sim, list == len of fps_list
    assert len(similarity) == len(fps_list)
    # Check if all elements in the list are floats
    for sim in similarity:
        assert isinstance(sim, float)
