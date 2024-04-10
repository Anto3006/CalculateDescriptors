from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
import subprocess
import os
import pandas as pd


def smiles2sdf(smiles, output_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    writer = Chem.SDWriter(output_file)
    writer.write(mol)
    writer.close()



def calculateConjugatedDescriptors(smiles):
    for i, smiles in enumerate(smiles):
        output_file = f'molecule_{i+1}.sdf'
        smiles2sdf(smiles, output_file)
    subprocess.call("Rscript conjugaR.R", shell=True)
    files = os.listdir()
    sdf_files = [file for file in files if file.split(".")[-1] == ".sdf"]
    for file in sdf_files:
        os.remove(file)
    descriptors_file_path = 'conjugated_descriptors.csv'
    descriptors = pd.read_csv(descriptors_file_path)
    os.remove(descriptors_file_path)
    return descriptors

calculateConjugatedDescriptors(["C","CC"])