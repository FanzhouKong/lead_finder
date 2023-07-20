import os

from rdkit import Chem
import sys
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from toolsets.API_gets import smiles_to_name
sdf = sys.argv[1]
suppl = Chem.SDMolSupplier(sdf)
name = []
inchikey = []
formula = []
mono_mass = []
smiles = []
for mol in suppl:
    # smiles.append(Chem.MolToSmiles(mol))
    name.append(smiles_to_name(Chem.MolToSmiles(mol)))
    inchikey.append(Chem.MolToInchiKey(mol))
    formula.append(CalcMolFormula(mol))
    mono_mass.append(ExactMolWt(mol))
    smiles.append(Chem.MolToSmiles(mol))
    # break
output = pd.DataFrame(zip(name, inchikey, formula, mono_mass, smiles),
                      columns=['name', 'inchikey', 'formula', 'mono_mass', 'smiles']
                      )
base_dir = os.path.dirname(sys.argv[1])
file_base = os.path.basename(sys.argv[1]).split['.'][0]
output.to_csv(os.path.join(base_dir, file_base+'.csv'), index = False)