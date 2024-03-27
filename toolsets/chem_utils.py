from rdkit import Chem
import re
import chemparse
from molmass import Formula
from pubchempy import Compound, get_compounds
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Descriptors import ExactMolWt
import cirpy# reference plz see https://www.resources.aropha.com/blog/get-chemical-smiles-by-cas-or-name/
import numpy as np
def desalter(input):
    if input != input:
        return np.NAN
    if is_mol(input) == False:
        smile = everything_to_smiles(input)
        mol = Chem.MolFromSmiles(smile)
    else:
        mol = input
    components = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(components)==1:
        if Chem.GetFormalCharge(components[0])==1:
            uncharged_smiles = remove_acidic_hydrogen(components[0])
        else:
            uncharged_smiles = Chem.MolToSmiles(components[0])
        return uncharged_smiles
    else:
        n_atoms = np.zeros(len(components))
        counter = 0
        for component in components:
            n_atoms[counter]=component.GetNumAtoms()
            counter = counter+1
        idx = np.argmax(n_atoms)
        # print(idx)
        charged = components[np.argmax(n_atoms)]
        un = rdMolStandardize.Uncharger()
        uncharged = un.uncharge(charged)
        if Chem.GetFormalCharge(uncharged)==1:
            uncharged_smiles = remove_acidic_hydrogen(uncharged)
        else:
            uncharged_smiles = Chem.MolToSmiles(uncharged)
        return( uncharged_smiles   )
def remove_acidic_hydrogen(molecule, is_smiles = False):
    # Convert the SMILES string to a RDKit molecule object
    if is_smiles==True:
        smiles = molecule
        molecule = Chem.MolFromSmiles(molecule)
    else:
        smiles = Chem.MolToSmiles(molecule)


    # Define the SMARTS pattern for carboxylic acids (includes the hydrogen in the hydroxyl group)
    carboxylic_acid_smarts = 'C(=O)[OH]'

    # Create a query molecule from the SMARTS pattern
    query = Chem.MolFromSmarts(carboxylic_acid_smarts)

    # Find substructures that match the query (carboxylic acids)
    matches = molecule.GetSubstructMatches(query)

    if not matches:
        # print("No carboxylic acid group found.")
        return smiles  # Return the original SMILES if no carboxylic acid group is found
    editable_mol = Chem.RWMol(molecule)

    # Assuming only one carboxylic group needs to be modified,
    # and focusing on the first match
    for match in matches:
        # The oxygen atom in the OH group is the last second in the matched pattern
        oxygen_idx = match[-1]

        # Get the oxygen atom
        oxygen_atom = editable_mol.GetAtomWithIdx(oxygen_idx)

        # Set the formal charge of the oxygen atom to -1
        oxygen_atom.SetFormalCharge(-1)

        # Set the implicit hydrogen count of the oxygen atom to 0
        # Assuming there's only one hydrogen bonded which we want to remove
        oxygen_atom.SetNumExplicitHs(0)

        # Break after the first modification, assuming only one modification is needed
        break

    # Convert back to a molecule
    modified_mol = editable_mol.GetMol()

    # Convert the modified molecule back to SMILES without sanitization
    modified_smiles = Chem.MolToSmiles(modified_mol)

    return modified_smiles
def everything_to_formula(input):
    if input != input:
        return np.NAN
    smiles = everything_to_smiles(input)
    mol = Chem.MolFromSmiles(smiles)
    formula_temp = CalcMolFormula(mol)
    formula = standarize_formula(formula_temp)
    return(formula_temp)
def standarize_formula(formula):
    if formula[-1] in ['+', '-']:
        main_body = formula[0:-1]
        polarity = formula[-1]
    else:
        main_body = formula
        polarity = ''
    formula_s = Formula(main_body).formula
    formula_s = formula_s+polarity
    return formula_s
def everything_to_smiles(input):
    if is_mol(input):
        smiles = Chem.MolToSmiles(input)
    elif is_smiles(input):
        smiles = input
    elif is_cas_number(input):
        smiles = cas_to_smiles(input)
    elif is_inchikey(input):
        smiles = inchikey_to_smiles(input)
    else:
        smiles = name_to_smiles(input)
    return(smiles)

#below are individual building blocks to smiles
def inchikey_to_smiles(inchikey):
    cc = get_compounds(inchikey, 'inchikey')
    if len(cc)>0:
        return (cc[0].isomeric_smiles)
    else:
        cc = get_compounds(inchikey[0:14], 'inchikey')
        if len(cc)>0:
            return (cc[0].isomeric_smiles)
        else:
            return (np.NAN)
def cas_to_smiles(cas):
    smile = cirpy.resolve(cas, 'smiles')
    if smile is None:
        smile = np.NAN
    return(smile)
def name_to_smiles(name):
    cc = get_compounds(name, 'name')
    if len(cc)>0:
        return (cc[0].isomeric_smiles)
    else:
        return (np.NAN)

def everything_to_image(molecule, savepath):
    from rdkit import Chem
    from rdkit.Chem import Draw
    if is_mol(molecule):
    # Create an RDKit molecule object
        mol = molecule

    elif is_smiles(molecule):
        # print('ttt')
        mol = Chem.MolFromSmiles(molecule)
    else:
        smiles = everything_to_smiles(molecule)
        mol = Chem.MolFromSmiles(smiles)
    # Generate the image of the molecule
    img = Draw.MolToImage(mol)
        # Save the image to a file
    img.save(savepath)
def calculate_precursormz(formula, adduct, if_smiles = False):
    if adduct[-1]=='+':
        adduct_polarity = 1
    elif adduct[-1]=='-':
        adduct_polarity = -1
    else:
        print('cannot determine adduct polarity!')
        return (np.NAN)
    if adduct[0]=='[':
        adduct_part = re.split(r'\[|]', adduct)[1]
        ind = re.split(r'(\+|-)', adduct_part)
    else:
        adduct_part = adduct
        ind = re.split(r'(\+|-)', adduct_part)
    coef,a  = parse_adduct(ind[0])

    if if_smiles == True:
        mol = Chem.MolFromSmiles(formula)
        # max_c_length = find_longest_element_chain(mol, 'C')
        formula = CalcMolFormula(mol)
    if formula[-1] in ['+','-'] and adduct_part not in ['M', 'Cat']:
        return 0
    elif formula[-1] not in ['+', '-'] and adduct_part in ['M', 'Cat']:
        return 0
    elif formula[-1] in ['+', '-'] and adduct_part in ['M', 'Cat']:
        if formula[-1]==adduct[-1]:
            accurate_mass = Formula(formula[:-1]).isotope.mass
            accurate_mass = accurate_mass-adduct_polarity*0.00054857990924
            return accurate_mass
        else:
            return 0

    if 'Hac' in adduct:
        adduct = adduct.replace("Hac", "C2H4O2")
    if 'FA' in adduct:
        adduct = adduct.replace("FA", "CH2O2")
    if 'DMSO' in adduct:
        adduct = adduct.replace("DMSO", "C2H6OS")

    accurate_mass = Formula(formula).isotope.mass

    accurate_mass = accurate_mass*coef


    for i in range(1, len(ind)):
        if (ind[i] not in ['+', '-']) and len(ind[i])>0:
            coef, a = parse_adduct(ind[i])
            if ind[i-1]=='+':
                polarity = 1
            elif ind[i-1]=='-':
                polarity = -1
            else:
                polarity = 0

            accurate_mass = accurate_mass+polarity*coef*Formula(a).isotope.mass
    accurate_mass = accurate_mass-adduct_polarity*0.00054857990924
    return accurate_mass
def parse_adduct(s):
    # Initialize an index variable to keep track of where the digits end
    index = 0
    coef = 1
    # Loop through each character in the string
    for char in s:
        # If the character is a digit, increment the index
        if char.isdigit():
            index += 1
        # If a non-digit character is found, break out of the loop
        else:
            break
    if index >0:
        coef = int(s[:index])
    # Return the part of the string without the leading integers
    return coef,s[index:]

#below are is_ section
def is_inchikey(string):
    # Define the regex pattern for InChIKeys
    pattern = r'^[A-Z]{14}-[A-Z]{10}-[A-Z0-9]$'

    # Use re.match to check if the pattern matches the entire string
    if re.match(pattern, string):
        return True
    else:
        return False

def is_mol(obj):
    return isinstance(obj, Chem.rdchem.Mol)
def is_smiles(smiles_string):
    # Attempt to create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)

    # If the molecule object is created successfully, the SMILES string is valid
    if mol is not None:
        return True
    else:
        # If the molecule object is None, the SMILES string is invalid
        return False
def is_cas_number(string):
    # Regex pattern for CAS numbers: one or more digits, followed by a hyphen, followed by two or more digits,
    # followed by a hyphen, and ending with a single digit
    pattern = r'^\d+-\d{2,}-\d$'

    # Check if the given string matches the pattern
    if re.match(pattern, string):
        return True
    else:
        return False
def everything_to_mw(mol):
    if is_mol(mol)==False:
        smiles = everything_to_smiles(mol)
        mol = Chem.MolFromSmiles(smiles)
    return(ExactMolWt(mol))
def is_formula(s):
    # Regular expression to match chemical formulas
    # Starts with an uppercase letter, optionally followed by a lowercase letter (for two-letter elements)
    # Optionally followed by a number (for the count of atoms)
    # This pattern repeats throughout the string
    pattern = r'^([A-Z][a-z]?\d*)+$'

    # Match the entire string against the pattern
    match = re.fullmatch(pattern, s)

    # If there's a match, the string is a valid chemical formula
    return bool(match)