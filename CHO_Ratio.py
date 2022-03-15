''' Module to calculate C/O and C/H ratio '''

from rdkit import Chem

# O/C ratio
def O_C_ratio(information):
    OCratio = information[7][8:-2].split(", ")
    for i in range(len(OCratio)):
        o = OCratio[i].count('O')
        c = OCratio[i].count('C') + OCratio[i].count('c')
        OCratio[i] = o / c if c != 0 else 'no C'
    OCratio = OCratio + [0] + [0]
    return OCratio


# H/C ratio
def H_C_ratio(information):
    HCratio = information[7][8:-2].split(", ")
    for i in range(len(HCratio)):
        my_smiles_string = HCratio[i][1:-1]
        my_mol = Chem.MolFromSmiles(my_smiles_string)
        try:
            Chem.AddHs(my_mol)
        except:
            HCratio[i] = 'error'
        else:
            my_mol_with_explicit_h = Chem.AddHs(my_mol)
            h = my_mol_with_explicit_h.GetNumAtoms() - my_mol_with_explicit_h.GetNumHeavyAtoms()
            c = my_smiles_string.count('C') + my_smiles_string.count('c')
            HCratio[i] = h / c if c != 0 else 'no C'
    HCratio = HCratio + [' '] + [' ']
    return HCratio

# error in index 549, 'O=\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tC1C=C(C)C(C)O1'
