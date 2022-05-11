#!/usr/bin/python

# Hailey Wallace
# October 25, 2021
# May 10, 2022 updated
# Parses out all functional groups that COMBS recognizes from ligand of interest

################################################################################

# Import time to track script runtime
import time
from datetime import timedelta
start_time = time.time()

# Import os, sys, etc... + rdkit!
import sys, argparse, os, re, subprocess, collections
import pandas as pd
from csv import reader
from rdkit import Chem, RDConfig
from rdkit.Chem import Draw, ChemicalFeatures, rdchem, rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D


################################################################################
#os.remove("ligand_CG_coords.txt")
scriptdir = os.path.abspath(os.path.dirname(sys.argv[0]))
#print(scriptdir)

# Parser and commandline options
parser = argparse.ArgumentParser()
parser.add_argument('--ligand', '-l', type=str, required=True)
parser.add_argument('--coords', '-c', type=str, required=False) # update to remove option
parser.add_argument('--output', '-o', type=str, required=False)
parser.add_argument('--image', '-i', type=str, required=False, default=False)

if len(sys.argv)==1:
    help = """
*********************************************************
*       Ligand Matcher parses COMBS-compatible          *
*        functional groups from input ligand            *
*      __________________________________________       *
*                     options:                          *
*     --ligand, -l: input protonated                    *
*                    PDB ligand (required)              *
*     --ouput, -o: output name for ligand.txt file      *
*                   "ligand.txt" if not specified;      *
*                    "none" for no ligand.txt file      *
*     --image, -i: ouputs png images with               *
*                   highlighted functional groups       *
*                                                       *
*********************************************************
    """
    print(help)
    parser.exit()
args = parser.parse_args()

print("\n")   # blank line

# If options are specified:
if args.output == "none":
    print('Input ligand:', args.ligand, '\n'"No output text file.")
elif args.output == None:
    print('Input ligand:', args.ligand, '\n'"Output text file: ligand.txt")
    PATH = scriptdir + '/ligand.txt'
    if os.path.isfile(PATH):
        print("File already exists. Moving previous ligand.txt to ./old_files")
        PATH2 = scriptdir + '/old_files/'
        isExist = os.path.exists(PATH2)
        # Directory does exist
        if isExist:
            os.rename(PATH, scriptdir + '/old_files/ligand.txt')
        # Create a new directory because it does not exist
        if not isExist:
            os.mkdir(PATH2)
            os.rename(PATH, scriptdir + '/old_files/ligand.txt')
else:
    print('Input ligand:', args.ligand, '\n'"Output text file:", args.output)

################################################################################
# Functional groups that COMBS recognizes
################################################################################
conh2 = Chem.MolFromSmarts("C(=O)([NH2])")
bb_cco = Chem.MolFromSmarts("[*][C,c](=O)") # to not include COO, use "[N,n,C,c][C,c](O[!H])[!O]". C=O goes from both C-C=O and N-C=O
ph = Chem.MolFromSmarts("[c,C]1[c,C][c,C][c,C][c,C][c,C]1")
bb_cnh = Chem.MolFromSmarts("[C,c][N][H]")
ccn = Chem.MolFromSmarts("[c,C,n,N,o,O][C,c][NH3]")
ccoh = Chem.MolFromSmarts("[C,c,N,n,O,o][C,c][OH]")
coh = Chem.MolFromSmarts("[C,c,N,n]O[H]")
coo = Chem.MolFromSmarts("[*][C,c](=O)O")
csc = Chem.MolFromSmarts("[cH2,CH2]S[CH3]")
csh = Chem.MolFromSmarts("[C,c][S][H]")
gn = Chem.MolFromSmarts("C(=[NH1,NH2])([NH2,NH1,NH3])[NH1][C,c]")
hid = Chem.MolFromSmarts("[CH,cH]1[N,n][CH,cH][NH,nH][C,c]1")
hie = Chem.MolFromSmarts("[CH,cH]1[NH,nH][CH,cH][N,n][C,c]1")
hip = Chem.MolFromSmarts("[CH,cH]1[NH,nH][CH,cH][NH,nH][C,c]1")
indole = Chem.MolFromSmarts("c21cc[nH](c1cccc2)")
phenol = Chem.MolFromSmarts("[c,C]1[c,C][c,C][c,C][c,C][c,C]1[OH]")
isopropyl = Chem.MolFromSmarts("[CH3][c,C][CH3]")
pro = Chem.MolFromSmarts("[CH2R][CH2R][CH2R]")
ch3 = Chem.MolFromSmarts("[C!H3,c!H3][C,c][CH3]")

# Sort this list with the largest CG first and smallest last.
# This is so the CG ligand coverage number starts with the largest CG.
func_groups = [indole, phenol, ph, hip, hid, hie, gn, isopropyl, pro, ch3, csc, csh, conh2, coo, ccoh, coh, bb_cco, ccn, bb_cnh]

################################################################################
# Setting up file and SMARTS strings
################################################################################
input = scriptdir + '/' + args.ligand

# Import mol2 ligand file
#file = Chem.MolFromMol2File(input, removeHs=False) # Preserve the original Hs
file = Chem.MolFromPDBFile(input, removeHs=False)


# Make sure connection values at end of PDB/mol2 are correct

# SMARTS pattern of your input mol2
ligand = Chem.MolToSmarts(file)
try:
  #print("Ligand SMARTS pattern:", ligand)
  print("\n","Finding COMBS groups in ligand...")
except NameError:
  print("PDB is not defined")

################################################################################
# Find matching substructures, output to new outfile
################################################################################
# To return the variable name
def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

# Iterate functional group SMARTS patterns through your ligand
for combs_groups in func_groups:
    substruct = file.GetSubstructMatches(combs_groups)

    # assign new variable with old variable name, aka the functional group name
    pre_var= str(namestr(combs_groups, globals())).split("'")
    var = pre_var[1]

    if substruct:
        # print out the matching atom names from the substructure
        for match in substruct:
            match_atoms = []
            for i in match:
                n_i = file.GetAtomWithIdx(i).GetPDBResidueInfo().GetName()
                match_atoms.append(n_i)
            #print(match_atoms,match)

            # if -c "none", do not make text file
            if args.coords == "none":
                pass

            # make ligand coord file with default names
            elif args.coords == None:
                f = open(scriptdir + "/ligand_CG_coords.txt", "a")

                ligand_groups = []
                for match_lines in match_atoms:
                    ligand_groups.append(match_lines)

                # get residue ID of the first atom in the ligand PDB file
                #ligand_code = file.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName()

################################################################################
# Match the ligand groups to the combs groups
################################################################################

                # indole COMBS groups
                if var == "indole":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()

                    # indole from TRP
                    df_I = pd.DataFrame({'combs_group': ['CD2','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3']})
                    df_I.insert(loc=0, column="aa", value='PHE')
                    df_I.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_I.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_I['var'] = var
                    string_both_CG_lig_atoms = df_I.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # phenol COMBS groups
                if var == "phenol":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # phenol from TYR, rotation 1
                    df_Y1 = pd.DataFrame({'combs_group': ['CE1','CD1','CG','CD2','CE2','CZ','OH']})
                    df_Y1.insert(loc=0, column="aa", value='TYR')
                    df_Y1.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Y1.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Y1['var'] = var
                    df_Y1a = df_Y1.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_Y1a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # phenol from TYR, rotation 2
                    df_Y2 = pd.DataFrame({'combs_group': ['CE2','CD2','CG','CD1','CE1','CZ','OH']})
                    df_Y2.insert(loc=0, column="aa", value='TYR')
                    df_Y2.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Y2.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Y2['var'] = var
                    df_Y2a = df_Y2.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_Y2a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # phenyl COMBS groups
                if var == "ph":  # only 'ph' and not 'phenol'
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()

                    # phenyl for PHE, rotation 1
                    df_F1 = pd.DataFrame({'combs_group': ['CD1','CG','CD2','CE2','CZ','CE1']})
                    df_F1.insert(loc=0, column="aa", value='PHE')
                    df_F1.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_F1.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_F1['var'] = var
                    df_F1a = df_F1.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_F1a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # phenyl for PHE, rotation 2
                    df_F2 = pd.DataFrame({'combs_group': ['CE1','CD1','CG','CD2','CE2','CZ']})
                    df_F2.insert(loc=0, column="aa", value='PHE')
                    df_F2.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_F2.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_F2['var'] = var
                    df_F2a = df_F2.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_F2a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # phenyl for PHE, rotation 3
                    df_F3 = pd.DataFrame({'combs_group': ['CZ','CE1','CD1','CG','CD2','CE2']})
                    df_F3.insert(loc=0, column="aa", value='PHE')
                    df_F3.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_F3.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_F3['var'] = var
                    df_F3a = df_F3.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_F3a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # phenyl for PHE, rotation 4
                    df_F4 = pd.DataFrame({'combs_group': ['CE2','CZ','CE1','CD1','CG','CD2']})
                    df_F4.insert(loc=0, column="aa", value='PHE')
                    df_F4.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_F4.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_F4['var'] = var
                    df_F4a = df_F4.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_F4a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # phenyl for PHE, rotation 5
                    df_F5 = pd.DataFrame({'combs_group': ['CD2','CE2','CZ','CE1','CD1','CG']})
                    df_F5.insert(loc=0, column="aa", value='PHE')
                    df_F5.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_F5.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_F5['var'] = var
                    df_F5a = df_F5.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_F5a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # phenyl for PHE, rotation 6
                    df_F6 = pd.DataFrame({'combs_group': ['CG','CD2','CE2','CZ','CE1','CD1']})
                    df_F6.insert(loc=0, column="aa", value='PHE')
                    df_F6.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_F6.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_F6['var'] = var
                    df_F6a = df_F6.sort_values('combs_group')  ## sorting isnt working
                    string_both_CG_lig_atoms = df_F6a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # hip, double protonated His COMBS groups
                if var == "hip":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # hip for HIS, double-proton His
                    df_HIP = pd.DataFrame({'combs_group': ['CD2','NE2','CE1','ND1','CG']})
                    df_HIP.insert(loc=0, column="aa", value='HIP')
                    df_HIP.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_HIP.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_HIP['var'] = var
                    string_both_CG_lig_atoms = df_HIP.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # hid, delta protonated His COMBS groups
                if var == "hid":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # hid for HIS, delta protonated
                    df_HID = pd.DataFrame({'combs_group': ['CD2','NE2','CE1','ND1','CG']})
                    df_HID.insert(loc=0, column="aa", value='HID')
                    df_HID.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_HID.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_HID['var'] = var
                    string_both_CG_lig_atoms = df_HID.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # hie, epsilon protonated His COMBS groups
                if var == "hie":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # hie for HIS, epsilon protonated
                    df_HIE = pd.DataFrame({'combs_group': ['CD2','NE2','CE1','ND1','CG']})
                    df_HIE.insert(loc=0, column="aa", value='HIE')
                    df_HIE.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_HIE.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_HIE['var'] = var
                    string_both_CG_lig_atoms = df_HIE.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # guanidine COMBS groups
                if var == "gn":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # gn for ARG, rotation 1
                    df_R = pd.DataFrame({'combs_group': ['CZ','NH2','NH1','NE','CD']})
                    df_R.insert(loc=0, column="aa", value='ARG')
                    df_R.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_R.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_R['var'] = var
                    string_both_CG_lig_atoms = df_R.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # gn for ARG, rotation 2
                    df_R = pd.DataFrame({'combs_group': ['CZ','NH1','NH2','NE','CD']})
                    df_R.insert(loc=0, column="aa", value='ARG')
                    df_R.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_R.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_R['var'] = var
                    string_both_CG_lig_atoms = df_R.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # isopropyl COMBS groups
                if var == "isopropyl":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # isopropyl for LEU, rotation 1
                    df_L1 = pd.DataFrame({'combs_group': ['CD1','CG','CD2']})
                    df_L1.insert(loc=0, column="aa", value='LEU')
                    df_L1.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_L1.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_L1['var'] = var
                    df_L1a = df_L1.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_L1a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # isopropyl for VAL, rotation 1
                    df_V1 = pd.DataFrame({'combs_group': ['CG2','CB','CG1']})
                    df_V1.insert(loc=0, column="aa", value='VAL')
                    df_V1.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_V1.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_V1['var'] = var
                    df_V1a = df_V1.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_V1a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # isopropyl for LEU, rotation 2
                    df_L2 = pd.DataFrame({'combs_group': ['CD2','CG','CD1']})
                    df_L2.insert(loc=0, column="aa", value='LEU')
                    df_L2.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_L2.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_L2['var'] = var
                    df_L2a = df_L2.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_L2a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # isopropyl for VAL, rotation 2
                    df_V2 = pd.DataFrame({'combs_group': ['CG1','CB','CG2']})
                    df_V2.insert(loc=0, column="aa", value='VAL')
                    df_V2.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_V2.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_V2['var'] = var
                    df_V2a = df_V2.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_V2a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

                # proline COMBS groups (C-C-C)
                if var == "pro":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # pro for PRO, rotation 1
                    df_P1 = pd.DataFrame({'combs_group': ['CB','CG','CD']})
                    df_P1.insert(loc=0, column="aa", value='PRO')
                    df_P1.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_P1.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_P1['var'] = var
                    df_P1a = df_P1.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_P1a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # pro for PRO, rotation 1
                    df_P2 = pd.DataFrame({'combs_group': ['CD','CG','CB']})
                    df_P2.insert(loc=0, column="aa", value='PRO')
                    df_P2.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_P2.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_P2['var'] = var
                    df_P2a = df_P2.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_P2a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # Methyl COMBS groups
                if var == "ch3":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # ch3 for ALA
                    df_A = pd.DataFrame({'combs_group': ['C','CA','CB']})
                    df_A.insert(loc=0, column="aa", value='ALA')
                    df_A.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_A.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_A['var'] = var
                    string_both_CG_lig_atoms = df_A.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # ch3 for ILE
                    df_I = pd.DataFrame({'combs_group': ['CB','CG1','CD1']})
                    df_I.insert(loc=0, column="aa", value='ILE')
                    df_I.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_I.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_I['var'] = var
                    string_both_CG_lig_atoms = df_I.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

                # C-S-C COMBS groups
                if var == "csc":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # csc for MET
                    df_M = pd.DataFrame({'combs_group': ['CG','SD','CE']})
                    df_M.insert(loc=0, column="aa", value='MET')
                    df_M.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_M.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_M['var'] = var
                    string_both_CG_lig_atoms = df_M.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

                # C-SH COMBS groups
                if var == "csh":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # csh for CYS
                    df_C = pd.DataFrame({'combs_group': ['CB','SG','HG']})
                    df_C.insert(loc=0, column="aa", value='CYS')
                    df_C.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_C.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_C['var'] = var
                    string_both_CG_lig_atoms = df_C.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

                # CONH2 COMBS groups
                if var == "conh2":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # conh2 for ASN
                    df_N = pd.DataFrame({'combs_group': ['CG','OD1','ND2']})
                    df_N.insert(loc=0, column="aa", value='ASN')
                    df_N.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_N.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_N['var'] = var
                    string_both_CG_lig_atoms = df_N.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # conh2 for GLN
                    df_Q = pd.DataFrame({'combs_group': ['CD','OE1','NE2']})
                    df_Q.insert(loc=0, column="aa", value='GLN')
                    df_Q.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Q.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Q['var'] = var
                    string_both_CG_lig_atoms = df_Q.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

                # coo COMBS groups
                if var == "coo":
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # coo for ASP, rotation 1
                    df_D1 = pd.DataFrame({'combs_group': ['CB','CG','OD2','OD1']})
                    df_D1.insert(loc=0, column="aa", value='ASP')
                    df_D1.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_D1.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_D1['var'] = var
                    df_D1a = df_D1.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_D1a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # coo for GLU, rotation 1
                    df_E1 = pd.DataFrame({'combs_group': ['CG','CD','OE2','OE1']})
                    df_E1.insert(loc=0, column="aa", value='GLU')
                    df_E1.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_E1.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_E1['var'] = var
                    df_E1a = df_E1.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_E1a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                    # coo for ASP, rotation 2
                    df_D2 = pd.DataFrame({'combs_group': ['CB','CG','OD1','OD2']})
                    df_D2.insert(loc=0, column="aa", value='ASP')
                    df_D2.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_D2.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_D2['var'] = var
                    df_D2a = df_D2.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_D2a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # coo for GLU, rotation 2
                    df_E2 = pd.DataFrame({'combs_group': ['CG','CD','OE1','OE2']})
                    df_E2.insert(loc=0, column="aa", value='GLU')
                    df_E2.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_E2.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_E2['var'] = var
                    df_E2a = df_E2.sort_values('combs_group')
                    string_both_CG_lig_atoms = df_E2a.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # ccoh COMBS groups
                if var == "ccoh":
                    # get residue ID of the first atom in the ligand PDB file
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()

                    # ccoh for SER
                    df_S = pd.DataFrame({'combs_group': ['CA','CB','OG']})
                    df_S.insert(loc=0, column="aa", value='SER')
                    df_S.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_S.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_S['var'] = var
                    string_both_CG_lig_atoms = df_S.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # ccoh for THR
                    df_T = pd.DataFrame({'combs_group': ['CA','CB','OG1']})
                    df_T.insert(loc=0, column="aa", value='THR')
                    df_T.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_T.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_T['var'] = var
                    string_both_CG_lig_atoms = df_T.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                #  coh COMBS groups
                if var == "coh":
                    # get residue ID of the first atom in the ligand PDB file
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()

                    # coh for SER
                    df_S = pd.DataFrame({'combs_group': ['CB','OG','HG']})
                    df_S.insert(loc=0, column="aa", value='SER')
                    df_S.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_S.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_S['var'] = var
                    string_both_CG_lig_atoms = df_S.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # coh for THR
                    df_T = pd.DataFrame({'combs_group': ['CB','OG1','HG1']})
                    df_T.insert(loc=0, column="aa", value='THR')
                    df_T.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_T.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_T['var'] = var
                    string_both_CG_lig_atoms = df_T.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)
                    print("", file=f)

                # bb_cco COMBS groups
                if var == "bb_cco":
                    # get residue ID of the first atom in the ligand PDB file
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # bb_cco for ALA
                    df_Abb = pd.DataFrame({'combs_group': ['CA','C','O']})
                    df_Abb.insert(loc=0, column="aa", value='ALA')
                    df_Abb.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Abb.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Abb['var'] = var
                    string_both_CG_lig_atoms = df_Abb.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # bb_cco for GLY
                    df_Gbb = pd.DataFrame({'combs_group': ['CA','C','O']})
                    df_Gbb.insert(loc=0, column="aa", value='GLY')
                    df_Gbb.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Gbb.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Gbb['var'] = var
                    string_both_CG_lig_atoms = df_Gbb.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # bb_cco for PRO
                    df_Pbb = pd.DataFrame({'combs_group': ['CA','C','O']})
                    df_Pbb.insert(loc=0, column="aa", value='PRO')
                    df_Pbb.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Pbb.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Pbb['var'] = var
                    string_both_CG_lig_atoms = df_Pbb.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

                # ccn COMBS groups
                if var == "ccn":
                    # get residue ID of the first atom in the ligand PDB file
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # ccn for LYS
                    df_K = pd.DataFrame({'combs_group': ['CD','CE','NZ']})
                    df_K.insert(loc=0, column="aa", value='LYS')
                    df_K.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_K.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_K['var'] = var
                    string_both_CG_lig_atoms = df_K.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

                # ccn COMBS groups
                if var == "bb_cnh":
                    # get residue ID of the first atom in the ligand PDB file
                    for i in match[:1]:
                        ligand_code=file.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueName()
                    # bb_cnh for LYS
                    df_Kbb = pd.DataFrame({'combs_group': ['CA','N','H']})
                    df_Kbb.insert(loc=0, column="aa", value='LYS')
                    df_Kbb.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Kbb.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Kbb['var'] = var
                    string_both_CG_lig_atoms = df_Kbb.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # bb_cnh for ALA
                    df_Abb = pd.DataFrame({'combs_group': ['CA','N','H']})
                    df_Abb.insert(loc=0, column="aa", value='ALA')
                    df_Abb.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Abb.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Abb['var'] = var
                    string_both_CG_lig_atoms = df_Abb.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    # bb_cnh for GLY
                    df_Gbb = pd.DataFrame({'combs_group': ['CA','N','H']})
                    df_Gbb.insert(loc=0, column="aa", value='GLY')
                    df_Gbb.insert(loc=0, column="ligand_code", value=ligand_code)
                    df_Gbb.insert(loc=1, column="ligand_atoms", value=ligand_groups)
                    df_Gbb['var'] = var
                    string_both_CG_lig_atoms = df_Gbb.to_string(index=False, header=False)
                    print(string_both_CG_lig_atoms, file=f)

                    print("", file=f)

            # make ligand coord file for specified filename
            # come back and edit this
            else:
                output = args.coords
                f = open(scriptdir + "/" + output, "a")
                for match_lines in match_atoms:
                    print(match_lines, var, file=f)
                f.close()
            f.close()

    # Count COMBS functional groups in ligand
    if len(substruct) != 0:
        print(var, "found", len(substruct), "matches.")

################################################################################
# Add group numbers to each chunk of text (seperated by newline)
# depicting combs groups
################################################################################

counter = 1
with open('ligand_CG_coords.txt') as infile:

    with open('out.txt', 'w') as outfile:

        for lines in infile:
            if lines == '\n':
                counter += 1
                print('\n', file=outfile)
            else:
                new_line = lines.rstrip('\n') +"\t" + str(counter)
                print(new_line, file=outfile)
    outfile.close()
infile.close()

################################################################################
# Matches these combs groups with other combs groups that share 2 atoms. (coverage number)
# Sorts groups starting with the largest combs group and descending in size
# Smaller combs groups will match into the larger combs groups' coverage numbers
################################################################################

# short bash script to place the unique atom names with each corresponding COMBS group number
subprocess.run(['bash','./support.sh'])

# dictionary
combs_group = {}
atoms = collections.defaultdict(list)

# assign combs group coverage number first pass
with open("sorted.csv") as h:
    for lines in h:
        number,items = lines.split(",")
        if items:
            combs_group[number] = items.split()
            for item in combs_group[number]:
                atoms[item].append(number)

# if there are 2 or more atoms that match another combs group,
# then assign new group the previous combs coverage number
with open('out_new.txt', 'w') as outfile_new:

    for atom,items in combs_group.items():
        for count in items:
            if len(atoms[count]) >= 2:
                first_number = f'{atom},{" ".join(items)},{atoms[count][0]}'
                print(first_number,file=outfile_new)
                break
            else:
                prior_number = f'{atom},{" ".join(items)}'
                print(prior_number,file=outfile_new)

outfile_new.close()
h.close()

################################################################################
# Output ligand.txt file and remove temp files
################################################################################

if args.output == None:
    with open('out3.txt', 'w') as final_file:
        with open('out.txt') as myfile:
            for line in myfile:
                if line == '\n':
                    print("", file=final_file)

                else:
                    with open('out_new.txt', 'r') as read_obj:
                        csv_reader = reader(read_obj)
                        for row in csv_reader:
                            group_number = [ row[0] ]
                            coverage_number = [ row[-1] ]
                            coverage_number_string = str(*coverage_number)
                            coverage_number_tab = "\t", coverage_number_string
                            line2 = line.replace('\n', '')
                            c = line2.split("\t")
                            v = [ c[-1] ]
                            if group_number == v:
                                final_line = "  ".join([line2, str(*coverage_number)])
                                print(final_line, file=final_file)
        myfile.close()
    final_file.close()

else:
    file = str(args.output)
    with open('out3.txt', 'w') as file_unique:
        with open('out.txt') as myfile2:
            for line in myfile2:
                if line == '\n':
                    print("", file=file_unique)

                else:
                    with open('out_new.txt', 'r') as read_obj:
                        csv_reader = reader(read_obj)
                        for row in csv_reader:
                            group_number = [ row[0] ]
                            coverage_number = [ row[-1] ]
                            coverage_number_string = str(*coverage_number)
                            coverage_number_tab = "\t", coverage_number_string
                            line2 = line.replace('\n', '')
                            c = line2.split("\t")
                            v = [ c[-1] ]
                            if group_number == v:
                                final_line = "  ".join([line2, str(*coverage_number)])
                                print(final_line, file=file_unique)
        myfile2.close()
    file_unique.close()

    with open ('out3.txt', 'r') as i, open (file, 'w') as respace:
        for line in i:
            respace.write(re.sub('[\t ]+',' ', line))

with open ('out3.txt', 'r') as i, open ('ligand.txt', 'w') as respace:
    for line in i:
        respace.write(re.sub('[\t ]+',' ', line))

os.remove("ligand_CG_coords.txt")
os.remove("out.txt")
os.remove("out_new.txt")
os.remove("sorted.csv")

################################################################################
# Optional -- Image creation --
# print substructure match to a JPG of your ligand structure only if -i is specified
################################################################################
# If -i is specified
if args.image:
    # Image options
    image = args.image
    IPythonConsole.drawOptions.addAtomIndices = True
    IPythonConsole.molSize = 300,300

    # Iterate over substructure matches
    for combs_groups_img in func_groups:
        substruct2 = file.GetSubstructMatches(combs_groups_img)
        if len(substruct2) != 0:
            for group in substruct2:
                i = 0
                while os.path.exists(scriptdir + "/" + image+"%s.png" % i):
                    i += 1

                # Make a label for the images using functional group variable names
                name = str(namestr(combs_groups_img, globals())).split("'")
                label = str(name[1])

                # Atom indices for highlight
                highlights = list(group)
                print("Atom indices",highlights,"highlighted in",image+"%s.png" % i,"for",label)

                # Drawing the image -- still has Hydrogens because of the H atom indices used in highlight
                d = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
                d.drawOptions().addAtomIndices = True
                d.DrawMolecule(file, highlightAtoms=highlights, legend=label)
                d.FinishDrawing()
                d.WriteDrawingText(scriptdir + "/" + image+"%s.png" % i)

################################################################################
# Parse the previously output file for matching atom names
# Precise atom names are crucial for COMBS
################################################################################
# Extract coordinates and match with group name

if args.coords == "none":
    elapsed_time_secs = time.time() - start_time
    msg = "...Finished in %s seconds." % timedelta(seconds=round(elapsed_time_secs))
    print(msg)

elif args.output =="none":
    elapsed_time_secs = time.time() - start_time
    msg = "...Finished in %s seconds." % timedelta(seconds=round(elapsed_time_secs))
    print(msg)

#elif args.coords == None:
#    if args.output == None:
#        process = subprocess.run(['bash','./ligand_mapping_CG_files/ligand_support.sh','-l',input,'-i','ligand_CG_coords.txt'])
#        elapsed_time_secs = time.time() - start_time
#        msg = "Finished making ligand.txt in %s seconds" % timedelta(seconds=round(elapsed_time_secs))
#        print(msg)
#    else:
#        process = subprocess.run(['bash','./ligand_mapping_CG_files/ligand_support.sh','-l',input,'-i','ligand_CG_coords.txt','-o',args.output])
#        elapsed_time_secs = time.time() - start_time
#        msg = "%s seconds" % timedelta(seconds=round(elapsed_time_secs))
#        print('Finished making', args.output, msg)

#else:
#    if args.output == None:
#        process = subprocess.run(['bash','./ligand_mapping_CG_files/ligand_support.sh','-l',input,'-i',args.coords])
#        elapsed_time_secs = time.time() - start_time
#        msg = "Finished making ___ in %s seconds" % timedelta(seconds=round(elapsed_time_secs))
#        print(msg)
#    else:
#        process = subprocess.run(['bash','./ligand_mapping_CG_files/ligand_support.sh','-l',input,'-i',args.coords,'-o',args.output])
#        elapsed_time_secs = time.time() - start_time
#        msg = "%s seconds" % timedelta(seconds=round(elapsed_time_secs))
#        print('Finished making', args.output, msg)

elapsed_time_secs = time.time() - start_time

if __name__ == "__main__":
    print("\n")
    print("ligand_matcher.py was ran directly and finished in %s seconds." % timedelta(seconds=round(elapsed_time_secs)))
else:
    print("Importing ligand_matcher...")


print("Make sure your model is protonated, your atom connectivity is correct, and atom names are unique in your input PDB file!","\n")
################################################################################
# Read this if you need to edit the code to add more COMBS-groups!
# Notations in progress...
################################################################################
