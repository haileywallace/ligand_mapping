# Ligand Mapping 
## Hailey Wallace 
Version 2.0 - May 10, 2022

 This project contains the ongoing development of a COMBS2-compatible application 
 that will automatically map out COMBS2 functional groups (CG). The output is a 
 'ligand map file' crucial for COMBS2 input and allows for the proper placement and 
 evaluation of van der Meers (vdMs).

 Convergent motifs for binding sites 2.0, or COMBS2, is an open-source package for 
 interfacial protein design created by Nick Polizzi. Download the full package [here](https://github.com/npolizzi/Combs2).

## Usage:

Install this application by either manual downloading via Github or by cloning this 
repository with the command:
```
git clone https://github.com/haileywallace/ligand_mapping.git
```
Keep this application up-to-date with:
```
git pull https://github.com/haileywallace/ligand_mapping.git
```

All of the required packages should be in the included environment. To create this environment, type:
```
conda env create -f ligand_map.yml 
```
To activate this environment, type:
```
conda activate ligand_map
```

Running this ligand_mapping python script requires an input PDB ligand. 
This ligand PDB must not have the receptor protein attached. There is an example
of this in the 'example' directory.
The default output for this script is a mapped ligand file named 'ligand.txt'
that can be input into COMBS2. To change the name of this output file, 
implement the --output option and enter a different filename. 
If you wish to output images of your ligand with the CGs highlighted, 
specify the --image option. Specifying nothing  will 
produce the following:
```
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
```

Example commandline input: 
```
python ligand_matcher.py -l 3fv2_ligand.pdb -o 3fv2_ligand.txt
```

