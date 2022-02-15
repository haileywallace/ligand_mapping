#!/bin/bash
# Hailey Wallace


# Command line options
while getopts ":hi:l:o:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # Enter the coords file
         i=$OPTARG;;
      l) # Enter the ligand PDB file
         l=$OPTARG;;
      o) # Enter the ouput file name
         o=$OPTARG;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

echo "Making ligand CG file..."

# Change the ligand coordinate file into a tab delimited tmp_LG.csv file
while read -r p; do
  group=$(echo $p | awk -F'(' '{print $1}')
  echo $p | awk -F'|' '{print $2}' | sed -e $'s/;/\\\n/g' | \
  awk -F '\t' -v group2=$group '{print $0","group2}'| \
  sed -e $'s/(/\\\n/g' | sed -e 's/)//' >> tmp_LG.csv
done < $i

# If no substructures are found, exit script here
[[ ! -f tmp_LG.csv ]] && echo "ERROR: No combs groups found." && exit 1

# Read the coordinate file by chunks separated by newline; feed into array
i=1
s=1
declare -a arr
while read -r line
do
    # Increase counter (i) if there is an empty line. Set s to 1.
    [[ $line == "" ]] && ((i++)) && s=1 && continue
    # If s=0, then it is the next block.
    # Value of array os the previous value + current line
    [[ $s == 0 ]] && arr[$i]="${arr[$i]}
$line" || {
            # Otherwise the value of array is current line & s=0
            arr[$i]="$line"
            s=0;
    }
done < tmp_LG.csv

# Reading in corresponding CG atoms
conh2=$(echo "ligand_mapping_CG_files/conh2_ASN.txt" "ligand_mapping_CG_files/conh2_GLN.txt")
bbcco=$(echo "ligand_mapping_CG_files/bbcco_GLY.txt" "ligand_mapping_CG_files/bbcco_ALA.txt" "ligand_mapping_CG_files/bbcco_PRO.txt")
ph1=$(echo "ligand_mapping_CG_files/ph_PHE_1.txt") # 6 rotations
ph2=$(echo "ligand_mapping_CG_files/ph_PHE_2.txt")
ph3=$(echo "ligand_mapping_CG_files/ph_PHE_3.txt")
ph4=$(echo "ligand_mapping_CG_files/ph_PHE_4.txt")
ph5=$(echo "ligand_mapping_CG_files/ph_PHE_5.txt")
ph6=$(echo "ligand_mapping_CG_files/ph_PHE_6.txt")
bb_cnh=$(echo "ligand_mapping_CG_files/bb_cnh_GLY.txt" "ligand_mapping_CG_files/bb_cnh_LYS.txt" "ligand_mapping_CG_files/bb_cnh_ALA.txt")
ccn=$(echo "ligand_mapping_CG_files/ccn_LYS.txt")
ccoh=$(echo "ligand_mapping_CG_files/ccoh_SER.txt" "ligand_mapping_CG_files/ccoh_THR.txt")
coh=$(echo "ligand_mapping_CG_files/coh_SER.txt" "ligand_mapping_CG_files/coh_THR.txt")
coo=$(echo "ligand_mapping_CG_files/coo_ASP.txt" "ligand_mapping_CG_files/coo_GLU.txt")
csc=$(echo "ligand_mapping_CG_files/csc_MET.txt")
csh=$(echo "ligand_mapping_CG_files/csh_CYS.txt")
gn=$(echo "ligand_mapping_CG_files/gn_ARG.txt")
hid=$(echo "ligand_mapping_CG_files/hid_HIS.txt")
hie=$(echo "ligand_mapping_CG_files/hie_HIS.txt")
hip=$(echo "ligand_mapping_CG_files/hip_HIS.txt")
indole=$(echo "ligand_mapping_CG_files/indole_TRP.txt")
phenol1=$(echo "ligand_mapping_CG_files/phenol_TYR_1.txt")
phenol2=$(echo "ligand_mapping_CG_files/phenol_TYR_2.txt") # phenol numbering in opposite direction
isopropyl1=$(echo "ligand_mapping_CG_files/isopropyl_LEU_1.txt" "ligand_mapping_CG_files/isopropyl_VAL_1.txt")
isopropyl2=$(echo "ligand_mapping_CG_files/isopropyl_LEU_2.txt" "ligand_mapping_CG_files/isopropyl_VAL_2.txt")  # isopropyl numbering in opposite direction
pro1=$(echo "ligand_mapping_CG_files/pro_PRO_1.txt")
pro2=$(echo "ligand_mapping_CG_files/pro_PRO_2.txt") # proline numbering in opposite direction
ch3=$(echo "ligand_mapping_CG_files/ch3_ALA.txt" "ligand_mapping_CG_files/ch3_ILE.txt")


# Loop through the array and run the ligand support script for corresponding CG atoms
for i in "${arr[@]}"
do
# bb_cco groups
   if [[ "$i" == *"bb_cco"* ]]; then
      echo "################ C -> CA "
      for a in $bbcco; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'blank'
      done
 # alternate bb_cco groups
      echo "################ alternate N -> CA "
      for a in $bbcco; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'CA' | sed 's/blank/CA/g'
      done

# conh2 groups
   elif [[ "$i" == *"conh2"* ]]; then
      echo "################"
      for a in $conh2; do
         printf '%s\n' $i | paste -d',' - $a
      done

# ph groups - 6 rotations of symmetry -- each one also has a mirror image
elif [[ "$i" == *"ph"* ]] && [[ "$i" != *"phenol"* ]]; then
      echo "################"
      for a1 in $ph1; do
         printf '%s\n' $i | paste -d',' - $a1
         echo "################"
         printf '%s\n' $i | tail -r | paste -d',' - $a1
      done
      echo "################"
      for a2 in $ph2; do
         printf '%s\n' $i | paste -d',' - $a2
         echo "################"
         printf '%s\n' $i | tail -r | paste -d',' - $a2
      done
      echo "################"
      for a3 in $ph3; do
         printf '%s\n' $i | paste -d',' - $a3
         echo "################"
         printf '%s\n' $i | tail -r | paste -d',' - $a3
      done
      echo "################"
      for a4 in $ph4; do
         printf '%s\n' $i | paste -d',' - $a4
         echo "################"
         printf '%s\n' $i | tail -r | paste -d',' - $a4
      done
      echo "################"
      for a5 in $ph5; do
         printf '%s\n' $i | paste -d',' - $a5
         echo "################"
         printf '%s\n' $i | tail -r | paste -d',' - $a5
      done
      echo "################"
      for a6 in $ph6; do
         printf '%s\n' $i | paste -d',' - $a6
         echo "################"
         printf '%s\n' $i | tail -r | paste -d',' - $a6
      done

# bb_cnh groups
   elif [[ "$i" == *"bb_cnh"* ]]; then
      echo "################"
      for a in $bb_cnh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# ccn groups
   elif [[ "$i" == *"ccn"* ]]; then
      echo "################"
      for a in $ccn; do
         printf '%s\n' $i | paste -d',' - $a
      done

# ccoh groups
   elif [[ "$i" == *"ccoh"* ]]; then
      echo "################"
      for a in $ccoh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# coh groups
   elif [[ "$i" == *"coh"* ]]; then
      echo "################"
      for a in $coh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# coo groups
   elif [[ "$i" == *"coo"* ]]; then
      echo "################"
      for a in $coo; do
         printf '%s\n' $i | paste -d',' - $a
      done

# csc groups
   elif [[ "$i" == *"csc"* ]]; then
      echo "################"
      for a in $csc; do
         printf '%s\n' $i | paste -d',' - $a
      done

# csh groups
   elif [[ "$i" == *"csh"* ]]; then
      echo "################"
      for a in $csh; do
         printf '%s\n' $i | paste -d',' - $a
      done

# gn groups
   elif [[ "$i" == *"gn"* ]]; then
      echo "################"
      for a in $gn; do
         printf '%s\n' $i | paste -d',' - $a
      done

# hip groups
   elif [[ "$i" == *"hip"* ]]; then
      echo "################"
      for a in $hip; do
         printf '%s\n' $i | paste -d',' - $a
      done

# hie and hid groups
   elif [[ "$i" == *"his"* ]]; then
      echo "################"
      for a in $hie; do
         printf '%s\n' $i | paste -d',' - $a
      done
      echo "################"
      for a in $hid; do
         printf '%s\n' $i | paste -d',' - $a
      done

# indole groups
   elif [[ "$i" == *"indole"* ]]; then
      echo "################"
      for a in $indole; do
         printf '%s\n' $i | paste -d',' - $a
      done

# phenol groups - 6 rotations
   elif [[ "$i" == *"phenol"* ]]; then
      echo "################"
      for a in $phenol1; do
         printf '%s\n' $i | paste -d',' - $a
      done
      echo "################"
      for a in $phenol2; do
         printf '%s\n' $i | paste -d',' - $a
      done

# isopropyl groups
   elif [[ "$i" == *"isopropyl"* ]]; then
      echo "################"
      for a in $isopropyl1; do
         printf '%s\n' $i | paste -d',' - $a
      done

# isopropyl groups - opposite direction
  # elif [[ "$i" == *"isopropyl"* ]]; then
      echo "################ - opposite direction"
      for a in $isopropyl2; do
         printf '%s\n' $i | paste -d',' - $a
      done

# pro groups
   elif [[ "$i" == *"pro"* ]]; then
      echo "################"
      for a in $pro1; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'blank'
      done

# pro groups - opposite direction
   #elif [[ "$i" == *"pro"* ]]; then
      echo "################ - opposite direction"
      for a in $pro2; do
         printf '%s\n' $i | paste -d',' - $a | grep -v 'blank'
      done

# ch3 groups
   elif [[ "$i" == *"ch3"* ]]; then
      echo "################"
      for a in $ch3; do
         printf '%s\n' $i | paste -d',' - $a
      done

# default
   else
      echo "################"
      echo "PARSING ERROR: no combs groups in CG coordinate file. Re-run ligand_matcher.py"
   fi

done > tmp_LG_cg_coordinates.csv

# turn this output csv into a second array to work with CG in 'chunks' aka multiple arrays
p=1
k=1
declare -a arr2
while read -r line
do
    [[ $line == *"################"* ]] && ((p++)) && k=1 && continue
    [[ $k == 0 ]] && arr2[$p]="${arr2[$p]}
$line" || {
            arr2[$p]="$line"
            k=0;
    }
done < tmp_LG_cg_coordinates.csv

echo " \n" > tmp_LG8.csv

#ligand_code=$(grep -A1 'ATOM' $l | grep -v 'ATOM' | awk -F' ' '{print $8}' | uniq) # for mol2
ligand_code=$(grep -e 'ATOM' -e 'HETATM' $l | awk -F' ' '{print $4}' | uniq) # for PDB

COUNTER=0
for combs_groups_array in "${arr2[@]}"; do
   CG_block=$(printf '%s\n' $combs_groups_array)
   let COUNTER=COUNTER+1
   for CG_line in $CG_block; do
      x_coord=$(echo $CG_line | awk -F',' '{print $1}') # X-coord from input coords file
      y_coord=$(echo $CG_line | awk -F',' '{print $2}') # Y-coord from input coords file
      z_coord=$(echo $CG_line | awk -F',' '{print $3}') # Z-coord from input coords file

      # If all three coordinates match your input PDB file, then print out the corresponding atom name
      #if [ -z "$x_coord" ] || [ -z "$y_coord" ] || [ -z "$z_coord" ]; then
        pdb_atom=$(grep -e $x_coord $l | grep -e $y_coord | awk -F' ' '{print $3}') # $2 for mol2, $3 for pdb
        # took out the grep -e $z_coord 

      # Print out the other info that the script has already found --
      other_info=$(echo $CG_line | awk -F',' '{print $6,$5,$4}')

      # if the coordinate is in the ligand file, then add it to the output file
      # this removes the excess hydrogens atoms not in original file  -- there shouldn't be excess from RDKit though
      if grep -Fqe $x_coord $l
      then
         echo $ligand_code $pdb_atom $other_info $COUNTER
      fi
   done > tmp_LG8_unsorted.csv
   cat tmp_LG8_unsorted.csv | sort -k3,4
   echo ""

done > tmp_LG8.csv

# Print out all of the atoms that correspond to the unique group number.
# This will be used to find the ligand coverage number.
# The next line of code prints out all of the atoms that correspond to the group number.
# After that, we will loop through them to see if there are 2+ atoms that match between groups.
awk '{print $6,$2}' tmp_LG8.csv | sort | uniq | awk '{if(a!=$1) {a=$1; printf "\n%s%s",$0,FS} else {a=$1;$1="";printf $0 }} END {printf "\n"}' | sed '/^$/d' | sed 's/ /,/' | sed 's/ //' > tmp_LG9.csv

# This starts with the large ligand groups first (the order from the python script)
# It assigns the unique atoms groups from largest to smallest, and as it reads through
# the atom list file from the prior step, it will assign a CG coverage number to each CGs.
# If it then comes across a new (smaller) group with 2 or more matching atoms,
# then it will assign it to the same CG coverage number
awk -F, '
{ # cycle 1 - assign a CG coverage number
  CGs[ $1 ] =  $2
  split( $2, atoms, " " )
  for (atom in atoms) {
    atoms_in_CGs[ atoms[atom] ] = atoms_in_CGs[ atoms[atom] ] " " $1
  }
}
END {  # cycle 2 - if there are >=2 atoms that match another CG coverage number,
        # then assign CG coverage number to the first CG coverage number
  for ( CG in CGs ) {
    split( CGs[CG], atoms, " " )
    for ( atom in atoms ) {
      split( atoms_in_CGs[ atoms[atom] ] , in_CGs, " " )
      if (length( in_CGs ) >= 2)  {
        print CG,CGs[CG],in_CGs[1]
        break
      }
    }
  }
}
' tmp_LG9.csv | sort -n > tmp_LG10.csv


# Creating the ligand.txt files

if [ -z "$o" ]; then # If your output filename already exists,
   if [ -f "ligand.txt" ]; then
      if [ -d "./old_files" ]; then   # then this will move it to "old_files/"
         echo "Ligand.txt already exists. Moving old file to old_files directory."
         mv ligand.txt ./old_files/
      else
         echo "Ligand.txt already exists. Moving old file to old_files directory."
         mkdir ./old_files  # make the old_files directory if it is not already there
         mv ligand.txt ./old_files/
      fi
   fi
   while read -r line; do
      group_num=$(echo $line | awk -F' ' '{print $1}')
      cluster_num=$(echo $line | awk -F' ' '{print $NF}')
      awk -v a="$group_num" -v b="$cluster_num" 'a == $6 {print $0,b}' tmp_LG8.csv
      echo ""
   done < tmp_LG10.csv > ligand.txt
else
  if [ -f "$o" ]; then
     if [ -d "./old_files" ]; then
        echo  $o "already exists. Moving old file to old_files directory."
        mv ligand.txt ./old_files/
     else
        echo $o "already exists. Moving old file to old_files directory."
        mkdir ./old_files
        mv ligand.txt ./old_files/
     fi
  fi
   while read -r line; do
      group_num=$(echo $line | awk -F' ' '{print $1}')
      cluster_num=$(echo $line | awk -F' ' '{print $NF}')
      awk -v a="$group_num" -v b="$cluster_num" 'a == $6 {print $0,b}' tmp_LG8.csv
      echo ""
   done < tmp_LG10.csv > $o
fi

# Removes the temporary files from this script that were created
rm tmp_LG*
