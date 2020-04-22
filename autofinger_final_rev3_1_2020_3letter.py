#import proper modules
import re
import pandas as pd
import csv
import os
import glob
#find and iterate all pdb files in the directory
pdbs = glob.glob('*.pdb')
i=0
for r in pdbs:
    ref_pdb = pdbs[i]
    i = i+1
    residue_key = {}
    with open(ref_pdb, 'r') as parse:
        value = "0"
        for line in parse:
            if line[:4] == "ATOM":
                #check to see if the residue number has increased, adds residues to dictonary
                if line[20:26] != value:
                    value = line[22:26]
                    residue = line[17:20]
                    residue_key[value.strip()] = residue
    #states the name of the csv file as the PDBname_finger.csv
    sourcefile = ref_pdb[:-4] + "_finger.csv"
    #transposes the csv file for easier parsing
    pd.read_csv(sourcefile, header=None).T.to_csv('finger_flip.csv', header=False, index=False)
    #parses through the transposed file for interactions
    with open('finger_flip.csv', 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            next(csv_reader)
            contact = []
            backbone = []
            sidechain = []
            polar = []
            hydrophobic = []
            acceptor = []
            donor = []
            aromatic = []
            charged = []
            for line in csv_reader:
                key = line[0]
                matchkey = ["charged", "contact", "backbone", "acceptor", "aromatic", "polar", "donor", "hydrophobic", "sidechain"]
                if key[-7:] in matchkey or key[-8:] in matchkey or key[-5:] in matchkey or key[-11:] in matchkey or key[-9:] in matchkey:
                    if line[1] == "1":
                        interxn= line[0]
                        _index = interxn.find("_")
                        resnum_raw = str(interxn[:_index]).strip()
                        resnum_edit = re.sub(r'[A-Z]+', '', resnum_raw, re.I)
                        intype = interxn[(_index+1):]
                        eval(intype).append(str(residue_key.get(resnum_edit)) + str(resnum_edit)+ " ")
                    else:
                        continue
                else:
                    continue
            final_dict = {"Backbone":backbone,"Polar":polar,"Hydrophobic":hydrophobic,"Hydrogen Bond Acceptor":acceptor,"Hydrogen Bond Donor":donor,"Aromatic":aromatic,"Charged":charged}
    str_dict = {"Backbone":"".join(backbone),"Polar":"".join(polar),"Hydrophobic":"".join(hydrophobic),"Hydrogen Bond Acceptor":"".join(acceptor),"Hydrogen Bond Donor":"".join(donor),"Aromatic":"".join(aromatic),"Charged":"".join(charged)}
    exp_file = ref_pdb[:-4] + "_table.csv"   
    with open(exp_file, "w", newline = '') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(str_dict.keys())
        writer.writerows([str_dict.values()])
#removes the temporary transposed file
os.remove("finger_flip.csv")