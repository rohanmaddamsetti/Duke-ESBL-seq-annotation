#!/usr/bin/env python

"""
tabulate-ESBL-proteins.py by Rohan Maddamsetti

Usage: python tabulate-ESBL-proteins.py > ../results/annotated_ESBL_proteins.csv

"""

import os
import re
from Bio import SeqIO


def get_prot_data(feature):
    try:
        prot_id = feature.qualifiers['protein_id'][0]
    except:
        prot_id = "NA"
    
    try:
        prot_seq = feature.qualifiers['translation'][0]
    except:
        prot_seq = "NA"
        
    try: ## replace all commas with semicolons for csv formatting.
        prot_product = feature.qualifiers['product'][0].replace(',',';')
    except:
        prot_product = "NA"

    try: ## replace all commas with semicolons for csv formatting.
        prot_note = feature.qualifiers['note'][0].replace(',',';')
    except:
        prot_note = "NA"
        
    prot_location = str(feature.location)
                        
    cur_prot = {
        "id" : prot_id,
        "seq" : prot_seq,
        "product" : prot_product,
        "note" : prot_note,
        "location" : prot_location }
    return cur_prot

    
def is_beta_lactamase_or_regulator(cur_prot):
    '''
    return False if the protein is not a beta-lactamase or regulator
    return True if the protein is a beta-lactamase or regulator
    '''
    beta_lactam_keywords = "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\S*"
    ## check if the gene is a beta-lactamase or a beta-lactamase regulator.
    is_betalactam = True if re.search(beta_lactam_keywords, cur_prot["product"]) or re.search(beta_lactam_keywords, cur_prot["note"]) else False
    if not is_betalactam:
        return False
    else:
        return True


def get_strains_of_interest(myfilepath="../data/strains-of-interest.csv"):
    strains_of_interest = list()
    with open(myfilepath, "r") as strains_fh:
        for i, line in enumerate(strains_fh):
            ## remove leading and lagging whitespace like "\n" at the ends of lines.
            line = line.strip() 
            if i == 0: continue ## skip the header
            prefix, suffix = line.split(",")
            my_strain = prefix + suffix
            strains_of_interest.append(my_strain)
    return strains_of_interest


def main():
    datadir = "../data"

    ## import the list of strains to consider.
    strains_of_interest = get_strains_of_interest("../data/strains-of-interest.csv")
    strains_of_interest_with_genomes = list() ## keep track of the strains we have genomes for.
    
    header = "NCBI_BioProject,AnnotationAccession,Organism,Strain,NCBI_ProteinID,ProteinProduct,ProteinNote,ProteinSequence"
    print(header)
    
    ## check for 'PRJ' prefix to ignore files like .DS_Store.
    for NCBI_BioProject in [x for x in os.listdir(datadir) if x.startswith("PRJ")]:
        bioproject_dir = os.path.join(datadir, NCBI_BioProject)
        ## iterate through the genomes in the BioProject
        for genome_id in [x for x in os.listdir(bioproject_dir) if x.startswith("GCA") or x.startswith("GCF")]:
            ## keep track of whether this genome is a strain of interest or not
            is_strain_of_interest = True ## default value is true
            ## now open the actual Genbank annotation file
            my_gbk_path = os.path.join(bioproject_dir, genome_id, "genomic.gbff")
            with open(my_gbk_path,'rt') as gbk_fh:
                for replicon in SeqIO.parse(gbk_fh, "gb"):
                    if not is_strain_of_interest:
                        break
                    else:
                        for feature in replicon.features:
                            if feature.type == "source":
                                try:
                                    my_strain = feature.qualifiers['strain'][0]
                                    if my_strain not in strains_of_interest:
                                        is_strain_of_interest = False
                                        break
                                    else:
                                        strains_of_interest_with_genomes.append(my_strain)
                                except:
                                    my_strain = "NA"
                                    is_strain_of_interest = False
                                    break
                                try:
                                    my_species = feature.qualifiers['organism'][0]
                                except:
                                    my_species = "NA"
                            elif feature.type != "CDS":
                                continue ## skip over everything that is not a protein-coding gene.
                            else: ## the feature.type == "CDS"
                                cur_prot = get_prot_data(feature)
                                if is_beta_lactamase_or_regulator(cur_prot) and is_strain_of_interest:
                                    my_row = ",".join([NCBI_BioProject, genome_id, my_species, my_strain, cur_prot["id"], cur_prot["product"], cur_prot["note"], cur_prot["seq"]])
                                    print(my_row)

    ## print out the set of strains of interest without genomes
    strains_of_interest_without_genomes = set(strains_of_interest) - set(strains_of_interest_with_genomes)
    with open("../results/strains-of-interest-without-genomes.txt", "w") as outfh:
        for s in strains_of_interest_without_genomes:
            outfh.write(s + "\n")
    
    return

                            
## run the script.
if __name__ == "__main__":
    main()

    
