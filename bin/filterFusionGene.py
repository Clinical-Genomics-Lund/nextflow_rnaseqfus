#!/usr/bin/env python

import re 
import os 
import argparse

def CreateDictionary (fusionPanel: str) -> list:
    """
    Creates a dictonary of fusion panel of genes such that anchor genes and paired genes 
    """
    FilePath = os.path.abspath(fusionPanel)
    try:
        fusion_dict = {}
        with open(FilePath) as f:
            for line in f:
                fusionline = line.strip().split("\t")
                keys = str(fusionline[0])
                values = fusion_dict.get(keys,'')
                values = fusionline[1]
                fusion_dict[keys] = values
        return fusion_dict
    except FileNotFoundError:
        print (f"Sorry, the file does not exit")  

def PrintFusion(gene1,gene2,annotation1,annotation2,line,fusion_dict,file,confidence):
    """
    Checks for the genes pairs and annotations which are not in intronic regions. Prints selected line in new file.
    """
    try:
        if confidence == 'high' and annotation1 != "intron" and annotation2 != "intron":
            print (gene1,gene2)
            file.write(line)

        if confidence == 'medium' and annotation1 != "intron" and annotation2 != "intron":
            print (gene1,gene2)
            file.write(line)

        if gene1 in fusion_dict.keys() and confidence == 'low':
            if gene2 in fusion_dict[gene1] and gene1 != gene2 and annotation1 != "intron" and annotation2 != "intron":
                print (gene1,gene2)
                file.write(line)

        if gene1 in fusion_dict.values() and confidence == 'low':
            x = [i for i, d in enumerate(fusion_dict.values()) if gene1 in d]
            for index in x:
                y =  [d for i, d in enumerate(fusion_dict.keys()) if index==i]
                gene = ("".join(y))
                if gene == gene2 and gene1 != gene2 and annotation1 != "intron" and annotation2 != "intron":
                    print (gene1,gene2)
                    file.write(line)
        return(file)
    except NameError:
        print (f"Sorry, the variables doesn't exit")  

def SelectFusion(fusion_dict, AllFusion):
    """
    Create the selected fusion list based on the gene list provided from dictionary
    """
    FilePath = os.path.abspath (AllFusion)
    try:
        file = open("Selected.fusion.tsv", "w")
        with open(FilePath) as fusionFile:
            next(fusionFile)
            for line in fusionFile:
                fusionline = line.strip().split("\t")
                gene1 = fusionline[0]
                gene2 = fusionline[1]
                annotation1 = fusionline[6]
                annotation2 = fusionline[7]

                confidence = fusionline[14]
                PrintFusion (gene1, gene2, annotation1,annotation2,line,fusion_dict,file,confidence)
                PrintFusion(gene2, gene1, annotation1,annotation2,line,fusion_dict,file,confidence)
        file.close()
        return file
    except FileNotFoundError:
            print (f"Sorry, the file does not exit")

def Main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--g',
                    dest='genelist',
                    default='fusion.panel.AL',
                    help='Fusion panel defined for gene selection'
                    )

    parser.add_argument('--f',
                    dest='fusion',
                    default='None',
                    help='Accepted and discarded gene list from arriba'
                    )

    args = parser.parse_args()
    fusion = args.genelist
    allFusion = args.fusion
    fusionDict = CreateDictionary(fusion)
    SelectFusion(fusionDict, allFusion)

if __name__ == "__main__":
    Main()


