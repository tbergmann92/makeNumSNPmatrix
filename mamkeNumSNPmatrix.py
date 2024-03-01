#!/usr/bin/python3

## work in progress !!!

## Script to convert SNP chip data into a numeric matrix (0,1,2,3) based on the MAF.
## The coding of the output matrix can be used as input for SNPRelate.
## 2 --> AA
## 1 --> AB
## 0 --> BB
## 3 --> -

import openpyxl
import pandas as pd
import numpy as np

Matrix = pd.read_excel("./Faba_Test/Faba_10k.xlsx", sheet_name="prime", index_col=0)
Matrix = Matrix.replace(["failed"], "-")

## prepare new data frame
col_names = Matrix.columns.values.tolist()
col_names.insert(0, "Marker")
num_df = pd.DataFrame(columns=col_names)

## iterate through rows
marker_stats = []
char_to_replace = {k:"" for k in "RYSWKMH-"}
for _, row in Matrix.iterrows():
    # define variables
    num_alleles, num_allele_list = [], []
    all_allele_counts, main_allele_counts = {}, {}
    # retrieve all snp calls for the marker
    all_snp_calls = list(row.values)
    # subset the row values only to the main SNP calls 
    main_snp_calls = list("".join(all_snp_calls).translate(str.maketrans(char_to_replace)))
    # get stats for each marker for the main alleles
    for main_allele in main_snp_calls:
        main_allele_counts[main_allele] = main_allele_counts.get(main_allele, 0) + 1
        # get stats for each marker for the main alleles
    #for all_calls in all_snp_calls:
    #    all_allele_counts[all_calls] = all_allele_counts.get(all_calls, 0) + 1
    #print(_, main_allele_counts)
    #print(_, main_allele_counts)
    # convert alleles into numeric values (use the subset list for this!)
    for value in row.values[:]:
        if value not in ("R","Y","S","W","K","M","-") and (value == max(main_allele_counts, key=main_allele_counts.get)):
            value = "2"
            num_allele_list.append(value)
        elif value not in ("R","Y","S","W","K","M","-"):
            value = "0"
            num_allele_list.append(value)
        elif value == "-":
            value = "3"
            num_allele_list.append(value)
        else:
            value = "1"
            num_allele_list.append(value)
    print(_, main_allele_counts, num_allele_list)
    num_alleles = ",".join(num_allele_list)
    new_num_row = (f"{_},{num_alleles}")
    num_df.loc[len(num_df)] = new_num_row.split(",")

df.to_excel("./Genotypes_0123_Matrix.xlsx", index=False)
