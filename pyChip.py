#!/usr/bin/python3

## Script to convert SNP chip data into a numeric matrix (0,1,2,3) based on the MAF.
## The coding of the output matrix can be used as input for SNPRelate.
## 2 --> AA (Major allele)
## 1 --> AB
## 0 --> BB (Minor allele)
## 3 --> -

import os
import sys
import openpyxl
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt  
import plotly.express as px

# Argument Parsing
parser = argparse.ArgumentParser(description = "Script for analyzing and converting SNP-chip data.",
                                 epilog="Note: input matrix should be in (n x m) format with markers as n-rows and genotypes as m-columns.",
                                 usage = "python3 pyChip.py -i ./path/to/input.xlsx -o ./path/to/outputDir -r ./path/to/color_coding_file",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", help = "SNP chip data to be analyzed (first excel sheet will be used!)", required=True)
parser.add_argument("-o", "--outputDir", help = "Define the output folder for results", required=True) 
parser.add_argument("-r", "--color_coding_file_pca", help = "Define the groups and the color coding for the pca plot", required=False) 

# Add validation checks for required arguments
args = parser.parse_args()

# Check input parameters
if not all ([args.input, args.outputDir]):
    sys.exit("\033[1;31;40m ERROR! All of the following arguments are required: -i/--input, -o/--outputDir")

# Input file validation function
def check_file_existence(file_path):
    if not os.path.isfile(file_path):
        sys.exit(f"{file_path} doesn't exist! Please check your input files.")

# Check whether input file exists
check_file_existence(args.input)

# create the matrix object
class DataMatrix:
    # 
    def __init__(self, input_file, sheet, output_dir):
        self.input_file = input_file
        self.sheet = sheet
        self.output_dir = output_dir
        self.raw_data = None
        self.marker_statistics = None
        self.marker_status = None
        self.numeric_matrix = None
        self.kinship_matrix = None
        self.maf_data = None
        self.hapmap = None
    #
    def create_output_directory(self, overwrite=True):
        """ Checks whether output already exists."""
        if os.path.exists(self.output_dir):
            if overwrite:
                os.makedirs(self.output_dir, exist_ok=True)
                return self.output_dir
            else:
                sys.exit("\033[1;31;40m Attention! Output directory already exists!")
        else:
            os.makedirs(self.output_dir)
            return self.output_dir
    #
    def load_data(self):
        """ Loads input and converts "failed" to "-". """
        try:
            self.raw_data = pd.read_excel(self.input_file, sheet_name=self.sheet, index_col=0)
            self.raw_data = self.raw_data.replace(["failed"],"-")
            #print("Loaded data succesfully!")
        except FileNotFoundError:
            sys.exit("Error: Input file not found! Please check your file path.")
        except pd.errors.EmptyDataError:
            sys.exit("Error: Input file is empty!")
    #
    def data_check(self):
        """ Checks whether unrecognized SNP calls are in the matrix. """
        unique_characters = list(np.unique(self.raw_data.values))
        allowed_characters = ["-","A","C","G","T","K","M","R","S","T","W","Y"]
        if any(char not in allowed_characters for char in unique_characters):
            print("\033[1;31;40m ATTENTION: Unknown characters identified in the matrix!")
            print("\033[1;37;40m Identified characters: " + str(unique_characters))
            sys.exit("Allowed characters: " + str(allowed_characters))
        else:
            pass
    #
    def calc_snp_stats(self):
        """ Calculates statistics about the SNP calls across the matrix. """
        snp_a, snp_t, snp_c, snp_g = (self.raw_data.values == "A").sum(), (self.raw_data.values == "T").sum(), (self.raw_data.values == "C").sum(), (self.raw_data.values == "G").sum()
        snp_r, snp_y, snp_s, snp_w, snp_k, snp_m  = (self.raw_data.values == "R").sum(), (self.raw_data.values == "Y").sum(), (self.raw_data.values == "S").sum(), (self.raw_data.values == "W").sum(), (self.raw_data.values == "K").sum(), (self.raw_data.values == "M").sum()
        snp_fails = (self.raw_data.values == "-").sum()
        sum_hets = snp_r + snp_y + snp_s + snp_w + snp_k + snp_m
        sum_snps = snp_a + snp_t + snp_c + snp_g + sum_hets + snp_fails
        # create dataframe of SNP stats
        snp_stats = pd.DataFrame({
            "Type":["A", "T", "C", "G", "HETs", "Failed", "Total"],
            "Count":[snp_a, snp_t, snp_c, snp_g, sum_hets, snp_fails, sum_snps]
        })
        snp_stats["Fraction"] = (snp_stats["Count"]/snp_stats["Count"].values[-1])*100
        # add to object
        self.snp_stats = snp_stats
        # define the output for the picture
        out_pict1 = os.path.join(self.output_dir, "SNPStats.png").replace("\\", "/")
        # export SNP stats
        self.make_barplots(data = snp_stats.iloc[:-1], 
            output = out_pict1, 
            title = "SNP calls",
            bar_color = "#6799a8", # Optional: Specifiy the color of the bars
            text_color = "#120f0f" # Optional: Specifiy the color of the text labels
            )
    #
    def calc_marker_status(self):
        """ Calculates whether a marker is polymorphic or not.
        
        Returns:
        
        pol_marker_fraction (float): Percentage of polymorphic markers.
        mon_marker_fraction (float): Percentage of monomorphic markers.
        fail_marker_fraction (float): Percentage of markers with no alleles.
        
        """
        # define variables
        char_to_replace = str.maketrans("","","RYSWKMH-")
        allele_a, allele_b = [], []
        mon_counter, fail_counter, pol_counter = 0,0,0
        #
        for _, row in self.raw_data.iterrows():
            counter = 0
            all_snps = list(row.values)
            main_snps = list("".join(all_snps).translate(char_to_replace))
            unique_snps = list(dict.fromkeys(main_snps)) 
            if len(unique_snps) == 2:
                allele_a.append(unique_snps[0])
                allele_b.append(unique_snps[1])
                pol_counter += 1
            ## if only one alleles is present construct the list as following:
            elif len(unique_snps) == 1:
                allele_a.append(unique_snps[0])
                allele_b.append("NA")
                mon_counter += 1
            ## if no allele is present construct the list as following:
            else:
                allele_a.append("NA")
                allele_b.append("NA")
                fail_counter += 1
        #
        total_markers = self.raw_data.shape[0]
        pol_marker_fraction = self.calculate_percentage(pol_counter, total_markers)
        mon_marker_fraction = self.calculate_percentage(mon_counter, total_markers)
        fail_marker_fraction = self.calculate_percentage(fail_counter, total_markers)
        #
        marker_status = pd.DataFrame({
            "Type":["Polymorph", "Monomorph", "Failed"],
            "Fraction":[pol_marker_fraction, mon_marker_fraction, fail_marker_fraction]
        })
        # add to object
        self.marker_status = marker_status
        # define the output for the picture
        out_pict2 = os.path.join(self.output_dir, "MarkerStats.png").replace("\\", "/")
        # export marker stats
        self.make_barplots(data = marker_status, 
                      output = out_pict2, 
                      title = "Marker calls",
                      bar_color = "#6799a8", # Optional: Specifiy the color of the bars
                      text_color = "#120f0f" # Optional: Specifiy the color of the text labels
                      )
    #
    def calculate_percentage(self, count, total):
        """ Calculate percentage given count and total. """
        return (count / total) * 100
    #
    def align_center(self, x):
        """ Function to center the cell values for export to excel."""
        return ['text-align: center' for x in x]
    #
    def make_barplots(self, data, output, title, bar_color, text_color):
        """ Takes dataframe with SNP counts, the name of the output file with png extension, and the title of the plot."""
        # figure size
        plt.rcParams['figure.figsize'] = (5, 6)
        fig, ax = plt.subplots()
        # create barplots
        bars = ax.bar(x = data["Type"],
            height = data["Fraction"],
            color=bar_color,
            edgecolor="black"
        )
        # plot the figures over the bars
        for bar in bars:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 1,
                round(bar.get_height(), 1),
                horizontalalignment='center',
                color=text_color,
                weight='bold'
            )
        # adjust style
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='#9c9a9a')
        ax.xaxis.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_color('#DDDDDD')
        ax.tick_params(bottom=False, left=False)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='#ffffff')
        ax.xaxis.grid(True, color='#ffffff')
        ax.set_ylabel("Fraction [%]", labelpad=15, color="#080404", weight="bold")
        ax.set_title(title, pad=15, color="#080404", weight="bold")
        ax.set_facecolor("#ebe8e8")
        ax.set_ylim(0,100)
        ## save figure
        fig.tight_layout()
        plt.savefig(output, bbox_inches="tight", dpi=300)        
    #
    def generate_numeric_matrix(self):
        """ Converts the raw SNP matrix into numeric matrix. 
        
        Coding:
        
        2 --> Major allele.
        1 --> Heterozygous allele.
        0 --> Minor allele.
        3 --> Failed SNP call.
        
        """
        ## prepare new data frame
        col_names = self.raw_data.columns.values.tolist()
        col_names.insert(0, "Marker")
        num_df = pd.DataFrame(columns=col_names)
        ## iterate through rows
        new_cols, row_id = [], []
        keys_to_keep = ["A","G","C","T"]
        #
        for idx, row in self.raw_data.iterrows():
            # define variables
            num_alleles, num_allele_list = [], []
            all_allele_counts = {}
            # retrieve all snp calls for the marker
            snp_calls = list(row.values)
            # get stats for each marker for allele calls
            for call in snp_calls:
                all_allele_counts[call] = all_allele_counts.get(call, 0) + 1
            # filter the dictionary for the main allele calls
            filtered_all_allele_counts = {key: all_allele_counts[key] for key in keys_to_keep if key in all_allele_counts}        
            # convert alleles into numeric values (use the subset list for this!)
            if not bool(filtered_all_allele_counts):
                for value in row.values[:]:
                    value = "3"
                    num_allele_list.append(value)
            if bool(filtered_all_allele_counts):    
                for value in row.values[:]:
                    if value not in ("R","Y","S","W","K","M","-") and (value == max(filtered_all_allele_counts, key=filtered_all_allele_counts.get)):
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
            row_id.append(idx)
            new_cols.append(all_allele_counts)
            num_alleles = ",".join(num_allele_list)
            new_num_row = (f"{idx},{num_alleles}")
            num_df.loc[len(num_df)] = new_num_row.split(",")    
        marker_statistics = pd.DataFrame(new_cols).fillna(0)
        marker_statistics.insert(0, "Marker", row_id)
        num_df_col = num_df.loc[:,"Marker"].tolist()
        num_df = num_df.drop(num_df.columns[0], axis=1)
        num_df.index = num_df_col
        num_df = num_df.astype(int)
        # add to object 
        self.numeric_matrix = num_df
        self.marker_statistics = marker_statistics
    #
    def calculate_pca(self):
        # retrieve matrix from object
        numeric_matrix_transposed = self.numeric_matrix.T
        # filtering and preparing the data
        filter_values = [3,2,1,0]
        for value in filter_values:
            numeric_matrix_transposed = numeric_matrix_transposed.loc[:,(numeric_matrix_transposed != value).any()]
        numeric_matrix_transposed["Marker"] = numeric_matrix_transposed.index
        numeric_matrix_transposed = numeric_matrix_transposed.reset_index(drop=True)
        extracted_cols = numeric_matrix_transposed.columns.values.tolist()
        # extract the array
        x = numeric_matrix_transposed.loc[:,extracted_cols[:-1]].values
        y = numeric_matrix_transposed.loc[:,["Marker"]].values
        # perform the standardization
        x = StandardScaler().fit_transform(x)
        # compute the PCA
        pca = PCA(n_components=3)
        principal_components = pca.fit_transform(x)
        # prepare the results as df
        principal_data_frame = pd.DataFrame(data = principal_components,
                                            columns = ["PC1","PC2","PC3"])
        PCA_data_frame_final = pd.concat([principal_data_frame, numeric_matrix_transposed[["Marker"]]], axis=1)
        PCA_data_frame_final = PCA_data_frame_final.rename(columns={"Marker":"Genotype"})
        # extract amount of variance
        explained_variance = pca.explained_variance_ratio_
        # set max axis scales
        x_range = [PCA_data_frame_final['PC1'].min() - 5, PCA_data_frame_final['PC1'].max() + 5]
        y_range = [PCA_data_frame_final['PC2'].min() - 5, PCA_data_frame_final['PC2'].max() + 5]
        # Load genotype-group mappings from Excel file
        genotype_group_df = pd.read_excel(args.color_coding_file_pca, sheet_name="Group_Coding", dtype=object)
        genotype_group_dict = dict(zip(genotype_group_df["Genotype"], genotype_group_df["Group"]))
        # Load group-color mappings from Excel file
        group_color_df = pd.read_excel(args.color_coding_file_pca, sheet_name="Color_Coding", dtype=object)
        group_color_dict = dict(zip(group_color_df["Group"], group_color_df["Color"]))
        # Assign colors to genotypes based on their groups
        genotype_colors = {genotype: group_color_dict.get(genotype_group_dict.get(genotype), "black")
                       for genotype in PCA_data_frame_final["Genotype"].unique()}
        # create the PCA plot
        fig=px.scatter(PCA_data_frame_final, x='PC1', y='PC2', title='Principal Component Analysis', hover_data=["Genotype"], color="Genotype", color_discrete_map=genotype_colors, size_max=25)
        fig.update_traces(textposition='top center')
        fig.update_layout(
            xaxis=dict(title=f'PC1:{explained_variance[0]:.2%}', title_font=dict(size=24),tickfont=dict(size=20), range=x_range),
            yaxis=dict(title=f'PC2:{explained_variance[1]:.2%}', title_font=dict(size=24),tickfont=dict(size=20), range=y_range),
            title=dict(x=0.5, xanchor="center", yanchor="top", y=0.95, font=dict(size=24)),
            showlegend=False,
            plot_bgcolor="lightgray"
        )
        # save as interactive html
        out_pict3 = os.path.join(self.output_dir, "PCA.html").replace("\\", "/")
        fig.write_html(out_pict3)
    #
    def calculate_kinship_matrix(self):
        # retrieve matrix from object
        genotype_matrix = self.numeric_matrix
        # transpose the matrix
        genotype_matrix_transposed = genotype_matrix.T
        # filter for monomorphic marker
        filter_values = ["3","2","1","0"]
        for value in filter_values:
            genotype_matrix_transposed = genotype_matrix_transposed.loc[:,(genotype_matrix_transposed != value).any()]
        # mark values as integers
        genotype_matrix_transposed = genotype_matrix_transposed.astype(int)
        # centering by subtracting the row mean
        centered_genotype_matrix = genotype_matrix_transposed.sub(genotype_matrix_transposed.mean(axis=1), axis=0)
        # calculate kinship matrix
        kinship_matrix = np.dot(centered_genotype_matrix, centered_genotype_matrix.T)/genotype_matrix_transposed.shape[1]
        # convert to dataframe
        kinship_dataframe = pd.DataFrame(kinship_matrix, index=genotype_matrix_transposed.index,  columns=genotype_matrix_transposed.index)
        self.kinship_matrix = kinship_dataframe
    #
    def calculate_maf(self):
        numeric_matrix_transposed = self.numeric_matrix.T
        # remove monomorphic markers
        filter_values = [3,2,1,0]
        for value in filter_values:
            numeric_matrix_transposed = numeric_matrix_transposed.loc[:,(numeric_matrix_transposed != value).any()]
        # re-transpose the matrix
        numeric_matrix = numeric_matrix_transposed.T
        # count maf    
        maf_results = {}
        snp_list, maf_list = [], []
        for idx, row in numeric_matrix.iterrows():
            allele_counts = {0: 0, 1: 0, 2: 0, 3: 0}
            for value in row.values:
                allele_counts[value] = allele_counts.get(value, 0) + 1
            # calculate MAF
            maf = allele_counts[0] / sum(allele_counts.values())
            maf_results[idx] = maf
            # append to list
            snp_list.append(idx)
            maf_list.append(f"{maf:.4f}")
        # merge into data frame
        maf_data = pd.DataFrame({"SNP": snp_list, "MAF": maf_list})
        maf_data["MAF"] = maf_data["MAF"].astype(float)
        self.maf_data = maf_data
    #
    def make_hapmap(self):
        hapmap = self.raw_data
        # replace fails and ambigous codes
        characters_to_replace = {"R" : "", 
                                 "Y" : "", 
                                 "S" : "", 
                                 "W" : "", 
                                 "K" : "", 
                                 "M" : "", 
                                 "-" : ""}
        # list for the allele column that will be added
        allele_list = []
        failed_call = "/-"
        # iterate through rows
        for idx, row in hapmap.iterrows():
            # get unique characters in each row as list, skip first two columns
            unique_nucleotides = set(row.values)
            # join the unique characters into one string and replace fails and ambigous bases
            nucleotides = "".join(unique_nucleotides).translate(str.maketrans(characters_to_replace))
            # check if only one allele information is provided and add "/-" and add to the list
            if len(nucleotides) == 1:
                alleles = nucleotides + failed_call
            elif len(nucleotides) < 1:
                alleles = "-/-"
            # if both allele information are provided add it to the list
            else:
                alleles = ("/".join(nucleotides))
            allele_list.append(alleles)
        # add as column to dataframe at first position
        hapmap.insert(0, "Alleles", allele_list)
        # change allele coding
        alleles_to_replace = {"A" : "AA",
                              "G" : "GG",
                              "C" : "CC",
                              "T" : "TT",
                              "R" : "AG",
                              "Y" : "CT",
                              "S" : "GC",
                              "W" : "AT",
                              "K" : "GT",
                              "M" : "AC",
                              "-" : "NN"
                             }
        # replace the values in the selected matrix
        hapmap = hapmap.replace(alleles_to_replace)
        ## load physical data (will be added in the future)
        #marker_physical_positions = pd.read_excel("Faba_Positions.xlsx", sheet_name="Tabelle1")
        #marker_physical_positions.set_index("Name", inplace=True)
        #hapmap = pd.merge(hapmap, marker_physical_positions[["CHR", "POS"]], left_index=True, right_index=True, how="left")
        # add hapmap columns
        hapmap.insert(1, "strand", "+")
        hapmap.insert(2, "assembly", "NA")
        hapmap.insert(3, "center", "NA")
        hapmap.insert(4, "protLSID", "NA")
        hapmap.insert(5, "assayLSID", "NA")
        hapmap.insert(6, "panelLSID", "NA")
        hapmap.insert(7, "QCcode", "NA")
        #hapmap.drop(["CHR","POS"], axis=1, inplace=True)
        self.hapmap = hapmap
    #
    def save_results_to_excel(self):
        out_excel = os.path.join(self.output_dir, "Results.xlsx").replace("\\", "/")
        with pd.ExcelWriter(out_excel) as writer:
            self.marker_statistics.to_excel(writer, sheet_name="marker_statistics", index=False)
            self.marker_status.to_excel(writer, sheet_name="marker_status", index=False)
            self.maf_data.to_excel(writer, sheet_name="maf", index=False)
            self.numeric_matrix.to_excel(writer, sheet_name="numeric_matrix", index=True)
            self.kinship_matrix.to_excel(writer, sheet_name="kinship_matrix")
            self.hapmap.to_excel(writer, sheet_name="hapmap", index=True)
#
if __name__ == "__main__":
    chip_data = DataMatrix(input_file=args.input, sheet=0, output_dir=args.outputDir)
    print("\nLoading data ...")
    chip_data.load_data()
    chip_data.data_check()
    print("\nLoaded data successfully!")
    chip_data.create_output_directory()
    print("\nCalculating SNP and marker stats ...")
    chip_data.calc_snp_stats()
    chip_data.calc_marker_status()
    print("\nGenerating numeric matrix ...")
    chip_data.generate_numeric_matrix()
    print("\nPerforming PCA ...")
    chip_data.calculate_pca()
    print("\nCalculating kinship ...")
    chip_data.calculate_kinship_matrix()
    print("\nCalculating MAF ...")
    chip_data.calculate_maf()
    print("\nConverting to HAPMAP ...")
    chip_data.make_hapmap()
    print("\nSaving results ...")
    chip_data.save_results_to_excel()
