## preprocessing cBioportal data
## input: "data_clinical_patient.txt"
## output data is to be used by R
import numpy as np
import pandas as pd
import difflib
import os.path, sys, getopt
pd.set_option("display.max_column", 100)


def clean_patient_data(dataFile):
    # read cBioPortal data file
    patient_file = os.path.join(dataFile, "data_clinical_patient.txt")
    patient_df = pd.read_csv(patient_file, sep="\t", skiprows=4, index_col=0)
    patient_df = patient_df.loc[:, ["DFS_MONTHS", "DFS_STATUS", "OS_MONTHS", "OS_STATUS"]]
    patient_df.replace({'LIVING': 0, 'DECEASED': 1}, inplace=True)
    patient_df.replace({'DiseaseFree': 0, 'Recurred/Progressed': 1}, inplace=True)
    patient_df.replace('[Not Available]', np.nan, inplace=True)
    patient_df.dropna(axis=0, how='any', inplace=True)
    return patient_df


def clean_expression_data(dataFile, gene):
    expression_file = os.path.join(dataFile, "data_RNA_Seq_v2_expression_median.txt")
    exp_df = pd.read_csv(expression_file, sep="\t", header=0, index_col=0)
    if gene in exp_df.index:
        gene_exp = exp_df.loc[exp_df.index == gene, :]
        return gene_exp.T
    else:
        print "{} is not found.".format(gene)
        return None


def jointables(patient_df, expression_df):
    patient_df.index = patient_df.index.map(lambda x: difflib.get_close_matches(x, expression_df.index)[0])
    df = patient_df.join(expression_df, how="left")
    return df


def main():
    # Read file and arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:d:o:", ["gene=", "data=", "output="])
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"wge
        sys.exit(2)
    for o, a in opts:
        if o in ["-g", "--gene"]:
            gene = a.upper()
        elif o in ("-d", "--data"):
            dataFile = a
        elif o in ("-o", "--output"):
            outputFile = a
        else:
            assert False, "unhandled option"

    patient_df = clean_patient_data(dataFile)
    expression_df = clean_expression_data(dataFile, gene)
    total_df = jointables(patient_df, expression_df)
    total_df.sort_values(gene, ascending = False, inplace=True)
    total_df.to_csv(outputFile, sep="\t")
    pass



if __name__ == "__main__":
    main()


# python src/survival_mRNA_cBio.py -g LILRB1 -d /Users/yiwenbu/Desktop/cBioPortal/Leukemia/laml_tcga_pub -o ./temp/LILRB1_survival_leukemia.xls
# python src/survival_mRNA_cBio.py -g CMTM8 -d /Users/yiwenbu/Desktop/cBioPortal/Leukemia/laml_tcga_pub -o ./temp/CMTM8_survival_leukemia.xls