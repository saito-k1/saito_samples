"""
File: clean_data.py
History: 14-Sep-2022
This program cleans the scRNA-seq output files.
"""

import os
import argparse
import pandas as pd
import shutil


def main():
    """Business logic"""

    args = get_cli_args()
    output_dir = args.output_dir.casefold()
    remove_insig_snv(output_dir)


def get_cli_args():
    """
    Just get the command line options using argparse
    @return: Instance of argparse arguments
    """
    parser = argparse.ArgumentParser(description='Give the name of the output directory')

    parser.add_argument('-output_dir',
                        type=str, help='Name of output directory', required=False,
                        default='outputs')

    return parser.parse_args()


def remove_insig_snv(output_dir):
    """ This function removes SNV directories if the TSV does not have any significant
    genes.
    :param: output_dir: name of the output directory where SNV directories are
    :return: filtered_tsv_list
    """

    # since tsv files are in subdirectories, use walk
    tsv_list = []
    no_sig_genes = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith(".tsv"):
                tsv_list.append(os.path.join(root, file))
        # loop through tsv list and find snvs that do not have any significant p-values
        for snv in tsv_list:
            try:
                df = pd.read_table(f'{snv}', index_col=0)
                remove_cols = [col for col in df.columns if 'Unnamed' in col]
                df.drop(remove_cols, axis='columns', inplace=True)
                if all(df['padj'] <= 0.05):
                    continue
                snv_filtered = df[df['padj'] <= 0.05]
                # if none of the genes have a padj less than 0.05, then remove the folder
                if len(snv_filtered) == 0:
                    shutil.rmtree(root)
                    no_sig_genes.append(snv)
                    print(f"{snv} does not have any significant genes, removed directory")
                else:
                    # rewrite tsv with only genes that have a padj less than 0.05
                    snv_filtered.to_csv(f'{snv}', sep="\t")
                    print(f"{snv} has sig genes: rewrote tsv with only sig genes.")
            except OSError:
                continue

    return tsv_list, no_sig_genes


if __name__ == '__main__':
    main()
