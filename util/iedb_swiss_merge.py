#!/usr/bin/env python
#
# Script which reads in an IEDB epitope dataset and cross references it
# with the Swiss-Prot database to create a training set.
import os, sys, argparse
import pandas as pd
import numpy as np
import subprocess
from itertools import groupby


# https://www.biostars.org/p/710/
def __fasta_iter(fasta_name):
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


def __create_swiss_prot_dataframe(output_path, swiss_db_path, protein_id_file):
    protein_fasta = output_path + "/" + "proteins.fasta"
    cmd = [
        "blastdbcmd", "-db", swiss_db_path, "-entry_batch", protein_id_file,
        "-out", protein_fasta
    ]
    print " ".join(cmd)
    with open(os.devnull, "w+") as devnull:
        subprocess.call(cmd, stdout=devnull, stderr=subprocess.STDOUT)
    data = []
    for (header, seq) in __fasta_iter(protein_fasta):
        header_fields = header.split('|')
        assert (header_fields[0] == 'gi')
        gi = header_fields[1]
        data.append((gi, seq))
    return pd.DataFrame(data, columns=['protein_id', 'protein_sequence'])


# Converts the merged IEDB/Swiss-Prot table into training data [(cleaved_sample, True), (uncleaved_sample, False), ...]
def __create_training_set(output_dir, iedb_swiss_merged, generated_sample_len):
    progress = ['/', '-', '\\', '-']
    data = []
    samples_seen = 0
    for (protein_id,
         row_ids) in iedb_swiss_merged.groupby("protein_id").groups.items():
        epitope_rows = iedb_swiss_merged.loc[row_ids]
        # Grab the protein sequence from the first row
        protein_sequence = epitope_rows['protein_sequence'].iloc[0]
        # Sort by the C-terminus ('end')
        sorted_epitopes = epitope_rows.sort_values(by="end")
        protein_sequence_len = len(protein_sequence)
        current_loc = 0
        for (i, epitope_sequence, start, end) in \
                sorted_epitopes[["epitope_sequence", "start", "end"]].itertuples():
            epitope_start = int(start) - 1
            epitope_end = int(end) - 1

            samples_seen += 1
            sys.stderr.write("\b%s" % (progress[samples_seen % 4]))

            epitope_sequence_len = len(epitope_sequence)
            epitope_midpoint = epitope_end - epitope_sequence_len / 2
            sample_start_end = lambda pos: (pos - (generated_sample_len / 2) + 1, pos + (generated_sample_len / 2) + 1)
            cleaved_sample_pos = sample_start_end(epitope_end)
            uncleaved_sample_pos = sample_start_end(epitope_midpoint)

            # Check if our samples are off the protein sequence or overlap (TODO: handle these)
            if (uncleaved_sample_pos[0] < current_loc or
                    cleaved_sample_pos[1] > protein_sequence_len):
                continue

            current_loc = cleaved_sample_pos[1]

            # Double check that the start and end positions are correct
            assert protein_sequence[epitope_start:epitope_end + 1] == epitope_sequence,\
                     "Epitope failed to align to protein"

            fetch_seq = lambda pos: protein_sequence[pos[0]:pos[1]]
            data.append((fetch_seq(cleaved_sample_pos), 1))
            data.append((fetch_seq(uncleaved_sample_pos), 0))

    cleavage_data = pd.DataFrame(data)
    cleavage_data.columns = ['sequence', 'is_cleaved']
    cleavage_data.to_csv(output_dir + "/" + "training_set.csv")


def create_iedb_swiss_dataset(epitope_csv, swiss_db_path, output_dir,
                              max_epitope_len, min_epitope_len,
                              generated_sample_len):
    # Load the CSV from IEDB (skipping the first line, [2:])
    epi = pd.DataFrame.from_csv(epitope_csv)[2:]
    # Rename the columns to remove whitespace
    epi.columns = [
        'type', 'epitope_sequence', 'start', 'end', 'chebi', 'syn', 'protein',
        'protein_id', 'organism', 'oid', 'comments'
    ]
    # Remove GI entries that start with prefix "SRC"
    epi = epi[epi.protein_id.str.startswith("SRC") == False]
    # Remove entries with '+' notation (note: looking into this, e.g. "PLNISLGDVVLY + DEAM(N3)")
    epi = epi[epi.epitope_sequence.str.find('+') == -1]
    # Remove the "GI:" prefix from the GIs provided by IEDB
    epi["protein_id"] = epi.protein_id.str.replace("GI:", "")
    # Drop any epitopes that are not the desired length
    iedb_epitopes = epi[(epi.epitope_sequence.str.len() >= min_epitope_len) \
                         & (epi.epitope_sequence.str.len() <= max_epitope_len)]\
             .loc[:, ["epitope_sequence", "protein_id", "start", "end"]]
    # Create a file with a list of unique protein (antigen) IDs for use with BLAST
    antigen_ids = iedb_epitopes["protein_id"].unique()
    protein_id_file = output_dir + "/" + "protein_ids.txt"
    np.savetxt(protein_id_file, antigen_ids, fmt="%s")
    num_iedb_epitopes = iedb_epitopes.shape[0]
    num_iedb_proteins = antigen_ids.shape[0]
    print "There are %d IEDB MHC-1 epitopes from %d unique proteins (antigens)." % (
        num_iedb_epitopes, num_iedb_proteins)
    swiss_prot_db = __create_swiss_prot_dataframe(output_dir, swiss_db_path,
                                                  protein_id_file)
    print "Found %d IEDB proteins in Swiss-Prot sequence database" % (
        swiss_prot_db.shape[0])
    # Merge the IEDB epitops with the SWISS protein data and drop any NaN values (proteins not in SWISSprot)
    iedb_swiss_merged = pd.merge(
        iedb_epitopes, swiss_prot_db, on="protein_id", how="left").dropna()
    # Since we're going to sort on position, convert to numeric columns
    iedb_swiss_merged["start"] = pd.to_numeric(iedb_swiss_merged[
        "start"]).astype(int)
    iedb_swiss_merged["end"] = pd.to_numeric(iedb_swiss_merged["end"]).astype(
        int)
    # Save the raw merged data for reference
    iedb_swiss_merged.to_csv(output_dir + "/" + "merged_iedb_swiss.csv")
    __create_training_set(output_dir, iedb_swiss_merged, generated_sample_len)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Read IEDB epitope, cross with Swiss-Prot, and output dataset'
    )
    parser.add_argument(
        '-e', '--epitope_csv', metavar="IEDB_CSV", required=True)
    parser.add_argument(
        '-o', '--output_dir', metavar='OUTPUT_DIR', required=True)
    parser.add_argument('-s', '--swissdb', metavar="PATH", required=True)
    parser.add_argument(
        '-a', '--max_epitope_len', metavar="LEN", default=12, type=int)
    parser.add_argument(
        '-i', '--min_epitope_len', metavar="LEN", default=5, type=int)
    parser.add_argument(
        '-g', '--generated_sample_len', metavar="LEN", default=28, type=int)

    args = parser.parse_args()
    create_iedb_swiss_dataset(args.epitope_csv, args.swissdb, args.output_dir,
                              args.max_epitope_len, args.min_epitope_len,
                              args.generated_sample_len)
