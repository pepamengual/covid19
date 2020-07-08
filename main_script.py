import numpy as np
from scripts.align_and_find_mutations import align_and_fix

def read_fasta_file(fasta_file_path):
    protein_data, protein_info = {}, {}
    with open(fasta_file_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                line = line[1:].split("|")
                protein_name, date, patient_id, country = line[0], line[2], line[3], line[-1]
                protein_data.setdefault(protein_name, {}).setdefault(patient_id, "")
                protein_info.setdefault(protein_name, {}).setdefault(patient_id, (date, country))
                continue
            sequence = line
            protein_data[protein_name][patient_id] += sequence
    return protein_data, protein_info

def sequence_checker(sequence):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    amino_acid_list = [i for i in amino_acids]
    sequence_list = [i for i in sequence]
    difference = set(sequence_list) - set(amino_acid_list)
    if len(difference) > 0:
        return False
    else:
        return True

def group_sequences(protein_data):
    sequences_all, sequences_filtered = {}, {}
    for protein_name, patient_id_dictionary in protein_data.items():
        for patient_id, sequence in patient_id_dictionary.items():
            sequences_all.setdefault(protein_name, set()).add(sequence)
            result = sequence_checker(sequence)
            if result == True:
                sequences_filtered.setdefault(protein_name, []).append(sequence)
    return sequences_all, sequences_filtered

def align_reference_query(sequences_all, sequences_filtered):
    for protein, sequence_set in sequences_filtered.items():
        reference = max(set(sequence_set), key = sequence_set.count)
        for query in sequences_all[protein]:
            mutations, aligned_reference, aligned_query, regenerated_sequence = align_and_fix(reference, query)
            if mutations:
                print(mutations)
                print(aligned_reference)
                print(aligned_query)
                print(regenerated_sequence)
                print("----")

def main():
    fasta_file_path = "../GISAID/allprot0707.fasta"
    protein_data, protein_info = read_fasta_file(fasta_file_path)
    sequences_all, sequences_filtered = group_sequences(protein_data)
    align_reference_query(sequences_all, sequences_filtered)
main()
