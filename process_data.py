import numpy as np

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

def identify_unique_sequences(protein_data):
    sequences = {}
    for protein_name, patient_id_dictionary in protein_data.items():
        for patient_id, sequence in patient_id_dictionary.items():
            result = sequence_checker(sequence)
            if result == True:
                sequences.setdefault(protein_name, []).append(sequence)
    for protein, sequence_set in sequences.items():
        mean = int(np.mean([len(i) for i in set(sequence_set)]))
        print("{} protein: {} raw sequences, {} unique sequences, {} average amino acids".format(protein, len(sequence_set), len(set(sequence_set)), mean))

def main():
    fasta_file_path = "../GISEAD/allprot0707.fasta"
    protein_data, protein_info = read_fasta_file(fasta_file_path)
    identify_unique_sequences(protein_data)
main()
