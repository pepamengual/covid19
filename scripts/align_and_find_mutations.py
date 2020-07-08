from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

def align_sequences(reference, query, **kwargs):
    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)

    alignments = pairwise2.align.globalds(reference, query,
                                          matrix, gap_open, gap_extend,
                                          penalize_end_gaps=(False, False))

    best_alignment = alignments[0]
    aligned_reference, aligned_query, score, begin, end = best_alignment
    return aligned_reference, aligned_query

def identify_and_regenerate(aligned_reference, aligned_query):
    mutation_list = []
    regenerated_sequence = []
    for i, (aa_r, aa_q) in enumerate(zip(aligned_reference, aligned_query)):
        if aa_q == "X":
            aa_q = aa_r
        if aa_r == aa_q:
            regenerated_sequence.append(aa_r)
        else:
            if aa_r == "-": #insertion
                name = "_{}{}".format(i, aa_q)
                regenerated_sequence.append("-") #this is wrong, please check me
            elif aa_q == "-": #deletion
                name = "{}{}_".format(aa_r, i)
            else: #substitution
                name = "{}{}{}".format(aa_r, i, aa_q)
                regenerated_sequence.append(aa_q)
            mutation_list.append(name)
    mutation_names = ";".join(mutation_list)
    regenerated_sequence = "".join(regenerated_sequence)
    return mutation_names, regenerated_sequence


def align_and_fix(reference, query):
    aligned_reference, aligned_query = align_sequences(reference, query)
    mutation_names, regenerated_sequence = identify_and_regenerate(aligned_reference, aligned_query)
    return mutation_names, aligned_reference, aligned_query, regenerated_sequence
