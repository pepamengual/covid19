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

def identify_mutations(aligned_reference, aligned_query):
    mutation_list = []
    for i, (aa_r, aa_q) in enumerate(zip(aligned_reference, aligned_query)):
        if aa_r != aa_q: #different residue
            if aa_r != "X" and aa_q != "X": #different than unknown residue
                mutation = ""
                if aa_r == "-":
                    mutation = "ins"
                    name = "_{}{}".format(i, aa_q)
                elif aa_q == "-":
                    mutation = "del"
                    name = "{}{}_".format(aa_r, i)
                else:
                    name = "{}{}{}".format(aa_r, i, aa_q)
                mutation_list.append(name)
    mutation_names = ";".join(mutation_list)
    print(mutation_names)
    return mutation_names

def main():
    reference = "AFGIHKLKLLQARSX"
    query = "AFGPIHKKLLQARA"
    aligned_reference, aligned_query = align_sequences(reference, query)
    mutation_names = identify_mutations(aligned_reference, aligned_query)

main()
