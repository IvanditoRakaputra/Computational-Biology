from collections import Counter
from itertools import product

def complement_dna(dna_seq):
    seqm = dna_seq.upper()
    conversion = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = ''.join([conversion[base] for base in seqm])
    return complement

def transcript_dna(dna_seq):
    seqm = dna_seq.upper()
    conversion = {'C': 'G', 'G': 'C', 'A': 'U', 'T': 'A'}
    mrna = ''.join([conversion[base] for base in seqm])
    return mrna

def translate_mrna(mrna_seq):
    amino_acids = []
    codon_to_amino_acid = {
        'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'STOP', 'UAG': 'STOP',
        'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'STOP', 'UGG': 'Trp',
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    for i in range(0, len(mrna_seq), 3):
        codon = mrna_seq[i:i+3]
        if len(codon) == 3:
            amino_acid = codon_to_amino_acid.get(codon, '')
            if amino_acid == 'STOP':
                break
            if amino_acid:
                amino_acids.append(amino_acid)
    
    return '-'.join(amino_acids)

# Function to find all possible codon combinations for a given amino acid sequence
def codon_combinations_for_amino_acids(amino_acids):
    amino_acid_to_codons = {
        'Tyr': ['UAU', 'UAC'],  
        'Trp': ['UGG'],         
    }
    codon_combinations = list(product(*[amino_acid_to_codons[aa] for aa in amino_acids]))
    codon_freq_dict = {}
    for codon_comb in codon_combinations:
        mrna_seq = ''.join(codon_comb)
        codons_list = [mrna_seq[i:i+3] for i in range(0, len(mrna_seq), 3)]
        codon_freq = Counter(codons_list)
        codon_freq_dict[mrna_seq] = dict(codon_freq)
    
    return codon_freq_dict

# PART 1: DNA -> mRNA -> Amino Acids (and complement)
print("=== Part 1: DNA to mRNA and Amino Acid Sequence ===")
dna_sequence = "TTACGA"
complement = complement_dna(dna_sequence)
mrna = transcript_dna(dna_sequence)
amino_acid_sequence = translate_mrna(mrna)

print(f"Input DNA: {dna_sequence}")
print(f"Complement: {complement}")
print(f"mRNA: {mrna}")
print(f"Amino Acids: {amino_acid_sequence}")

# Separator 
print("\n" + "="*50 + "\n")

# PART 2: Amino Acid Sequence -> mRNA -> Codon Frequency
print("=== Part 2: Amino Acid Sequence and Codon Frequencies ===")
amino_acid_input = ['Trp', 'Tyr', 'Trp']  # Corresponds to 'WYW'
codon_freqs_for_combinations = codon_combinations_for_amino_acids(amino_acid_input)

for mrna_seq, codon_freq in codon_freqs_for_combinations.items():
    print(f"mRNA Sequence: {mrna_seq}")
    print(f"Codon Frequencies: {codon_freq}")
    print("-" * 30)
