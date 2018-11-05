
def protein_translation(data):
    res = ""
    codons = {
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '', 'UAG': '', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UGU': 'C', 'UGC': 'C', 'UGA': '', 'UGG': 'W', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M'}
    for i in range(0, len(data), 3):
        if data[i:i+3] in codons.keys():
            res += codons[data[i:i+3]]
    return res

def reverse_complement(genome):
    res = ''
    for i in genome[::-1]:
        if i == 'A':
            res += 'T'
        elif i == 'T':
            res += 'A'
        elif i == 'G':
            res += 'C'
        elif i == 'C':
            res += 'G'
    return res


def peptide_encoding(dna, peptide):
    rna = dna.replace("T", "U")

    for i in range(0, len(dna) - len(peptide)*3, 1):
        if protein_translation(rna[i:i + len(peptide)*3]) == peptide:
            print(dna[i:i + len(peptide)*3])
        elif protein_translation(reverse_complement(dna[i:i + len(peptide) * 3]).replace("T", "U")) == peptide:
            print(dna[i:i + len(peptide)*3])


def subpepetides_count(data):
    _len = int(data)
    subpeptides = _len * (_len - 1)
    print(subpeptides)


def spectrum(peptide):
    result = [0]
    mass = {
        'A': 71, 'I': 113, 'N': 114, 'D': 115, 'C': 103, 'Q': 128, 'E': 129,
        'G': 57, 'H': 137, 'L': 113, 'K': 128, 'M': 131, 'F': 147, 'P': 97,
        'S': 87, 'T': 101, 'W': 186, 'Y': 163, 'V': 99, 'R': 156}

    total_sum = 0
    for i in peptide:
        total_sum += mass[i]
        result.append(mass[i])
    result.append(total_sum)

    for i in range(len(peptide)):
        for j in range(1, len(peptide) - 1):
            loc_sum = 0
            if i + j > (len(peptide) - 1):
                tmp = peptide[i:len(peptide)] + peptide[:j - (len(peptide) - i) + 1]
            else:
                tmp = peptide[i:i + j + 1]
            for k in range(len(tmp)):
                loc_sum += mass[tmp[k:k + 1]]
            result.append(loc_sum)

    result.sort()
    return result

def peptide_count(mass):
    mass = int(mass)
    acid_mass = [57, 71, 87, 97, 99, 101, 103, 113,
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    res = {}
    res[0] = 1
    for i in range(57, mass + 1):
        res[i] = 0
        for j in acid_mass:
            if(i - j) in res:
                res[i] = res[i - j] + res[i]
    return res[mass]

