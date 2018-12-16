# Task 1

def CountOfPatterns(pattern, genome):
    res = int()
    for i in range(len(genome) - len(pattern)):
        if pattern == genome[i:i + len(pattern)]:
            res += 1
    return res

def FrequentWords(genome, k):
    patterns = {}
    frequent_patterns = []
    for i in range(len(genome) - k):
        if genome[i: i + k] not in patterns:
            patterns.update({genome[i: i + k]: count_of_patterns(genome[i: i + k], genome)})
    maxcount = max(patterns.values())
    for key, value in patterns.items():
        if value == maxcount:
            frequent_patterns.append(key)
    return frequent_patterns

def ReverseComplement(genome):
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

# Task 2

def ProteinTranslation(data):
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

def PeptideEncoding(dna, peptide):
    rna = dna.replace("T", "U")

    for i in range(0, len(dna) - len(peptide)*3, 1):
        if protein_translation(rna[i:i + len(peptide)*3]) == peptide:
            print(dna[i:i + len(peptide)*3])
        elif protein_translation(reverse_complement(dna[i:i + len(peptide) * 3]).replace("T", "U")) == peptide:
            print(dna[i:i + len(peptide)*3])


def Subpepetides_count(data):
    _len = int(data)
    subpeptides = _len * (_len - 1)
    print(subpeptides)


def Spectrum(peptide):
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
            for k in tmp:
                loc_sum += mass.get(k)
            result.append(loc_sum)
    result.sort()
    return result


def Peptide_count(mass):
    mass = int(mass)
    acid_mass = [57, 71, 87, 97, 99, 101, 103, 113,
                 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    res = dict()
    res[0] = 1
    for i in range(57, mass + 1):
        res[i] = 0
        for j in acid_mass:
            if(i - j) in res:
                res[i] = res[i - j] + res[i]
    return res[mass]

# Task 3

def Score(peptide, input_spectrum):
    score_ = 0
    spectrum_ = Spectrum(peptide)
    for i in spectrum_:
        if i in input_spectrum:
            score_ += 1
    print(score_)

def LinearScore(peptide, input_spectrum):
    score = 0
    spectrum = LinearSpectrum(peptide)
    for i in input_spectrum:
        if i in spectrum:
            score += 1
            spectrum.remove(i)
    return score

def Consistent(peptide, spectrum):
    masses = {
        'A': 71, 'I': 113, 'N': 114, 'D': 115, 'C': 103, 'Q': 128, 'E': 129,
        'G': 57, 'H': 137, 'L': 113, 'K': 128, 'M': 131, 'F': 147, 'P': 97,
        'S': 87, 'T': 101, 'W': 186, 'Y': 163, 'V': 99, 'R': 156}
    linear_spec = LinearSpectrum(peptide)
    spectrum_ = list(spectrum)
    for i in linear_spec:
        if i in spectrum_:
            spectrum_.remove(i)
        else:
            return False
    return True


def Mass(peptide):
    mass = 0
    masses = {
        'A': 71, 'I': 113, 'N': 114, 'D': 115, 'C': 103, 'Q': 128, 'E': 129,
        'G': 57, 'H': 137, 'L': 113, 'K': 128, 'M': 131, 'F': 147, 'P': 97,
        'S': 87, 'T': 101, 'W': 186, 'Y': 163, 'V': 99, 'R': 156}
    for i in peptide:
        mass += masses[i]
    return mass

def ParentMass(spectrum):
    return spectrum[-1]

def Expand(peptides):
    masses = {
        'A': 71, 'I': 113, 'N': 114, 'D': 115, 'C': 103, 'Q': 128, 'E': 129,
        'G': 57, 'H': 137, 'L': 113, 'K': 128, 'M': 131, 'F': 147, 'P': 97,
        'S': 87, 'T': 101, 'W': 186, 'Y': 163, 'V': 99, 'R': 156}
    res = []
    for peptide in peptides:
        for el in masses:
            tmp = list(peptide)
            tmp.append(el)
            res.append(tmp)
    return res

def LinearSpectrum(peptide):
    spectrum = [0]
    for i, val in enumerate(peptide):
        for j in range(len(peptide) - i):
            spectrum.append(Mass(peptide[i:i + j + 1]))
    spectrum.sort()
    return spectrum

def Trim(leaderBoard, spectrum, n):
	boardLen = len(leaderBoard)
	if boardLen <= n:
		return leaderBoard
	scoreBoard = [0] * boardLen
	for index, peptide in enumerate(leaderBoard):
		scoreBoard[index] = LinearScore(peptide, spectrum)
	scoreBoard, leaderBoard = (list(t) for t in zip(*sorted(zip(scoreBoard, leaderBoard), reverse = True)))
	retBoard = leaderBoard[:n]; i = n
	while i < boardLen and scoreBoard[i] == scoreBoard[n-1]:
		retBoard.append(leaderBoard[i])
		i += 1
	return retBoard 	

def PeptideSequencing(spectrum):
    masses = {
        'A': 71, 'I': 113, 'N': 114, 'D': 115, 'C': 103, 'Q': 128, 'E': 129,
        'G': 57, 'H': 137, 'L': 113, 'K': 128, 'M': 131, 'F': 147, 'P': 97,
        'S': 87, 'T': 101, 'W': 186, 'Y': 163, 'V': 99, 'R': 156}
    peptides = ['']
    res_pept = []
    while len(peptides) != 0:
        peptides = Expand(peptides)
        tmp = 0
        for i in range(len(peptides)):
            i -= tmp
            if Mass(peptides[i]) == ParentMass(spectrum):
                if Spectrum(peptides[i]) == spectrum:
                    res_pept.append(peptides[i])
                    peptides.remove(peptides[i])
                    tmp += 1
            elif not Consistent(peptides[i], spectrum):
                peptides.remove(peptides[i])
                tmp += 1
    res = []
    for pept in res_pept:
        tmp = ''
        for acid in pept:
            tmp += str(Mass(acid)) + '-'
        tmp = tmp[:len(tmp) - 1]
        res.append(tmp)
    return set(res)


def LeaderboardSequencing(spectrum, n):
    leaderboard = ['']
    leaderpeptide = []
    while len(leaderboard) != 0:
        leaderboard = Expand(leaderboard)
        tmp = 0
        for i in range(len(leaderboard)):
            i -= tmp
            if Mass(leaderboard[i]) == ParentMass(spectrum):
                if LinearScore(leaderboard[i], spectrum) > LinearScore(leaderpeptide, spectrum):
                    leaderpeptide = list(leaderboard[i])
            elif Mass(leaderboard[i]) > ParentMass(spectrum):
                    leaderboard.remove(leaderboard[i]); tmp += 1
        leaderboard = Trim(leaderboard, spectrum, n)
    pmass = [0] * len(leaderpeptide) 
    res = ''
    for pept in leaderpeptide[::-1]:
        res += str(Mass(pept)) + '-'
    return res[:len(res) - 1]


def main():
    in_pept = 'NQEL'
    in_spec = '0 99 113 114 128 227 257 299 355 356 370 371 484'
    inp = '0 113 128 186 241 299 314 427'.split(" ")
    inp = list(map(int, inp))
    list_spec = in_spec.split(" ")
    int_list = [0]*len(list_spec)
    for i, val in enumerate(list_spec):
        int_list[i] = int(val)

    print(Peptide_sequencing(inp))
    
    
main()
