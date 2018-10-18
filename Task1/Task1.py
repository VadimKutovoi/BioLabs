def count_of_patterns(pattern, genome):
    res = int()
    for i in range(len(genome) - len(pattern)):
        if pattern == genome[i:i + len(pattern)]:
            res += 1
    return res

def frequent_words(genome, k):
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

def reverse_complement(genome):
    res = ''
    for i in range(len(genome)-1, -1, -1):
        if genome[i] == 'A':
            res += 'T'
        elif genome[i] == 'T':
            res += 'A'
        elif genome[i] == 'G':
            res += 'C'
        elif genome[i] == 'C':
            res += 'G'
    return res


print(count_of_patterns('ATAT','GATATATGCATATACTT'))
print(frequent_words('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4))
print(reverse_complement('AAAACCCGGT'))
