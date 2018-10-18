def count_of_patterns(pattern, genome)
	res = int()
	for i in range(len(genome) - len(pattern)):
    if pattern == genome[i:i + len(pattern)]:
        res += 1


