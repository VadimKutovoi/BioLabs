codons = {"AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R"}
data = "aabaacaad"

def protein_translation(data):
    res = ""
    for i in range(0, len(data), 3):
        print(i)
        print(data[i:i+3])
        if data[i:i+3] in codons.keys():
            res += codons[data[i:i+3]]
    return res

print(protein_translation(data))