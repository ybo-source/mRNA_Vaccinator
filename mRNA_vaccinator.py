def translate(seq):
    genetic_code = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += genetic_code[codon]
    return protein
def codon_optimizer(seq):
    codon_optimized_genetic_code = {'_': 'TGA', 'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 'E': 'GAG', 'F': 'TTC',
                                    'G': 'GGC', 'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 'L': 'CTG', 'M': 'ATG',
                                    'N': 'AAC', 'P': 'CCC', 'Q': 'CAG', 'R': 'AGA', 'S': 'AGC', 'T': 'ACC',
                                    'V': 'GTG', 'W': 'TGG', 'Y': 'TAC'}
    protein_sequence = translate(seq)
    optimized_sequence = ""
    for i in protein_sequence:
        optimized_sequence += codon_optimized_genetic_code[i]
    return optimized_sequence
def gc_optimizer(sequence):
    gc_optimized_genetic_code = {'I': ['ATC'], 'M': ['ATG'], 'T': ['ACC', 'ACG'], 'N': ['AAC'], 'K': ['AAG'],
                              'S': ['AGC', 'TCC'], 'R': ['CGC', 'CGG'], 'L': ['CTC', 'CTG'], 'P': ['CCC', 'CCG'],
                              'H': ['CAC'], 'Q': ['CAG'], 'V': ['GTC', 'GTG'], 'A': ['GCC', 'GCG'], 'D': ['GAC'],
                              'E': ['GAG'], 'G': ['GGC', 'GGG'], 'F': ['TTC'], 'Y': ['TAC'], '_': ['TAG', 'TGA'],
                              'C': ['TGC'], 'W': ['TGG']}
    protein_sequence = translate(sequence)
    optimized_sequence = ""
    for i in protein_sequence:
        optimized_sequence += gc_optimized_genetic_code[i][0]
    return optimized_sequence
def gc_content_finder(sequence):
    sequence = sequence.upper()
    frequency = 0
    for i in sequence:
        if i == "G" or i == "C":
            frequency += 1
    return (frequency/len(sequence)) * 100
def len_sorter(dict):
    len_sorted = sorted(dict.items(),reverse=True,key = lambda x: len(x[1]))
    return len_sorted
def gc_sorter(dict):
    gc_sorted = sorted(dict.items(), reverse=True, key=lambda x: gc_content_finder(x[1]))
    return gc_sorted
def fatsa_reader(file):
    with open(file, "r") as new_file:
        sequence = ""
        for i in new_file:
            sequence += i.strip()
    return sequence
def mRNA_maker(five_utr,covid_sequence,three_utr,Boolean = False):
    if Boolean:
        covid_sequence = covid_sequence.replace("T","U")
    return five_utr + (covid_sequence) + three_utr
def start_stop_checker (seq):
    if not seq.startswith("ATG"):
        seq = "ATG" + seq
    if not seq.endswith("TGA") and not seq.endswith("TAA") and not seq.endswith("TAG"):
        seq = seq + "TGA"
    return seq

input_file = "/Users/yusufberkoruc/Desktop/PycharmProjects/Yusuf/mRNA Project/Main/cds.fna"
dict_1 = dict()
with open(input_file, "r") as new_file:
    for i in new_file:
        if i.startswith(">"):
            splitted = i.split(" ",maxsplit = 1)[1].split("[",maxsplit= 1)[0].strip()
            dict_1[splitted] = str()
        else:
            dict_1[splitted]+=str(i.strip())


dict_2 = {key:start_stop_checker(value) for (key,value) in dict_1.items()}
gc_optimized_dict_1 = {key:gc_optimizer(value) for (key,value) in dict_1.items()}
codon_optimized_dict_1 = {key:codon_optimizer(value) for (key,value) in dict_1.items()}


trial_sequence = fatsa_reader("/Users/yusufberkoruc/Desktop/PycharmProjects/Yusuf/mRNA Project/Trial_sequence")
luciferase_sequence = fatsa_reader("/Users/yusufberkoruc/Desktop/PycharmProjects/Yusuf/mRNA Project/Luciferase_sequence")
fiveutr = fatsa_reader("/Users/yusufberkoruc/Desktop/PycharmProjects/Yusuf/mRNA Project/Main/five_utr")
threeutr = fatsa_reader("/Users/yusufberkoruc/Desktop/PycharmProjects/Yusuf/mRNA Project/Main/three_utr")

mRNA_library = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in dict_1.items()}
mRNA_library_gc_optimized = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in gc_optimized_dict_1.items()}
mRNA_library_codon_optimized = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in codon_optimized_dict_1.items()}
mRNA_library_len_sorted = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in len_sorter(dict_1)}
mRNA_library_gc_sorted = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in gc_sorter(dict_1)}
mRNA_library_gc_optimized_len_sorted = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in len_sorter(gc_optimized_dict_1)}
mRNA_library_gc_optimized_gc_sorted = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in gc_sorter(gc_optimized_dict_1)}
mRNA_library_codon_optimized_gc_sorted = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in len_sorter(codon_optimized_dict_1)}
mRNA_library_codon_optimized_len_sorted = {key:mRNA_maker(fiveutr,start_stop_checker(value),threeutr) for (key,value) in gc_sorter(codon_optimized_dict_1)}


print(mRNA_library)
print(mRNA_library_gc_optimized)
print(mRNA_library_codon_optimized)
print(mRNA_library_len_sorted)
print(mRNA_library_gc_sorted)
print(mRNA_library_gc_optimized_len_sorted)
print(mRNA_library_gc_optimized_gc_sorted)
print(mRNA_library_codon_optimized_len_sorted)
print(mRNA_library_codon_optimized_gc_sorted)


"""for i,k in dict_2.items():
    if k.startswith("ATG"):
        print(i,"have start codon",k)
    else:
        print(i,"doesn't have the start codon",k)

print("\n")
for i,k in dict_2.items():
    if k.endswith("TGA") or k.endswith("TAA") or k.endswith("TAG"):
        print(i,"have stop codon",k[::-1])
    else:
        print(i,"doesn't have the stop codon",k[::-1])"""