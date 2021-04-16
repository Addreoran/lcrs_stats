# 1. Jakie typy LCRów istnieją.
# 2. Ile tego jest i jak się rozkłada filogenetycznie?
# 3. Jaki stanowią procent wszystkich białek zarówno pod względem sumy
# długości jak i zajęcia procentowego długości białek?
# 4. Czy są aminokwasy, które się uwielbiają i te, które nie lubią ze sobą być? - brak
# 5. Policzyć entropię shannona na całych sekwencjach i
# sprawdzić gdzie są LCRy na wykresie - jeszcze brak
# 6. procentowość aminokwasów w LCRach vs cała baza
# 7. w ilu % białek są LCRy
# 8. które narzędzia są dobre - porównać ze sobą metody do LCRów - brak
import matplotlib.pyplot as plt
import numpy as np
from ete3 import NCBITaxa


def find_all_protein(file):
    proteins = []
    with open(file) as f:
        for line in f.readlines():
            if line.startswith(">"):
                proteins.append(line.split("|")[1])
    return list(set(proteins))


def len_by_protein(file_in):
    hist_len = {}
    lcrs=0
    with open(file_in) as f:
        for line in f:
            if line.startswith(">"):
                last = line.split("|")[1]
            else:
                if last in hist_len:
                    hist_len[last] += len(line.strip())
                else:
                    hist_len[last] = len(line.strip())
                lcrs+=1
    print(lcrs)
    return hist_len


def find_repeats(repeats_in):
    proteins = {}
    with open(repeats_in) as f:
        for i in f.readlines():
            print(i)

            if i.split("|")[1] not in proteins:
                proteins[i.split("|")[1]] = []
            proteins[i.split("|")[1]].append(i.rsplit(";", 2)[-2])
    return {i: set(j) for i, j in proteins.items()}


def get_tax(lcrs, repeats_in):
    tax = {}
    repeats = find_repeats(repeats_in)
    with open(lcrs) as f:
        for line in f.readlines():
            if line.startswith(">"):
                protein = line.split("|")[1]
                taxon = line.split("OX=")[1].split()[0]
                if taxon not in tax.keys():
                    tax[taxon] = []
            else:
                if protein in repeats:
                    tax[taxon].append(repeats[protein])
    return tax


def research_2(file_lcrs, repeats_in):
    ncbi = NCBITaxa()
    # ncbi.update_taxonomy_database()
    tax = get_tax(file_lcrs, repeats_in)
    tree = ncbi.get_topology(tax.keys())
    names = ncbi.get_taxid_translator(tax.keys())
    # lieage = {t: ncbi.get_rank(ncbi.get_lineage(t)) for t in tax.keys()}
    tax_tmp = {i: {";".join(k): j.count(k) for k in j if len(k) == 2} for i, j in tax.items()}
    for leaf in tree.get_leaves():
        try:
            leaf.name = names[int(leaf.name)] + " " + str(tax_tmp[leaf.name])
        except:
            print(type(leaf.name), leaf.name, leaf.name in names)
            # print(names.keys())
    # print(tree)
    tree.show()
    pass


def research_5(file_lcrs, file_database):
    pass


def research_3(file_lcrs, file_database):
    lcrs_protein_lenths = len_by_protein(file_lcrs)
    proteins_lenths = len_by_protein(file_database)
    print("W swissprot jest ", len(proteins_lenths), " białek o długości ", sum(proteins_lenths.values()))
    print("W LCRach jest ", len(lcrs_protein_lenths), " białek o długości ", sum(lcrs_protein_lenths.values()))
    print("LCRy stanowią ", 100 * sum(lcrs_protein_lenths.values()) / sum(proteins_lenths.values()), " długości białek")
    proteins_with_lcrs = [proteins_lenths[i] for i in lcrs_protein_lenths.keys()]
    print("LCRy stanowią ", 100 * sum(lcrs_protein_lenths.values()) / sum(proteins_with_lcrs),
          " białek w których występują")
    # print(lcrs_protein_lenths)


def research_7(file_lcrs, file_database):
    lcrs = find_all_protein(file_lcrs)
    database = find_all_protein(file_database)
    print("Liczba białek z LCRami: ", len(lcrs))
    print("Liczba białek: ", len(database))
    print("Procent białek z LCRami: ", 100 * len(lcrs) / len(database))


def calculate_hist(file):
    hist = {}
    with open(file) as f:
        for line in f:
            if not line.startswith(">"):
                for aa in set(line.strip()):
                    if aa in hist.keys():
                        hist[aa] += line.count(aa)
                    else:
                        hist[aa] = line.count(aa)
    return hist


def tidy_aa(labels, aa_dict):
    return [100 * aa_dict.get(i, 0) / sum(aa_dict.values()) for i in labels]


def research_2_1(file_lcrs):
    repeats = find_repeats2(file_lcrs)
    tmp = {}
    for i, j in repeats.items():
        if len(i) in tmp.keys():
            tmp[len(i)] += j
        else:
            tmp[len(i)] = j
    labels = sorted(list(tmp.keys()))
    x = np.arange(len(labels))
    width = 0.2
    # for i in labels:
    #     print(i, tmp[i])
    plt.bar(x, [tmp[i] for i in labels], width)
    # plt.legend(['LCRs', 'Swissprot'])
    plt.xticks(x, labels)
    # plt.title("Percentage of all types aminoacids in LCRs and swissprot proteins.")
    # plt.xlabel("Aminoacids")
    # plt.ylabel("Percentage of aminoacids")
    plt.show()


def research_6(file_lcrs, file_database):
    lcrs_aa = calculate_hist(file_lcrs)
    database_aa = calculate_hist(file_database)
    # print(lcrs_aa)
    # print(database_aa)
    labels = list(set(list(lcrs_aa.keys()) + list(database_aa.keys())))
    x = np.arange(len(labels))
    width = 0.2
    plt.bar(x - 0.1, tidy_aa(labels, lcrs_aa), width)
    plt.bar(x + 0.1, tidy_aa(labels, database_aa), width)
    plt.legend(['LCRs', 'Swissprot'])
    plt.xticks(x, labels)
    plt.title("Percentage of all types aminoacids in LCRs and swissprot proteins.")
    plt.xlabel("Aminoacids")
    plt.ylabel("Percentage of aminoacids")
    plt.show()


def find_repeats2(file_in):
    proteins = {}
    with open(file_in) as f:
        for i in f.readlines():
            if i.split("|")[1] not in proteins:
                proteins[i.split("|")[1]] = []
            # print(i)
            proteins[i.split("|")[1]].append(i.rsplit(";", 2)[-2])
    tmp = []
    cos = {}
    for i, j in proteins.items():
        tmp += j
    for i in set(tmp):
        cos[i] = tmp.count(i)
    return cos


def research_1(file_in, aa=1, aa_min=None, aa_max=None):
    repeats_tmp = find_repeats2(file_in)

    if aa >= 7:
        labels = [i for i in repeats_tmp.keys() if len(i) >= aa]
    else:
        labels = [i for i in repeats_tmp.keys() if len(i) == aa]
    if aa_min and aa_max:

        if aa >= 7:
            labels = [i for i, j in repeats_tmp.items() if len(i) >= aa and j >= aa_min and j < aa_max]
        else:
            labels = [i for i, j in repeats_tmp.items() if len(i) == aa and j >= aa_min and j < aa_max]
    x = np.arange(len(labels))
    width = 0.2
    plt.bar(x, [repeats_tmp[i] for i in labels], width)
    plt.xticks(x, labels)
    plt.xticks(rotation=90)
    if aa > 4 and not aa_max:
        plt.tick_params(axis='x', labelsize=3)
    if aa > 7:
        plt.title(f"Plot with repeats at least lenth {aa}")
    else:
        plt.title(f"Plot with repeats of lenth {aa}")
    plt.xlabel("Type of repeats")
    plt.ylabel("Number of occurence")
    plt.show()
    print(repeats_tmp)


if __name__ == "__main__":
    # research_1("./data/repeats.csv", aa=1)
    # research_1("./data/repeats.csv", aa=2)
    # research_1("./data/repeats.csv", aa=2, aa_min=10, aa_max=500)
    # research_1("./data/repeats.csv", aa=2, aa_min=1, aa_max=10)
    # research_1("./data/repeats.csv", aa=3)
    # research_1("./data/repeats.csv", aa=4)
    # research_1("./data/repeats.csv", aa=5)
    # research_1("./data/repeats.csv", aa=6)
    # research_1("./data/repeats.csv", aa=7)
    # research_2("./data/seg_tmp.fasta", "./data/repeats.csv")
    # research_2_1("./data/repeats.csv")
    research_3("./data/seg_tmp.fasta", "./data/swiss.fasta")
    # research_7("./data/seg_tmp.fasta", "./data/swiss.fasta")
    # research_6("./data/seg_tmp.fasta", "./data/swiss.fasta")
