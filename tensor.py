import pandas as pd
import numpy as np

df = pd.read_csv("cfs.csv") # reads the file
#print(df.shape)
#print(df.head())

taxa = set()

# iterates through the 4 columns and adds them into taxa set
for col in ["t1", "t2", "t3", "t4"]:
    for name in df[col]:
        taxa.add(name)

taxa = sorted(taxa, key=lambda x: int(x[1:]))
n = len(taxa)
q = len(df)
#print(taxa) # Current taxas present in the csv
#print("n =", len(taxa))

# Creating dictionary and assigning each taxa its index
taxon_to_i = {}
for i, taxon in enumerate(taxa):
    taxon_to_i[taxon] = i
#print(taxon_to_i)

T = np.zeros((n, n, q), dtype=float) #rows, columns, quartets
#print(T.shape)

#row = df.iloc[0]   # first quartet row
#k = 0              # tensor layer 0
# a = row["t1"]
# b = row["t2"]
# c = row["t3"]
# d = row["t4"]
#
# print(a, b, c, d)
#
# ia = taxon_to_i[a]
# ib = taxon_to_i[b]
# ic = taxon_to_i[c]
# id = taxon_to_i[d]
# print(ia, ib, ic, id)
#
# cf12_34 = row["CF12_34"]
# cf13_24 = row["CF13_24"]
# cf14_23 = row["CF14_23"]
# print(cf12_34, cf13_24, cf14_23, cf14_23)
#
# T[ia,ib,k] += cf12_34
# T[ib,ia,k] += cf12_34
# T[ic,id,k] += cf12_34
# T[ic,id,k] += cf12_34
#
# T[ia,ic,k] += cf13_24
# T[ic,ia,k] += cf13_24
# T[ib,id,k] += cf13_24
# T[id,ib,k] += cf13_24
#
# T[ia,id,k] += cf14_23
# T[id,ia,k] += cf14_23
# T[ib,ic,k] += cf14_23
# T[ic,ib,k] += cf14_23
#
# print("Quartet:", a, b, c, d)
# print(T[:, :, k])


for k,row in df.iterrows():
    a = row["t1"]
    b = row["t2"]
    c = row["t3"]
    d = row["t4"]

    #print(a,b,c,d)
    ia = taxon_to_i[a]
    ib = taxon_to_i[b]
    ic = taxon_to_i[c]
    id = taxon_to_i[d]
    #print(ia, ib, ic, id)

    cf12_34 = row["CF12_34"]
    cf13_24 = row["CF13_24"]
    cf14_23 = row["CF14_23"]
    #print(cf12_34)

    T[ia,ib,k] += cf12_34
    T[ib,ia,k] += cf12_34
    T[ic,id,k] += cf12_34
    T[ic,id,k] += cf12_34

    T[ia,ic,k] += cf13_24
    T[ic,ia,k] += cf13_24
    T[ib,id,k] += cf13_24
    T[id,ib,k] += cf13_24

    T[ia,id,k] += cf14_23
    T[id,ia,k] += cf14_23
    T[ib,ic,k] += cf14_23
    T[ic,ib,k] += cf14_23

    print("Quartet:", a, b, c, d)
    print(T[:, :, k])