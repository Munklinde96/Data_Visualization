modifications = df["PTM"]
differentModifications = {}
for x in modifications:
    if pd.isnull(x):
        continue
    x = x.strip()
    mods = x.split(";")
    for mod in mods:
        if mod not in differentModifications:
            differentModifications[mod] = 1
        else:
            differentModifications[mod] +=1
for mod, count in differentModifications.items():
    print("mod is "+mod+" count is "+str(count))

#Column chart (a.k.a. Bar Chart) showing distribution of modifications per peptide
plt.figure(figsize=(16,8))
plt.bar(differentModifications.keys(), differentModifications.values())
plt.title("Distribution of modifications over all proteins")
plt.xlabel('Name of Modifications')
plt.ylabel('Total Count')
plt.xticks(rotation='vertical')

modificationsAndProteins = df[["PTM", "Protein ID"]]
differentModificationsByProtein = {}
for x, proteinId in modificationsAndProteins.itertuples(index=False):
    if pd.isnull(x):
        continue
    if proteinId not in differentModificationsByProtein:
        differentModificationsByProtein[proteinId] = {}
    mods = x.split(";")
    for mod in mods:
        if mod not in differentModificationsByProtein[proteinId]:
             differentModificationsByProtein[proteinId][mod] = 1
        else:
             differentModificationsByProtein[proteinId][mod] += 1
for protein, mods in differentModificationsByProtein.items():
    for mod, count in mods.items():
        print("protein is: " +str(protein)+ "mod is "+mod+" count is "+str(count))

plt.figure(figsize=(16,8))
plt.bar(differentModificationsByProtein[1].keys(), differentModificationsByProtein[1].values())
plt.title("Distribution of modifications over single protein")
plt.xlabel('Name of Modifications')
plt.ylabel('Total Count')
plt.xticks(rotation='vertical')