import pandas as pd

# load file in as a tsv and convert to a csv 
df = pd.read_csv('jacques_wines.tsv', sep='\t')
df.to_csv('jacques_wines.csv', index=False)

print("First 5 rows of the dataset: \n", df.head(), "\n\n")

# print out the number of lines in the files (not including headers)
print("2) Number of lines in the dataset: ", len(df), "\n")

# print out the number of unique varieties in the dataset
print("3) Number of unique varieties in the dataset: ", df['variety'].nunique(), "\n")

# print the 12th most frequency variety in the dataset
print("4) 12th most frequency variety in the dataset: ", df['variety'].value_counts().index[11], "\n")

# print vineyard that procudes the greatest number of wines of different varieties in the dataset
print("5) Vineyard that procudes the greatest number of wines of different varieties in the dataset: ", df['winery'].value_counts().index[0], "\n")

variety_value_counts = df['variety'].value_counts()
winery_value_counts = df['winery'].value_counts()

print("\nvariety_value_counts: \n", variety_value_counts, "\n")
print("\nwinery_value_counts: \n", winery_value_counts, "\n")