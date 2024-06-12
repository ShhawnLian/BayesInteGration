import pandas as pd
import os

csv_path = 'BayesInteGration-main/Data/csvFile/'
save_path = 'BayesInteGration-main/Data/'

file_list = [file for file in os.listdir(csv_path) if file.endswith('.csv')]

dfs = []

# Iterate over each CSV file and append its data to the merged data DataFrame
for file in file_list:
    file_path = os.path.join(csv_path, file)
    data = pd.read_csv(file_path)
    dfs.append(data)

merged_data = pd.concat(dfs, ignore_index=True)
merged_data.to_csv(save_path + 'merged.csv', index=False)