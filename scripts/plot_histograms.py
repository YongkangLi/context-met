import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the datasets
data_true = pd.read_csv('example_output/true_number_of_mutations.txt', header=None, names=['Values'])
data_mh_hobolth = pd.read_csv('example_output/mh_hobolth_number_of_mutations.txt', header=None, names=['Values'])
data_mh_indpt = pd.read_csv('example_output/mh_indpt_number_of_mutations.txt', header=None, names=['Values'])
data_hobolth = pd.read_csv('example_output/uncorrected_number_of_mutations.txt', header=None, names=['Values'])

# Combine all datasets into one DataFrame with an additional column for the labels
data_true['Label'] = 'Exact Samples'
data_mh_hobolth['Label'] = 'MH by Hobolth Proposal'
data_mh_indpt['Label'] = 'MH by Independent Proposal'
data_hobolth['Label'] = 'Hobolth Sampling Algorithm'

# Concatenate all data into a single DataFrame
data_all = pd.concat([data_true, data_mh_indpt, data_mh_hobolth , data_hobolth])

# Plot the histograms
plt.figure(figsize=(8, 6))
bin_edges = np.arange(1, data_all['Values'].max() + 1, 1)  # Start bins from 1 explicitly
hist_plot = sns.histplot(data_all, x='Values', hue='Label', multiple='dodge', bins=bin_edges, shrink=0.8)

plt.xlabel('Number of Mutations', fontsize=20)
plt.ylabel('Frequency', fontsize=20)

# Set xlim to start from 1
# plt.xlim(left=1, right=data_all['Values'].max() - 2)
plt.xlim(1, 8)

# Set xticks to exclude 0 and start from 1
# plt.xticks(np.arange(1, data_all['Values'].max() + 1, 1), fontsize=18)  # Force x-ticks to start from 1
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Access the automatically created legend by Seaborn and remove the title
legend = hist_plot.get_legend()
legend.set_title(None)  # This removes the "Label" title
for text in legend.get_texts():
    text.set_fontsize(16)  # Increase the font size of legend labels

plt.tight_layout()  # Ensure the layout fits tightly
plt.savefig('example_output/number_of_mutations_histograms.pdf')
plt.show()
