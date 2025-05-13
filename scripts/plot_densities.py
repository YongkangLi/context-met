import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the datasets
data_true = pd.read_csv('example_output/true_first_jump_time.txt', header=None, names=['Values'])
data_mh_indpt = pd.read_csv('example_output/mh_indpt_first_jump_time.txt', header=None, names=['Values'])
data_mh_hobolth = pd.read_csv('example_output/mh_hobolth_first_jump_time.txt', header=None, names=['Values'])
data_hobolth = pd.read_csv('example_output/uncorrected_first_jump_time.txt', header=None, names=['Values'])

# Plot the density plots
plt.figure(figsize=(8, 6))
sns.kdeplot(data_true['Values'], color='blue', label='Exact Samples')
sns.kdeplot(data_mh_indpt['Values'], color='orange', label='MH by Independent Proposal')
sns.kdeplot(data_mh_hobolth['Values'], color='green', label='MH by Hobolth Proposal')
sns.kdeplot(data_hobolth['Values'], color='red', label='Hobolth Sampling Algorithm')
# sns.kdeplot(data_indpt['Values'], color='purple', label='JC Model')

plt.xlabel('First Jump Time', fontsize=20)
plt.ylabel('Density', fontsize=20)
plt.legend(fontsize=14)
plt.xlim(left=0.0, right=0.06)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.savefig('example_output/first_jump_time_densities.pdf')
plt.show()
