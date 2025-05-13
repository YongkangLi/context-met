import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Define column names since the header is removed
col_names = ['lag', 'autocorrelation']

# Load data from CSV files without headers, assigning column names
indpt_first_jump_site = pd.read_csv(
    'example_output/mh_indpt_auto_correlation.txt',
    header=None,
    names=col_names
)
hobolth_first_jump_site = pd.read_csv(
    'example_output/mh_hobolth_auto_correlation.txt',
    header=None,
    names=col_names
)

max_lag = 128

lag1_m1 = indpt_first_jump_site['lag'][:max_lag]
AC1_m1 = indpt_first_jump_site['autocorrelation'][:max_lag]

lag1_m2 = hobolth_first_jump_site['lag'][:max_lag]
AC1_m2 = hobolth_first_jump_site['autocorrelation'][:max_lag]

# Create a single figure
plt.figure(figsize=(8, 6))

# Plot for First Jump Time (manual dashed lines with markers using plt.plot)
plt.plot(lag1_m1, AC1_m1, 'o-', color='orange', label='Autocorrelation for first jump time (Independent proposal)')
plt.plot(lag1_m2, AC1_m2, 'o-', color='green', label='Autocorrelation for first jump time (Hobolth proposal)')

# Set labels, grid, and legend
plt.xlabel('Lag (Sweeps)', fontsize=20)
plt.ylabel('Autocorrelation', fontsize=20)
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

plt.xlim(0, 32)
plt.ylim(0, 0.002)

plt.xticks(fontsize=18)
plt.yticks(fontsize=16)

# Create custom legend with both line styles (solid and dashed)
legend_elements = [
    Line2D([0], [0], marker='o', color='orange', linestyle='-', label='Autocorrelation for first jump time (Independent proposal)'),
    Line2D([0], [0], marker='o', color='green', linestyle='-', label='Autocorrelation for first jump time (Hobolth proposal)'),
]

plt.legend(handles=legend_elements, loc='best', fontsize=12)

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig('example_output/autocorrelation_plot.pdf')
plt.show()