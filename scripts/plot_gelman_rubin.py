import pandas as pd
import matplotlib.pyplot as plt

# Define column names since the header is removed
col_names = ['step', 'GR']

# Load data from CSV files without headers, assigning column names
indpt_first_jump_time = pd.read_csv(
    'example_output/mh_indpt_gelman_rubin.txt', # Ensure this path is correct
    header=None,
    names=col_names
)
hobolth_first_jump_time = pd.read_csv(
    'example_output/mh_hobolth_gelman_rubin.txt', # Ensure this path is correct
    header=None,
    names=col_names
)

# --- For debugging: Print the first few rows to check data loading ---
# print("Independent Proposal Data:")
# print(indpt_first_jump_time.head())
# print("\nHobolth Proposal Data:")
# print(hobolth_first_jump_time.head())
# print("\n")
# --------------------------------------------------------------------

step_m1 = indpt_first_jump_time['step']
GR_m1 = indpt_first_jump_time['GR'] # This column will contain floats, NaN is handled

step_m2 = hobolth_first_jump_time['step']
GR_m2 = hobolth_first_jump_time['GR'] # This column will contain floats, NaN is handled

# Create a figure
plt.figure(figsize=(8, 6))

# Plot Gelman Rubin Shrink Factor
plt.plot(step_m1, GR_m1, color='orange', linestyle='--', label='GR for first jump time (Independent proposal)')
plt.plot(step_m2, GR_m2, color='green', linestyle='--', label='GR for first jump time (Hobolth proposal)')

# Set labels, limits, and legend
plt.xlabel('Step', fontsize=20)
plt.ylabel('Gelman Rubin Shrink Factor', fontsize=20)
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

# --- Important: Check these limits based on your data ---
plt.xlim(3, 64)
# The Y-limit below might be causing issues if your GR values are outside this range.
# For example, your sample data has a GR value of ~190 at step 2.
plt.ylim(1.00, 4.00)
# If lines are not appearing, try commenting out plt.ylim() to see auto-scaled limits
# or adjust the limits to encompass your actual data range.
# For example: plt.ylim(0.9, max(GR_m1.max(), GR_m2.max()) * 1.1) if you want to see all points
# Or, if you only care about values near 1:
# indpt_first_jump_time = indpt_first_jump_time[indpt_first_jump_time['GR'] <= 8]
# hobolth_first_jump_time = hobolth_first_jump_time[hobolth_first_jump_time['GR'] <= 8]
# Re-extract step_m1, GR_m1 etc. if you filter the data.
#----------------------------------------------------------

plt.legend(fontsize=14)

plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# Adjust layout and save the figure
plt.tight_layout()
# Make sure the directory 'experiments/gelman_rubin/' exists
plt.savefig('example_output/gelman_rubin_plot.pdf')
plt.show()