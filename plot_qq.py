import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Read files
study_code = sys.argv[1]
file1 = sys.argv[2]
file2 = sys.argv[3]
df1 = pd.read_csv(file1, header=0)
df2 = pd.read_csv(file2, header=0)
df1['-log10(P)'] = np.negative(np.log10(df1['P']))
df2['-log10(P)'] = np.negative(np.log10(df2['P']))
# Get quantiles 1-100 for each file
np.quantile(df1['P'], 0.10)
x = np.arange(101) / 100
y1 = []
y2 = []
for quantile in x:
    y1.append(np.quantile(df1['-log10(P)'], quantile))
    y2.append(np.quantile(df2['-log10(P)'], quantile))

df = pd.DataFrame({
    'Quantile': x,
    'P - Broad': y1,
    'P - Symptom Cluster': y2
})
# Plot the two quantile scatterplots together
fig, ax = plt.subplots()
sns.regplot(x='Quantile', y='P - Broad', data=df, label='P - Broad', color='crimson', ax=ax) # lightpink, crimson
sns.regplot(x='Quantile', y='P - Symptom Cluster', data=df, label='P - Symptom Cluster', color='steelblue', ax=ax) # lightskyblue, steelblue

# Set axis labels, title, legend
ax.set_xlabel('Quantile')
ax.set_ylabel('-log10(P value)')
ax.set_title('-log10(P value) quantiles for Broad vs Symptom Cluster GWAS\nStudy: ' + study_code)
ax.legend()

# Save png to plots/{study_code}
fig.savefig('plots/' + study_code + '_qq_scatterplot.png')
