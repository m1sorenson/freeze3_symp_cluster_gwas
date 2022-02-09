#/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df1 = pd.read_csv('results_cat/ongb_ongb_pcl4_b_5q_eur_pcs_scaled.assoc.gz',
                  header=0, sep='\t')[['ID', 'P']]
df2 = pd.read_csv('/home/maihofer/freeze3_gwas/results_cat/ongb_ongb_eur_pcs.Current_PTSD_Continuous.assoc.gz',
                  header=0, sep='\t')[['ID', 'P']]

snps1 = df1.rename({'P': 'P_sympclust'}, axis=1)
snps2 = df2.rename({'P': 'P_broad'}, axis=1)

print(snps1.head())
print(snps2.head())
print(snps1.dtypes)
print(snps2.dtypes)

snps1['P_sympclust'] = np.negative(np.log10(snps1['P_sympclust']))
snps2['P_broad'] = np.negative(np.log10(snps2['P_broad']))

snps1 = snps1.drop_duplicates(subset='ID', keep=False).set_index('ID', drop=True)
snps2 = snps2.drop_duplicates(subset='ID', keep=False).set_index('ID', drop=True)

joined = snps1.join(snps2, on='ID', how='inner')
fig, ax = plt.subplots()
sns.regplot(x='P_broad', y='P_sympclust', color='coral', scatter_kws={'s':1}, data=joined, ax=ax)
fig.savefig('plots/ongb_p_scatter.png')
