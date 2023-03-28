import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr
import time

result_path=' '
res_df=pd.read_csv(result_path,index_col=0)
dfs={"sub1":res_df[res_df['sub']==1],"sub2":res_df[res_df['sub']==2]}
nrows=8
ncols=6
fig, axs = plt.subplots(nrows=nrows,ncols=ncols ,figsize=(6*nrows, 8*ncols))
i=0
row=0
# Loop over the dataframes and plot each one
for (id,df) in (dfs.items()):
    # Select the subplot to plot on
    
    if i>ncols-1:
        i=i-ncols
        row+=1
    ax = axs[row,i]
    # Plot the data using Seaborn
    sns.regplot(x="original_{}".format(task), y="predicted_{}".format(task), data=df,ax=ax).set(title='{} | {} set | total point: {}| spearman: {} | pearson: {}'.format(split, task, df.shape[0], np.round(metrics_dict[f'{id}-spearman'], 3), np.round(metrics_dict[f'{id}-pearson'], 3) ))
    i+=1
plt.savefig(f'./figures/fold_{fold}/all.png', bbox_inches='tight')