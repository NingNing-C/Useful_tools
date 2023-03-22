import numpy as np
import random
import pandas as pd

def divide_clusters_into_folds(dataset, clusters, fold_size):
    # Step 1: Get list of data points for each cluster
    cluster_data = {}
    for i, cluster in enumerate(clusters):
        if cluster not in cluster_data:
            cluster_data[cluster] = []
        cluster_data[cluster].append(dataset[i])

    # Step 2: Compute size of each cluster's data list
    cluster_sizes = {}
    for cluster in cluster_data:
        cluster_sizes[cluster] = len(cluster_data[cluster])

    # Step 3: Sort clusters in descending order of data list size
    # sorted_clusters = sorted(cluster_sizes, key=cluster_sizes.get, reverse=True)
    random.seed(3)
    sorted_clusters =list(cluster_sizes.keys())
    random.shuffle(sorted_clusters)
    # Step 4: Initialize empty folds dictionary
    folds = {}

    # Step 5: Compute average data points per fold
    avg_data_points = len(dataset) // fold_size

    # Step 6: Initialize list of used clusters
    used_clusters = []

    # Step 7: Initialize fold ID
    fold_id = 0

    # Step 8: Loop through sorted clusters
    for i in sorted_clusters:
        # Step 9: If enough folds have been created, exit loop
        if len(folds) == fold_size:
            break

        # Step 10: Compute threshold for current fold
        threshold = avg_data_points + 50

        # Step 11: Loop through remaining clusters
        for j in sorted_clusters:
            # Step 12: If cluster has already been used, skip it
            if j in used_clusters:
                continue

            # Step 13: If cluster i is larger than cluster j, skip it
            if cluster_sizes[j] >= cluster_sizes[i]:
                continue

            # Step 14: If adding cluster j to current fold won't exceed threshold, add it
            if len(folds.get(fold_id, [])) + cluster_sizes[j] < threshold:
                folds.setdefault(fold_id, []).extend(cluster_data[j])
                used_clusters.append(j)
            # Step 15: Otherwise, move on to next fold
            else:
                fold_id += 1
                break

    # Step 19: Add any unused clusters to last fold
    for j in sorted_clusters:
        if j not in used_clusters:
            folds.setdefault(-1, []).extend(cluster_data[j])

    # Step 22: Return folds dictionary
    return folds

def divide_clusters_into_folds(dataset, clusters, N):
    
    cluster_data = {}
    for i, cluster in enumerate(clusters):
        if cluster not in cluster_data:
            cluster_data[cluster] = []
        cluster_data[cluster].append(dataset[i])
    
    cluster_sizes = {}
    for cluster in cluster_data:
        cluster_sizes[cluster] = len(cluster_data[cluster])
    
    unused_clusters =list(cluster_data.keys())
    random.shuffle(unused_clusters)
    folds = {}    
    avg_data_points = len(dataset) // N    
    used_clusters = []
        
    for fold_id in range(N):
        print(fold_id)
        threshold = avg_data_points + np.random.randint(0,20)
        while(len(folds.get(fold_id, []))<threshold)&(len(unused_clusters)>0):
            seen_clusters=[]
            cluster_id=unused_clusters.pop()
            seen_clusters.append(cluster_id)
            if len(folds.get(fold_id, [])) + cluster_sizes[cluster_id] < threshold:
                folds.setdefault(fold_id, []).extend(cluster_data[cluster_id])
            else:
                unused_clusters.append(cluster_id)
                break

    smallest_clsuter = min(folds, key=lambda k: len(folds[k]))
    for j in unused_clusters:
        folds.setdefault(smallest_clsuter, []).extend(cluster_data[j])
    return folds


## Test codes
# df=pd.read_csv('./data/binding_affinity/SKP1402m.csv')
# fold_dir=divide_clusters_into_folds(dataset, clusters, 5)

# sum=0
# for k,v in fold_dir.items():
#     print(k, len(v))
#     sum+=len(v)

# for i in fold_dir.keys():
#     print(i)
#     print(df.iloc[fold_dir[i]].cluster_pair.value_counts())
