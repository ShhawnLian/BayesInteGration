# Index the results from different methods
import pandas as pd
import numpy as np
import itertools
import networkx as nx

Span_loc = 500
Ratio_len = 1.5
csv_path = 'Data/merged.csv'


def merge(intervals):
    if len(intervals) == 0 or len(intervals) == 1:
        return intervals
    intervals.sort(key=lambda x: x[0])
    result = [intervals[0]]
    for interval in intervals[1:]:
        if interval[0] <= result[-1][1]:
            result[-1][1] = max(result[-1][1], interval[1])
        else:
            result.append(interval)
    return result


def svim_dist(start1, end1, start2, end2, len1, len2, n=900):
    SD = abs(len1 - len2)/max(len1, len2)
    PD = min(abs(end1 - end2), abs(start1 - start2), abs((end1 + start1)/2 - (end2 + start2)/2))
    return SD + PD/n


dfs = pd.read_table(csv_path, sep=',')
dfs = dfs.sort_values(by=['chro', 'start', 'length'])
dfs = dfs.reset_index(drop=True)

dfs['upper_start'] = dfs['start'] + Span_loc
dfs['upper_length'] = dfs['length'] * Ratio_len

method_list = list(set(dfs['Method']))
num_method = len(method_list)

dfs_res = pd.DataFrame(columns=['ID', 'chro', 'start', 'end', 'length', 'type',
                                'method_list', 'method_number', 'score_list'])

type_list = list(set(dfs['type']))
chr_list = list(set(dfs['chro']))
for _chr in chr_list:
    for _subtype in type_list:
        sub_dfs = dfs[(dfs['type'] == _subtype) & (dfs['chro'] == _chr)]
        sub_dfs = sub_dfs.reset_index(drop=True)

        collection_var_df = pd.DataFrame(columns=['ID', 'chro', 'start', 'end', 'length', 'type',
                                                  'method_list', 'method_number', 'score_list'])

        # In fact, len(sub_dfs) should be large than 1
        if len(sub_dfs) > 0:
            # start part
            start_range = [[sub_dfs.loc[i, 'start'], sub_dfs.loc[i, 'upper_start']] for i in range(len(sub_dfs))]
            start_merge = merge(start_range)

            tmp_idx_collection = []
            for sub_list in start_merge:
                tmp_start_dfs = sub_dfs[(sub_dfs['start'] >= sub_list[0]) & (sub_dfs['upper_start'] <= sub_list[1])]

                # length part
                length_range = [[tmp_start_dfs.loc[i, 'length'], tmp_start_dfs.loc[i, 'upper_length']] for i in
                                list(tmp_start_dfs.index)]
                length_merge = merge(length_range)
                for sub_sub_list in length_merge:
                    tmp_length_dfs = tmp_start_dfs[(tmp_start_dfs['length'] >= sub_sub_list[0]) & (
                                tmp_start_dfs['upper_length'] <= sub_sub_list[1])]

                    # method part
                    order_method = list(tmp_length_dfs['Method'].value_counts().reset_index().sort_values(['count'])['Method'])
                    tmp_length_dfs.insert(tmp_length_dfs.shape[1], 'tmp_index', tmp_length_dfs.index)
                    tmp_length_dfs.index = tmp_length_dfs['Method']
                    sort_tmp_length_dfs = tmp_length_dfs.loc[order_method]
                    sort_tmp_length_dfs.index = sort_tmp_length_dfs["tmp_index"]

                    coverd_index = []
                    for idx1 in sort_tmp_length_dfs.index:
                        if idx1 not in coverd_index:
                            tmp_i_collection = [idx1]
                            tmp_method = [method for method in order_method if method != sort_tmp_length_dfs.loc[idx1, "Method"]]
                            for m in tmp_method:
                                tmp_length_method_dfs = sort_tmp_length_dfs[sort_tmp_length_dfs['Method']==m]
                                tmp_idx = [i for i in tmp_length_method_dfs.index if i not in coverd_index]
                                if (len(tmp_idx)) > 0:
                                    dist = []
                                    for idx2 in tmp_idx:
                                        dist.append(svim_dist(sort_tmp_length_dfs.loc[idx1, 'start'],
                                                              sort_tmp_length_dfs.loc[idx1, 'end'],
                                                              tmp_length_method_dfs.loc[idx2, 'start'],
                                                              tmp_length_method_dfs.loc[idx2, 'end'],
                                                              sort_tmp_length_dfs.loc[idx1, 'length'],
                                                              tmp_length_method_dfs.loc[idx2, 'length']))
                                    tmp_res = tmp_idx[dist.index(min(dist))]
                                    tmp_i_collection.append(tmp_res)
                                    coverd_index.append(tmp_res)
                            tmp_idx_collection.append(tmp_i_collection)
                            coverd_index.append(idx1)

            connected_components = tmp_idx_collection
            collection_var_df_loop = 0

            for sub_component in connected_components:
                # tmp_sub_dfs = sub_dfs.loc[sub_component, :]
                tmp_method_collection = list(sub_dfs.loc[sub_component, 'Method'])
                tmp_start_collection = list(sub_dfs.loc[sub_component, 'start'])
                tmp_end_collection = list(sub_dfs.loc[sub_component, 'end'])
                tmp_score_collection = list(sub_dfs.loc[sub_component, 'qual'])
                collection_len = list(sub_dfs.loc[sub_component, 'length'])

                mid_start = np.round(np.median(tmp_start_collection))
                mid_end = np.round(np.median(tmp_end_collection))
                if _subtype == 'INS':
                    mid_len = np.round(np.median(collection_len))
                else:
                    # For DEL, DUP, INV, take length = end - start
                    mid_len = mid_end - mid_start
                collection_var_df.loc[collection_var_df_loop, :] = \
                    [collection_var_df_loop, _chr, mid_start, mid_end, mid_len, _subtype,
                     tmp_method_collection, len(tmp_method_collection), tmp_score_collection]
                collection_var_df_loop = collection_var_df_loop + 1

        dfs_res = pd.concat([dfs_res, collection_var_df])

dfs_res = dfs_res.sort_values(by=['chro', 'start', 'length'])
dfs_res = dfs_res.reset_index(drop=True)
dfs_res['ID'] = list(range(len(dfs_res)))



# indexed matrix and score matrix
index_matrix = pd.DataFrame(np.zeros([num_method, len(dfs_res)]), index=method_list)
score_matrix = pd.DataFrame(np.full([num_method, len(dfs_res)], np.nan), index=method_list)
for location in range(len(dfs_res)):
    tmp_method_involved = dfs_res.loc[location, 'method_list']
    tmp_score_involved = dfs_res.loc[location, 'score_list']
    for m in tmp_method_involved:
        index_matrix.iloc[method_list.index(m), location] = 1
        score_matrix.iloc[method_list.index(m), location] = np.mean([tmp_score_involved[i] for i, x in enumerate(tmp_method_involved) if x == m])


index_matrix.to_csv('Indexd.csv', sep=',')
score_matrix.to_csv('Scored.csv', sep=',')

dfs_res = dfs_res.drop(columns=['method_list', 'method_number', 'score_list'])
dfs_res.to_csv('IndexSVInfo.csv', sep=',', index=False)
