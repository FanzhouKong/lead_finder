import numpy as np
from scipy.stats import ttest_ind
from toolsets.file_io import get_list_idx
def biodata_prep(bio_data):
    pos_idx = get_list_idx(bio_data['mix'], 'pos')
    neg_idx = get_list_idx(bio_data['mix'], 'neg')
    pos_ctr = []
    neg_ctr = []
    for idx in pos_idx:
        pos_ctr.extend(bio_data.iloc[idx][2:4].values)
    for idx in neg_idx:
        neg_ctr.extend(bio_data.iloc[idx][2:4].values)
    ctr_idx = pos_idx+neg_idx
    bio_data_fractions = bio_data.drop(index=ctr_idx, axis = 0)
    bio_data_fractions.reset_index(inplace=True, drop=True)
    pos_p = []
    neg_p = []
    label =[]
    average_adjusted = []
    for index, row in bio_data_fractions.iterrows():
        pos_p_value = ttest_ind(pos_ctr, row[2:4].tolist(), alternative='less').pvalue
        pos_p.append(pos_p_value)
        neg_p_value = ttest_ind(neg_ctr, row[2:4].tolist(), alternative='greater').pvalue
        neg_p.append(neg_p_value)
        if neg_p_value <= 0.05:
            label.append(True)
        else:
            label.append(False)
        if row['Average'] >= np.median(neg_ctr):
            average_adjusted.append(0)
        else:
            average_adjusted.append(bio_data_fractions['Average'].max()-row['Average'])
    bio_data_fractions['average_adjusted']=average_adjusted
    bio_data_fractions['pos_p_value']=pos_p
    bio_data_fractions['neg_p_value']=neg_p
    bio_data_fractions['peak_label']=label
    return(bio_data_fractions)
def find_edge(lst, peak_idx):
    diff = np.diff(lst)
    for i in range(peak_idx, len(diff)):
        if diff[i] >=0:
            break
    right_edge= i
    for i in range(1, peak_idx):
        if diff[peak_idx-i] <=0:
            break
    left_edge= peak_idx-i+1
    return(left_edge, right_edge)

