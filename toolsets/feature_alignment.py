import numpy as np
from toolsets.search import quick_search_values, string_search
import pandas as pd
from tqdm import tqdm
import os
def align(file_paths, master_list_input):
    master_list = master_list_input.copy()
    isotopic_state = np.repeat('unknown', len(master_list))
    for file in tqdm(file_paths):
        file_temp = pd.read_csv(file)
        intensity_temp = []
        # suspicious_df = pd.DataFrame()
        for index, row in master_list.iterrows():
            feature_temp = find_feature(file_temp, mz = row['pmz'], rt =row['rt'])

            # feature_temp.drop_duplicates(subset = ['PeakID'])
            if len(feature_temp)>0:
                feature_temp.sort_values(by = ['Height'], inplace = True)
                intensity_temp.append(feature_temp.iloc[0]['Height'])
                if 'M + 0' in feature_temp['Isotope'].unique() and isotopic_state[index]!= 'M + 0':
                    # update isotopic state
                    isotopic_state[index]= 'M + 0'
            else:
                intensity_temp.append(0)
        master_list[os.path.basename(file).split('.')[0]]=intensity_temp
    master_list['isotopic_state']=isotopic_state
    return(master_list)

    # pass
# file_paths = [os.path.join(working_dir, x+'.csv') for x in file_list]
# alignment_result= pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
# for file in tqdm(file_paths):
#     file_temp = pd.read_csv(file)
#     intensity_temp = []
#     # pmz_temp = []
#     for i in (range(0, len(pmz_list))):
#         feature_temp = find_feature(file_temp, mz = pmz_list[i], rt =rt_list[i])
#         if len(feature_temp)>0:
#             intensity_temp.append(feature_temp['Height'].sum())
#             # pmz_temp.append()
#             # break
#         else:
#             intensity_temp.append(0)
#     alignment_result[os.path.basename(file).split('.')[0]]=intensity_temp
def initilize_pmz_rt_list(qc_paths,
                          # seed_idx
                          # istd_pmz_col = '[M+H]+',
                          # istd_rt_col = 'RT_suggested',
                          # istd_name_col = 'compound_name'
                          mode='inclusive'
                          ):

    # print('i am in current!')
    if mode == 'inclusive':
        seed_idx = determine_seed_idx(qc_paths)
        qc_seed =  pd.read_csv(qc_paths[seed_idx])
        qc_seed = drop_duplicate_features(qc_seed)
        # return(qc_seed)
        qc_rest = qc_paths[0:seed_idx]+(qc_paths[seed_idx+1:])
        pmz_list = qc_seed['Precursor m/z'].tolist()
        rt_list = qc_seed['RT_adjusted'].tolist()
        # iso_state = qc_seed['Isotope'].tolist()
        master_list = pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
        print('length_of_master list: ', str(len(master_list)))
        # return(qc_seed)
        if len(qc_rest)>0:
            for qc in qc_rest:
                qc_temp = pd.read_csv(qc)
                # break
                for index, row in tqdm(qc_temp.iterrows(), total = len(qc_temp)):
                    feature_temp = find_feature(master_list, mz = row['Precursor m/z'],
                                                rt = row['RT_adjusted'], mz_column='pmz',
                                                rt_column= 'rt'
                                                )
                    if len(feature_temp)<1:
                        pmz_list.append(row['Precursor m/z'])
                        rt_list.append(row['RT_adjusted'])

        master_list = pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
        print('length_of updated master list: ', str(len(master_list)))

    elif mode == 'exclusive':
        seed_idx = determine_seed_idx(qc_paths, mode = 'exclusive')
        qc_seed =  pd.read_csv(qc_paths[seed_idx])
        qc_seed = drop_duplicate_features(qc_seed)
        # return(qc_seed)
        qc_rest = qc_paths[0:seed_idx]+(qc_paths[seed_idx+1:])
        pmz_list = qc_seed['Precursor m/z'].tolist()
        rt_list = qc_seed['RT_adjusted'].tolist()

        master_list = pd.DataFrame(zip(pmz_list, rt_list), columns=['pmz', 'rt'])
        # return(master_list)
        print('length_of_master list: ', str(len(master_list)))
        if len(qc_rest)>0:
            for qc in qc_rest:
                qc_temp = pd.read_csv(qc)
                drop_idx = []
                for index, row in tqdm(master_list.iterrows(), total=len(master_list)):
                    # score.append(distance_score())
                    feature_temp = find_feature(qc_temp, mz = row['pmz'],
                                                rt = row['rt'], mz_column='Precursor m/z',
                                                rt_column= 'RT_adjusted'
                                                )
                    if len(feature_temp)<1:
                        drop_idx.append(index)
                master_list.drop(index = drop_idx, axis = 0, inplace=True)
        print('length_of updated master list: ', str(len(master_list)))

    else:
        print('wrong mode information is passed, accepted values: inclusive/exclusive')
        return()
    master_list.reset_index(inplace=True, drop = True)
    return master_list
def drop_duplicate_features(peak_list):
    peak_list['key']=peak_list['RT_adjusted'].astype(str)+"_"+peak_list['Precursor m/z'].astype(str)
    peak_list.drop_duplicates(subset = ['key'], inplace = True)
    peak_list.reset_index(inplace=True, drop = True)
    return(peak_list)
def determine_seed_idx(qc_paths, mode = 'inclusive'):
    feature_num = []
    for qc in qc_paths:
        temp= pd.read_csv(qc)
        feature_num.append(len(temp))
    if mode == 'inclusive':
        return(np.argmax(feature_num))
    else:
        return(np.argmin(feature_num))
def find_istd(alignmentt, istd_info, mz_column = 'pmz', rt_column = 'rt'):
    # print('i am in new')
    alignment = alignmentt.copy()
    istd_idx = []
    for index, row in istd_info.iterrows():
        # mz =
        # rt = row['RT_suggested']
        feature_mz_search = quick_search_values(alignment, mz_column, row['Precursor m/z']-0.002, row['Precursor m/z']+0.002, ifsorted=False)
        feature_mzrt_search = quick_search_values(feature_mz_search, rt_column, row['RT_suggested']-10/60, row['RT_suggested']+10/60, ifsorted=False)
        istd_idx.extend(feature_mzrt_search.index)
    # return(istd_idx)
    annotation = np.repeat('unknown', len(alignment))
    feature_type = np.repeat('compound', len(alignment))
    feature_type[istd_idx]=np.repeat('istd', len(istd_idx))
    annotation[istd_idx]=istd_info['compound_name']

    alignment.insert(2, 'annotation', annotation)
    alignment.insert(3, 'feature_type', feature_type)
    # alignment['feature_type'].loc[istd_idx]='istd'
    # alignment['annotation'].loc[istd_idx]=

    return (alignment)
def find_feature(feature_table, mz, rt, mz_column = 'Precursor m/z', rt_column = 'RT_adjusted', intensity_column =None, rt_offset = 2):
    # print('i am nin new')
    mz_step = 0.005
    # print(mz_step)
    feature_mz_search = quick_search_values(feature_table, mz_column, mz-mz_step, mz+mz_step, ifsorted=False)
    # return(feature_mz_search)
    feature_mzrt_search = quick_search_values(feature_mz_search, rt_column, rt-rt_offset/60, rt+rt_offset/60, ifsorted=False)
    if intensity_column != None:
        feature_mzrt_search.sort_values(by = intensity_column, inplace=True, ascending=False)
    # print(feature_mzrt_search)
    return (feature_mzrt_search)
def clean_bad_features(alignment_result):
    alignment_result_refined = alignment_result.copy()
    bad_idx = []
    for index, row in alignment_result.iterrows():
        if np.sum(row[4:]) == np.max(row[4:]):
            bad_idx.append(index)
    alignment_result_refined.drop(bad_idx, inplace = True)
    alignment_result_refined.reset_index(inplace = True, drop = True)
    return(alignment_result_refined)
from toolsets.file_io import get_list_idx
def filter_with_blank(alignment_result, qc_labes, blank_labels, threshold = 1):
    alignment = alignment_result.copy()
    alignment_result_qcs = alignment[qc_labes]
    alignment_result_blk = alignment[blank_labels]
    bad_feature_idx = []
    for index, row in alignment.iterrows():
        if alignment_result_qcs.loc[index].median()<alignment_result_blk.loc[index].max()*100/threshold:
            bad_feature_idx.append(index)
    alignment.drop(bad_feature_idx, inplace = True)
    return(alignment)
    # check if blank labels is given
    # if blank_lables is None:

    # blank_df = alignment_result[blank_lables]