import numpy as np
from toolsets.search import quick_search_values, string_search
import pandas as pd
from tqdm import tqdm
import os
def initilize_pmz_rt_list(qc_paths,
                          # seed_idx
                          # istd_pmz_col = '[M+H]+',
                          # istd_rt_col = 'RT_suggested',
                          # istd_name_col = 'compound_name'
                          ):

    # print('i am in current!')
    seed_idx = determine_seed_idx(qc_paths)
    print(seed_idx)
    qc_seed =  pd.read_csv(qc_paths[seed_idx])
    qc_rest = qc_paths[0:seed_idx]+(qc_paths[seed_idx+1:])
    pmz_list = qc_seed['Precursor m/z'].tolist()
    rt_list = qc_seed['RT_adjusted'].tolist()
    iso_state = qc_seed['Isotope'].tolist()
    master_list = pd.DataFrame(zip(pmz_list, rt_list, iso_state), columns=['pmz', 'rt', 'iso_state'])
    # return(qc_seed)
    if len(qc_rest)>0:
        for qc in qc_rest:
            qc_temp = pd.read_csv(qc)
            # break
            for index, row in master_list.iterrows():
                pass
            # for index, row in tqdm(qc_temp.iterrows(), total = len(qc_temp)):
            #     feature_temp = find_feature(master_list, mz = row['Precursor m/z'],
            #                                 rt = row['RT_adjusted'], mz_column='pmz',
            #                                 rt_column= 'rt'
            #                                 )
            #     if len(feature_temp)<1:
            #         pmz_list.append(row['Precursor m/z'])
            #         rt_list.append(row['RT_adjusted'])
            #         iso_state.append(row['Isotope'])
            #     else:
            #         if
            # master_list = pd.DataFrame(zip(pmz_list, rt_list, iso_state), columns=['pmz', 'rt', 'isostate'])

    return(master_list)
def determine_seed_idx(qc_paths, mode = 'inclusive'):
    feature_num = []
    for qc in qc_paths:
        temp= pd.read_csv(qc)
        feature_num.append(len(temp))
    return(np.argmax(feature_num))
    # pass
    # qc_seed =  pd.read_csv(qc_path[seed_idx])
    # pmz_list = qc_seed['Precursor m/z'].tolist()
    # rt_list = qc_seed['RT_adjusted'].tolist()
    # # qc_list_excluded = qc_path[0:seed_idx]
    # for qc in qc_path[]:
    #     qc_temp = pd.read_csv(qc)
    #     # break
    #     for index, row in tqdm(qc_temp.iterrows(), total = len(qc_temp)):
    #         feature_temp = find_feature(qc_seed, mz = row['Precursor m/z'],
    #                                     rt = row['RT_adjusted'], mz_column='Precursor m/z',
    #                                     rt_column= 'RT_adjusted'
    #                                     )
    #         if len(feature_temp)<1:
    #             pmz_list.append(row['Precursor m/z'])
    #             rt_list.append(row['RT_adjusted'])
    # ini_peak_list = pd.read_csv(seed_path)
    # # return(ini_peak_list)
    # master_list = ini_peak_list[['RT_adjusted', 'Precursor m/z']]
    #
    # master_list[os.path.basename(seed_path).split('.')[0]]=ini_peak_list['Height']
    # comments = np.repeat(np.NAN, len(master_list))
    # annotation = np.repeat('unknown', len(master_list))
    # feature_type = np.repeat('compound', len(master_list))
    # # return(master_list)
    # master_list.insert(2, 'comments', comments)
    # master_list.insert(2, 'feature_type', feature_type)
    # master_list.insert(4, 'annotation', annotation)
    #
    # istd_idx = []
    # for index, row in istd_info.iterrows():
    #     istd_idx.extend(find_istd(master_list, mz = row['[M+H]+'], rt = row['RT_suggested']).index)
    # # return(master_list)
    # master_list['feature_type'][istd_idx]='istd'
    # master_list['annotation'][istd_idx]=istd_info['compound_name']
    # return(master_list)
def find_istd(alignment, istd_info, mz_column = 'pmz', rt_column = 'rt'):
    # print('i am in new')
    istd_idx = []
    for index, row in istd_info.iterrows():
        # mz =
        # rt = row['RT_suggested']
        feature_mz_search = quick_search_values(alignment, mz_column, row['[M+H]+']-0.002, row['[M+H]+']+0.002, ifsorted=False)
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
def find_feature(feature_table, mz, rt, mz_column = 'Precursor m/z', rt_column = 'RT_adjusted', intensity_column =None):
    # print('i am nin new')
    mz_step = 0.005
    # print(mz_step)
    feature_mz_search = quick_search_values(feature_table, mz_column, mz-mz_step, mz+mz_step, ifsorted=False)
    # return(feature_mz_search)
    feature_mzrt_search = quick_search_values(feature_mz_search, rt_column, rt-2/60, rt+2/60, ifsorted=False)
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