import toolsets.raw_data_scaffold as rds
import numpy as np
import pandas as pd
from toolsets.search import string_search, quick_search_sorted, quick_search_values
import toolsets.helpers as helpers
import pybaselines
from tqdm import tqdm
import os
import toolsets.helpers as helpers
import toolsets.spectra_operations as so
from collections import Counter

def find_features(mzml_path, parent_dir,intensity_threshold = 30000, istd_info=None,
                  n_neighbor = 2):
    ms1, ms2 = rds.read_mzml(mzml_path, parent_dir = parent_dir)
    mass_sorted, intensity_sorted, index_sorted, rt_list = build_index(ms1)
    features_df = get_features(mass_sorted, intensity_sorted, index_sorted, rt_list, base_name=ms1.iloc[0]['base_name'], intensity_threshold=intensity_threshold,
                                n_neighbor = n_neighbor)
    iso_state_lst = [None]*len(features_df)
    if_dut_lst = [None]*len(features_df)
    if_halo_lst = [None]*len(features_df)
    idx = 0
    # ms1s = [None]*len(features_df)
    # for index, row in (features_df.iterrows()):
    #     iso_state, if_dut, if_halo = _determine_iso_state(row, mass_sorted, intensity_sorted, index_sorted)
    #     iso_state_lst[idx]=iso_state
    #     if_dut_lst[idx]=if_dut
    #     if_halo_lst[idx]=if_halo
    #     mass = ms1.iloc[row['ms1_idx']]['peaks'].T[0]
    #     inte = ms1.iloc[row['ms1_idx']]['peaks'].T[1]
    #
    #     iso_idx_start, iso_idx_end = mass.searchsorted([(row['precursor_mz'])-0.005,row['precursor_mz']+3.2])
    #     mass = mass[iso_idx_start:iso_idx_end]
    #     inte = inte[iso_idx_start:iso_idx_end]
    #     inte = [x/np.max(inte) for x in inte]
    #     ms1s[idx]=so.pack_spectra(mass, inte)
    #     idx = idx+1
    # features_df['iso_state']=iso_state_lst
    # features_df['if_dut']=if_dut_lst
    # features_df['if_halo']= if_halo_lst
    # features_df['ms1']=ms1s
    # if istd_info is not None:
    #     features_df = find_istd_info(istd_info, features_df)
    # features_df_with_ms2, unmapped_ms2 = map_ms2(features_df, ms2)
    return(features_df)
def get_features(mass_sorted, intensity_sorted, index_sorted, rt_list, base_name=None, intensity_threshold=30000, n_neighbor = 2):
    pmz = np.zeros(len(mass_sorted))
    ion_trace_center = np.zeros(len(mass_sorted))
    # ion_trace_offset = np.zeros(len(mass_sorted))
    rt = np.zeros(len(mass_sorted))
    rt_start =np.zeros(len(mass_sorted))
    rt_end =np.zeros(len(mass_sorted))
    apex_intensity = np.zeros(len(mass_sorted))
    n_scans = np.zeros(len(mass_sorted))
    peak_range = [None]*(len(mass_sorted))
    ms1_idx = [None]*(len(mass_sorted))
    reci_snr_all = np.zeros((len(mass_sorted)))
    idx = 0
    intensity_sorted_indexing = intensity_sorted.copy()
    # num_grater_than_zero = np.sum(intensity_sorted_indexing>0)
    while np.max(intensity_sorted_indexing)>30000:

        seed_idx = np.argmax(intensity_sorted_indexing)
        seed_mass = mass_sorted[seed_idx]
        ion_trace = flash_eic(seed_mass, mass_sorted, intensity_sorted, index_sorted)
        peaks_all, reci_snrs_all, raw_apex_idx_all = detect_all_peaks(ion_trace,
                                                                      n_neighbor=n_neighbor,
                                                                      intensity_threshold=intensity_threshold)
        idx_start, idx_end = mass_sorted.searchsorted([seed_mass-0.005, seed_mass+0.005])
        intensity_sorted_indexing[idx_start:idx_end]=0
        # print(np.max(intensity_sorted_indexing))
        # if np.sum(intensity_sorted_indexing>0)==num_grater_than_zero:
        #     print('not zeroing out')
        #     return(seed_mass)
        # else:
        #     num_grater_than_zero = np.sum(intensity_sorted_indexing>0)
        if len(peaks_all)>0:
            for p, r, a in zip(peaks_all, reci_snrs_all, raw_apex_idx_all):
                pmz_statistics = guess_pmz(seed_mass, mass_sorted,
                                           intensity_sorted, index_sorted, idx_start, idx_end, int(a))
                if pmz_statistics[0]==pmz_statistics[0]:
                    pmz[idx]=pmz_statistics[0]
                    apex_intensity[idx]=pmz_statistics[1]
                    ms1_idx[idx]=a
                    reci_snr_all[idx]=r
                    ion_trace_center[idx]=seed_mass
                    rt_start[idx]=rt_list[p[0]]
                    rt_end[idx]=rt_list[p[2]]
                    n_scans[idx]=p[2]-p[0]-1
                    peak_range[idx]=[p[0],a,p[2]]
                    try:
                        rt[idx]=gaussian_estimator(tuple([int(a)-1,int(a), int(a)+1]),rt_list, ion_trace)
                        if rt[idx]!=rt[idx]:
                            rt[idx]=rt_list[int(a)]
                    except:
                        rt[idx]=rt_list[int(a)]
                    idx = idx+1


        # print(np.max(intensity_sorted_indexing))

    pmz = pmz[0:idx]
    rt = rt[0:idx]
    rt_start = rt_start[0:idx]
    rt_end = rt_end[0:idx]
    apex_intensity = apex_intensity[0:idx]
    n_scans = n_scans[0:idx]
    peak_range = peak_range[0:idx]
    ion_trace_center = ion_trace_center[0:idx]
    reci_snr_all=reci_snr_all[0:idx]
    df = pd.DataFrame(zip(pmz, rt, rt_start, rt_end,
                          apex_intensity,
                          n_scans,peak_range,ion_trace_center,reci_snr_all),
                      columns=['precursor_mz','rt_apex', 'rt_start', 'rt_end',
                               'ms1_intensity',
                               'n_scnas', 'ms1_scan_range','ion_trace_center', 'reci_snr'])
    df['base_name']=base_name
    return(df)
def detect_all_peaks(intensity_list, intensity_threshold = 30000, n_neighbor=2, return_all = False):
    intensity_list_smoothed = np.array(moving_average(intensity_list, n= n_neighbor))
    peaks_all =get_peaks(intensity_list_smoothed)
    peaks_return = [None]*len(peaks_all)
    reci_snrs = np.zeros(len(peaks_all))
    n_scans = np.zeros(len(peaks_all))
    # smoothed_apex_intensity =np.zeros(len(peaks_all))
    idx = 0
    raw_apex_idx = np.zeros(len(peaks_all), dtype=int)
    all_apex_intensity = np.array([intensity_list_smoothed[x[1]]for x in peaks_all])
    while np.max(all_apex_intensity)>intensity_threshold:
        current_peak_idx = np.argmax(all_apex_intensity)
        apex_peak = np.array(peaks_all[current_peak_idx])
        apex_peak_range = [current_peak_idx,current_peak_idx]
        l = 1
        r = 1
        while current_peak_idx-l>0:
            left_peak = peaks_all[current_peak_idx-l]
            if apex_peak[0]-left_peak[2]<=2 and intensity_list_smoothed[apex_peak[0]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[left_peak[1]])/1.3:
                apex_peak[0]=left_peak[0]
                apex_peak_range[0]=current_peak_idx-l
                l = l+1
            else:
                break
        while current_peak_idx+r<len(peaks_all):
            right_peak = peaks_all[current_peak_idx+r]
            if right_peak[0]-apex_peak[2]<=2 and intensity_list_smoothed[apex_peak[2]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[right_peak[1]])/1.3:
                apex_peak[2]=right_peak[2]
                apex_peak_range[1]=current_peak_idx+r
                r = r+1
            else:
                break

        # max_half_window = int(np.ceil((apex_peak[2]-apex_peak[0])/2))
        # baseline = pybaselines.smooth.snip(
        #     intensity_list_smoothed,
        #     max_half_window=40,
        #     decreasing=True)[0]
        # baseline = [0 if x <0 else x for x in baseline]
        # reci_snr = np.median(baseline[apex_peak[0]:apex_peak[2]+1])/intensity_list_smoothed[apex_peak[1]]

        peaks_return[idx]=(apex_peak)
        n_scans[idx]=apex_peak[2]-apex_peak[0]
        # smoothed_apex_intensity[idx]=intensity_list_smoothed[apex_peak[1]]
        raw_apex_idx[idx]=(apex_peak[0]+np.argmax(intensity_list[apex_peak[0]:apex_peak[2]+1]))
        idx = idx+1
        # all_apex_intensity = np.concatenate((all_apex_intensity[0:apex_peak_range[0]],all_apex_intensity[apex_peak_range[1]+1:]), axis=None)
        peaks_all = peaks_all[0:apex_peak_range[0]]+peaks_all[apex_peak_range[1]+1:]
        all_apex_intensity = np.array([intensity_list_smoothed[x[1]]for x in peaks_all])

        if len(peaks_all)==0:
            break
    # n_scans=n_scans[0:idx]

    max_half_window = int(np.ceil(np.max(n_scans)/2))
    # print(max_half_window)
    baseline = pybaselines.smooth.snip(
        intensity_list_smoothed,
        max_half_window=max_half_window,
        decreasing=False)[0]
    baseline = np.array([0 if x <0 else x for x in baseline])
    for i in range(0, idx):
        reci_snrs[i]=np.median(baseline[peaks_return[i][0]:peaks_return[i][2]+1])/intensity_list_smoothed[peaks_return[i][1]]
    peaks_return=peaks_return[0:idx]
    reci_snrs = reci_snrs[0:idx]
    raw_apex_idx=raw_apex_idx[0:idx]
    if return_all==True:
        return (peaks_return, reci_snrs, raw_apex_idx)



    high_quality_idx = reci_snrs<1/3
    peaks_return=np.array(peaks_return)

    peaks_return = peaks_return[high_quality_idx]
    reci_snrs=reci_snrs[high_quality_idx]
    raw_apex_idx = raw_apex_idx[high_quality_idx]
    return(peaks_return,reci_snrs , raw_apex_idx)
def guess_pmz(target_mass,mass_sorted, intensity_sorted,index_sorted,idx_start, idx_end, peak_apex, guess_step = 1):
    # idx_start, idx_end = mass_sorted.searchsorted([seed_mass-mass_tolerance,seed_mass+mass_tolerance ])
    mass_range = mass_sorted[idx_start:idx_end]
    index_range = index_sorted[idx_start:idx_end]
    intensity_range = intensity_sorted[idx_start:idx_end]
    if np.max(intensity_range)==0:
        return(np.NAN, np.NAN)
    pmz_candidates = np.zeros(2*guess_step+1)
    intensity_candidates = np.zeros(2*guess_step+1)
    pmz_idx = 0
    for i in range(peak_apex-guess_step, peak_apex+guess_step+1):
        if len(intensity_range[index_range==i])>0:
            difference_array = np.absolute(mass_range[index_range==i]-target_mass)

            mass_anchor = difference_array.argmin()
            pmz_candidates[pmz_idx]=mass_range[index_range==i][mass_anchor]
            intensity_candidates[pmz_idx]=intensity_range[index_range==i][mass_anchor]

        # if i == peak_apex:
        #     apex_intensity = intensity_candidates[pmz_idx]
        pmz_idx = pmz_idx+1
    pmz_candidates=pmz_candidates[0:pmz_idx]
    # print(pmz_candidates, peak_apex)
    intensity_candidates=intensity_candidates[0:pmz_idx]
    weighted_pmz = np.sum([x*y/np.sum(intensity_candidates) for x, y in zip(pmz_candidates, intensity_candidates)])
    apex_intensity = flash_eic(weighted_pmz, mass_sorted, intensity_sorted, index_sorted)[peak_apex]
    return (weighted_pmz, apex_intensity)
# def get_features_alt(mass_sorted, intensity_sorted, index_sorted, rt_list, base_name=None, intensity_threshold=30000, n_neighbor = 2):
#     pmz = np.zeros(len(mass_sorted))
#     print('updated!')
#     ion_trace_center = np.zeros(len(mass_sorted))
#     ion_trace_offset = np.zeros(len(mass_sorted))
#     rt = np.zeros(len(mass_sorted))
#     rt_start =np.zeros(len(mass_sorted))
#     rt_end =np.zeros(len(mass_sorted))
#     apex_intensity_raw = np.zeros(len(mass_sorted))
#     n_scans = np.zeros(len(mass_sorted))
#     peak_range = [None]*(len(mass_sorted))
#     reci_snr_all = np.zeros((len(mass_sorted)))
#     idx = 0
#     intensity_sorted_indexing = intensity_sorted.copy()
#     num_grater_than_zero=np.sum(intensity_sorted_indexing>30000)
#     while np.max(intensity_sorted_indexing)>intensity_threshold:
#         # print(np.max(intensity_sorted_indexing))
#         seed_idx = np.argmax(intensity_sorted_indexing)
#         seed_mass = mass_sorted[seed_idx]
#         idx_start, idx_end = mass_sorted.searchsorted([seed_mass-0.005,seed_mass+0.005 ])
#         mass_offset = np.array([abs(x-seed_mass) for x in mass_sorted[idx_start:idx_end]])
#         mass_tol = np.mean(mass_offset)+3*np.std(mass_offset)
#         if mass_tol <0.005:# this is a clean ion trace
#             intensity_list = flash_eic(seed_mass, mass_sorted, intensity_sorted, index_sorted,
#                                        mass_error=mass_tol)
#             idx_start, idx_end = mass_sorted.searchsorted([seed_mass-mass_tol,seed_mass+mass_tol ])
#             peaks_all, reci_snrs, apex_intensity_smoothed_all= detect_all_peaks(intensity_list,
#                                                                                 n_neighbor=n_neighbor,
#                                                                                 intensity_threshold=intensity_threshold)
#             for p, r in zip(peaks_all,reci_snrs):
#                 apex_range_left = p[0]
#                 apex_range_right = p[2]
#                 apex_offset = np.argmax(intensity_list[apex_range_left:apex_range_right+1])
#                 apex_index = apex_range_left+apex_offset
#                 try:
#                     pmz_statistics = guess_pmz(mass_sorted, intensity_sorted, index_sorted,
#                                                idx_start, idx_end, apex_index)
#                     pmz[idx]=pmz_statistics[0]
#                     rt[idx]=gaussian_estimator(tuple([apex_index-1,apex_index, apex_index+1]),rt_list, intensity_list)
#                     apex_intensity_raw[idx]=pmz_statistics[1]
#                 except:
#                     mass_range = mass_sorted[idx_start:idx_end]
#                     intensity_range = intensity_sorted[idx_start:idx_end]
#                     index_range = index_sorted[idx_start:idx_end]
#                     intensity_anchor = np.argmax(intensity_range[index_range==apex_index])
#                     pmz[idx]=mass_range[index_range==apex_index][intensity_anchor]
#                     rt[idx]=rt_list[apex_index]
#                     apex_intensity_raw[idx]=np.max(intensity_list[apex_range_left:apex_range_right+1])
#                 peak_range[idx]= [p[0], apex_index, p[2]]
#                 ion_trace_center[idx]=seed_mass
#                 ion_trace_offset[idx]=mass_tol
#                 rt_start[idx]=rt_list[p[0]]
#                 rt_end[idx]=rt_list[p[2]]
#                 n_scans[idx]=p[2]-p[0]-1
#                 reci_snr_all[idx]=r
#                 idx = idx+1
#             intensity_sorted_indexing[idx_start:idx_end]=0#update index intensity
#         else:
#             intensity_list = flash_eic(seed_mass, mass_sorted, intensity_sorted, index_sorted,
#                                        mass_error=0.005)
#             peak, reci_snr, apex_intensity_smoothed = detect_one_peak(intensity_list,
#                                                                       target_idx=index_sorted[seed_idx],
#                                                                       n_neighbor=10)
#             if peak != []:
#                 apex_range_left = p[0]
#                 apex_range_right = p[2]
#                 apex_index = index_sorted[seed_idx]
#                 try:
#                     pmz_statistics = guess_pmz(mass_sorted, intensity_sorted, index_sorted, idx_start, idx_end, apex_index)
#                     pmz[idx]=pmz_statistics[0]
#                     rt[idx]=gaussian_estimator(tuple([apex_index-1,apex_index, apex_index+1]),rt_list, intensity_list)
#                     apex_intensity_raw[idx]=pmz_statistics[1]
#                 except:
#                     mass_range = mass_sorted[idx_start:idx_end]
#                     intensity_range = intensity_sorted[idx_start:idx_end]
#                     index_range = index_sorted[idx_start:idx_end]
#                     intensity_anchor = np.argmax(intensity_range[index_range==apex_index])
#                     pmz[idx]=mass_range[index_range==apex_index][intensity_anchor]
#                     rt[idx]=rt_list[apex_index]
#                     apex_intensity_raw[idx]=np.max(intensity_list[apex_range_left:apex_range_right+1])
#                 idx = idx+1
#                 #update intensity index
#                 intensity_sorted_indexing[idx_start:idx_end][np.logical_and(index_sorted[idx_start:idx_end]>=peak[0],
#                                                                             index_sorted[idx_start:idx_end]<=peak[2])]=0
#             else:
#                 intensity_sorted_indexing[np.argmax[intensity_sorted_indexing]]=0
#         if np.sum(intensity_sorted_indexing>30000)==num_grater_than_zero:
#             print('not zeroing out')
#             return(seed_mass, mass_tol)
#         else:
#             num_grater_than_zero = np.sum(intensity_sorted_indexing>30000)
#     pmz = pmz[0:idx]
#     rt = rt[0:idx]
#     rt_start =  rt_start[0:idx]
#     rt_end = rt_end[0:idx]
#     apex_intensity_raw = apex_intensity_raw[0:idx]
#     n_scans = n_scans[0:idx]
#     peak_range = peak_range[0:idx]
#     ion_trace_center = ion_trace_center[0:idx]
#     reci_snr_all=reci_snr_all[0:idx]
#     ion_trace_offset= ion_trace_offset[0:idx]
#     df = pd.DataFrame(zip(pmz, rt, rt_start, rt_end,
#                           apex_intensity_raw,
#                           n_scans,peak_range,ion_trace_center,reci_snr_all, ion_trace_offset),
#                       columns=['precursor_mz','rt_apex', 'rt_start', 'rt_end',
#                                'ms1_intensity',
#                                'n_scnas', 'ms1_scan_range','ion_trace_center', 'reci_snr', 'ion_trace_offset'])
#     df['base_name']=base_name
#     return(df)
from scipy.signal import find_peaks
# def detect_one_peak(intensity_list, target_idx, n_neighbor = 2):
#     intensity_list_smoothed = np.array(moving_average(intensity_list, n= n_neighbor))
#     peaks_all =get_peaks(intensity_list_smoothed)
#
#     offsets = [abs(target_idx-x[1]) for x in peaks_all]
#     current_peak_idx = np.argmin(offsets)
#     apex_peak = np.array(peaks_all[current_peak_idx])
#     apex_peak_range = [current_peak_idx,current_peak_idx]
#     l = 1
#     r = 1
#     while current_peak_idx-l>0:
#         left_peak = peaks_all[current_peak_idx-l]
#         if apex_peak[0]-left_peak[2]<=2 and intensity_list_smoothed[apex_peak[0]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[left_peak[1]])/1.3:
#             apex_peak[0]=left_peak[0]
#             apex_peak_range[0]=current_peak_idx-l
#             l = l+1
#         else:
#             break
#     while current_peak_idx+r<len(peaks_all):
#         right_peak = peaks_all[current_peak_idx+r]
#         if right_peak[0]-apex_peak[2]<=2 and intensity_list_smoothed[apex_peak[2]]>min(intensity_list_smoothed[apex_peak[1]], intensity_list_smoothed[right_peak[1]])/1.3:
#             apex_peak[2]=right_peak[2]
#             apex_peak_range[1]=current_peak_idx+r
#             r = r+1
#         else:
#             break
#     max_half_window = int(np.ceil((apex_peak[2]-apex_peak[0])/2))
#     baseline = pybaselines.smooth.snip(
#         intensity_list_smoothed,
#         max_half_window=max_half_window,
#         decreasing=True)[0]
#     baseline = [0 if x <0 else x for x in baseline]
#     reci_snr = np.median(baseline[apex_peak[0]:apex_peak[2]+1])/intensity_list_smoothed[apex_peak[1]]
#     return(apex_peak, reci_snr, np.max(intensity_list_smoothed[apex_peak[0]:apex_peak[2]+1]))

def map_ms2(features_df, ms2):
    msms = []
    ms2_scan_idx = []
    ms2_working = ms2.copy()
    features_df.sort_values(by = 'ms1_intensity', ascending=False, inplace=True)
    ms2_working.sort_values(by = 'precursor_mz', ascending=True, inplace=True)
    for index, row in features_df.iterrows():
        if len(ms2_working)==0:
            break
        pmz_matched = quick_search_sorted(ms2_working, 'precursor_mz', row['precursor_mz']-0.005, row['precursor_mz']+0.005)
        ms2_under_peak = quick_search_values(pmz_matched, 'ms1_rt', row['rt_start'], row['rt_end'])
        ms2_under_peak.sort_values(by = 'ms1_precursor_intensity', ascending=False, inplace = True)
        if len(ms2_under_peak)>0:

            msms.append(so.convert_array_to_string(ms2_under_peak.iloc[0]['peaks']))
            ms2_scan_idx.append(ms2_under_peak.iloc[0]['scan_idx'])
            ms2_working.drop(ms2_under_peak.index.tolist(), inplace=True)
        else:
            msms.append(np.NAN)
            ms2_scan_idx.append(np.NAN)
    features_df['msms']=msms
    features_df['ms2_scan_idx']=ms2_scan_idx
    return(features_df,ms2_working)


import itertools
def build_index(ms1):
    # ms1.reset_index(inplace = True, drop = True)
    mass_nested = [None]*len(ms1)
    intensity_nested = [None]*len(ms1)
    index_nested = [None]*len(ms1)
    rt_list = np.zeros(len(ms1))
    for index, row in (ms1.iterrows()):
        mass_temp, intensity_temp = row['peaks'].T
        mass_nested[index]=(mass_temp)
        intensity_nested[index]=(intensity_temp)
        index_nested[index]=([index]*len(mass_temp))
        rt_list[index]=(row['rt'])
    mass_flatten = np.array(list(itertools.chain.from_iterable(mass_nested)))
    intensity_flatten = np.array(list(itertools.chain.from_iterable(intensity_nested)))
    index_flatten = np.array(list(itertools.chain.from_iterable(index_nested)))
    order = np.argsort(mass_flatten)
    mass_sorted = mass_flatten[order]
    intensity_sorted = intensity_flatten[order]
    index_sorted = index_flatten[order]
    return(mass_sorted, intensity_sorted, index_sorted, rt_list)
def flash_eic(pmz, mass_sorted, intensity_sorted, index_sorted, mass_error=0.005):
    index_start, index_end = mass_sorted.searchsorted([pmz-mass_error, pmz+mass_error+1E-9])
    index_range = index_sorted[index_start:index_end]

    intensity_range = intensity_sorted[index_start:index_end]

    intensity_list = np.zeros(np.max(index_sorted)+1)
    for idx in range(0,len(index_range)):
        intensity_list[index_range[idx]]= intensity_list[index_range[idx]]+intensity_range[idx]
    return(intensity_list)
def moving_average( x, n=2):
    window = np.ones(2*n + 1) / (2*n + 1)

    # Use 'same' mode in convolve to ensure the output array has the same size as the input array
    # 'valid' mode would reduce the size of the array to avoid edge effects
    # 'full' mode would increase the size to include all possible overlaps
    smoothed_arr = np.convolve(x, window, mode='same')

    return smoothed_arr



def get_peaks(intensity_list: np.ndarray) -> list:
    """Detects peaks in an array.

    Args:
        int_array (np.ndarray): An array with intensity values.

    Returns:
        list: A regular Python list with all peaks.
            A peak is a triplet of the form (start, center, end)

    """
    apex, _ = find_peaks(intensity_list)
    peak_list = []
    for cur_apex_idx in apex:
        peak_list.append(get_edges(intensity_list, cur_apex_idx))
    return(peak_list)
def get_edges(intensity_list, cur_apex_idx):
    intensity_list = np.array(intensity_list)
    # gradient = np.diff(intensity_list)
    left_edge_idx = cur_apex_idx-1
    right_edge_idx = cur_apex_idx+1
    while left_edge_idx>0:
        if intensity_list[left_edge_idx-1]<=intensity_list[left_edge_idx] and intensity_list[left_edge_idx]>0:
            left_edge_idx = left_edge_idx-1
        else:
            break
    while right_edge_idx <len(intensity_list)-1:
        if intensity_list[right_edge_idx+1]<=intensity_list[right_edge_idx] and intensity_list[right_edge_idx]>0:
            right_edge_idx = right_edge_idx+1
        else:
            break

    return((left_edge_idx, cur_apex_idx, right_edge_idx))
#%%
def gaussian_estimator(
        peak: tuple,
        mz_array: np.ndarray,
        int_array: np.ndarray
) -> float:
    """Three-point gaussian estimator.

    Args:
        peak (tuple): A triplet of the form (start, center, end)
        mz_array (np.ndarray): An array with mz values.
        int_array (np.ndarray): An array with intensity values.

    Returns:
        float: The gaussian estimate of the center.
    """
    start, center, end = peak

    m1, m2, m3 = mz_array[center - 1], mz_array[center], mz_array[center + 1]
    i1, i2, i3 = int_array[center - 1], int_array[center], int_array[center + 1]
    # print(m1,m2,m3)
    if i1 == 0:  # Case of sharp flanks
        m = (m2 * i2 + m3 * i3) / (i2 + i3)
    elif i3 == 0:
        m = (m1 * i1 + m2 * i2) / (i1 + i2)
    else:
        l1, l2, l3 = np.log(i1), np.log(i2), np.log(i3)
        m = (
                ((l2 - l3) * (m1 ** 2) + (l3 - l1) * (m2 ** 2) + (l1 - l2) * (m3 ** 2))
                / ((l2 - l3) * (m1) + (l3 - l1) * (m2) + (l1 - l2) * (m3))
                * 1
                / 2
        )

    return m
# def find_roi(ion_trace):
#     apex_index = np.argmax(ion_trace)
#     ion_trace_smoothed = trx.moving_average(ion_trace, 10)
#     if np.min(ion_trace_smoothed)>0:
#         return([0, apex_index, len(ion_trace)-1])
#     apex_smoothed = np.argmax(ion_trace)
#     left_side = ion_trace_smoothed[0:apex_smoothed]
#     right_side = ion_trace_smoothed[apex_smoothed+1:]
#     if np.min(left_side)>0:
#         left_bound = 0
#     else:
#         # for i in range(1, apex_smoothed+1):
#         #     if ion_trace_smoothed[apex_smoothed-i]==0:
#         #         left_bound= apex_smoothed-i
#         #         break
#         zero_position_left = np.where(left_side == 0)[0]
#         # index = (apex_smoothed-zero_position_left).argmin()
#         left_bound = zero_position_left[-1]
#     if np.max(right_side)>0:
#         right_bound = len(ion_trace)-1
#     else:
#         # for l in range(1, len(apex_smoothed)-apex_index):
#         #     if ion_trace_smoothed[apex_smoothed+l]==0:
#         #         right_bound = apex_smoothed+l
#         #         break
#         zero_position_right = np.where(right_side == 0)[0]
#         # index = (zero_position_right-apex_smoothed).argmin()
#         right_bound = zero_position_right[0]
#     return(np.array([left_bound, apex_index, right_bound]))
