from feature_finding.load_mzml_data import load_mzml_data
import numpy as np
import pandas as pd
from alphapept.feature_finding import extract_hills, remove_duplicate_hills, split_hills, filter_hills, get_hill_data, get_pre_isotope_patterns, feature_finder_report, get_isotope_patterns
from alphapept.constants import averagine_aa, isotopes, mass_dict
from toolsets.search import quick_search_values
from numba.typed import Dict
# from alphapept.feature_finding import grow_trail, get_trails, int_list_to_array, truncate
from alphapept.feature_finding import int_list_to_array, get_minpos, correlate, check_isotope_pattern_directed
DELTA_M = mass_dict['delta_M']
DELTA_S = mass_dict['delta_S']
maximum_offset = DELTA_M + DELTA_S
def find_feature_alphapept(mzml_file, debug = False):
    # print('i am in new')
    query_data = load_mzml_data(mzml_file)
    f_settings = {'centroid_tol': 25,
                  'hill_check_large': 40,
                  'hill_length_min': 3,
                  'hill_nboot': 150,
                  'hill_nboot_max': 100,
                  'hill_smoothing': 1,
                  # 'hill_smoothing':2,
                  'hill_split_level': 1.1,
                  'iso_charge_max': 1,
                  'iso_charge_min': 1,
                  'iso_corr_min': 0.6,
                  'iso_mass_range': 1,
                  'iso_n_seeds': 100,
                  'iso_split_level': 1.1,
                  'map_mob_range': 0.3,
                  # 'map_mz_range': 0.5,
                  'map_mz_range': 1,
                  'map_n_neighbors': 5,
                  'map_rt_range': 0.5,
                  'max_gap': 2,
                  'search_unidentified': True}
    # print(f_settings['iso_corr_min'])
    max_gap = f_settings['max_gap']
    centroid_tol = f_settings['centroid_tol']
    hill_split_level = f_settings['hill_split_level']
    iso_split_level = f_settings['iso_split_level']
    window = f_settings['hill_smoothing']
    hill_check_large = f_settings['hill_check_large']
    iso_charge_min = f_settings['iso_charge_min']
    iso_charge_max = f_settings['iso_charge_max']
    iso_n_seeds = f_settings['iso_n_seeds']
    hill_nboot_max = f_settings['hill_nboot_max']
    hill_nboot = f_settings['hill_nboot']
    iso_mass_range = f_settings['iso_mass_range']
    iso_corr_min = f_settings['iso_corr_min']
    query_data["int_list_ms1"] = np.array(query_data["int_list_ms1"]).astype(int)
    int_data = query_data['int_list_ms1']

    # logging.info('Feature finding on {}'.format(file_name))
    # logging.info(f'Hill extraction with centroid_tol {centroid_tol} and max_gap {max_gap}')

    hill_ptrs, hill_data, path_node_cnt, score_median, score_std = extract_hills(
        query_data, max_gap, centroid_tol)
    # if centroid_tol <= score_median + score_std * 3:
    # hill_ptrs, hill_data, path_node_cnt, score_median, score_std = extract_hills(
    #     query_data, max_gap, score_median + score_std * 3)
    print(score_median + score_std * 3)
    hill_ptrs, hill_data = remove_duplicate_hills(hill_ptrs, hill_data, path_node_cnt)
    hill_ptrs = split_hills(hill_ptrs, hill_data, int_data, hill_split_level=hill_split_level,
                            window=window)
    hill_data, hill_ptrs = filter_hills(
        hill_data, hill_ptrs, int_data, hill_check_large=hill_check_large, window=window)
    stats, sortindex_, idxs_upper, scan_idx, hill_data, hill_ptrs = get_hill_data(
        query_data, hill_ptrs, hill_data, hill_nboot_max=hill_nboot_max, hill_nboot=hill_nboot)
    pre_isotope_patterns = get_pre_isotope_patterns(
        stats, idxs_upper, sortindex_, hill_ptrs, hill_data, int_data, scan_idx, maximum_offset,
        iso_charge_min=iso_charge_min, iso_charge_max=iso_charge_max, iso_mass_range=iso_mass_range, cc_cutoff=iso_corr_min)
    stats_df = pd.DataFrame(stats, columns = ['mz','mz_std','int1', 'intensity_apex', 'rt_left', 'rt_end'])
    if debug == True:
        return(pre_isotope_patterns, stats_df, stats )

    isotope_patterns, iso_idx, isotope_charges = get_isotope_patterns(
        pre_isotope_patterns, hill_ptrs, hill_data, int_data, scan_idx, stats, sortindex_, averagine_aa, isotopes,
        iso_charge_min=iso_charge_min, iso_charge_max=iso_charge_max, iso_mass_range=iso_mass_range,
        iso_n_seeds=iso_n_seeds, cc_cutoff=iso_corr_min, iso_split_level=iso_split_level, callback=None)
    feature_table, lookup_idx = feature_finder_report(
        query_data, isotope_patterns, isotope_charges, iso_idx, stats, sortindex_, hill_ptrs, hill_data)
    return(feature_table, lookup_idx)


import numpy as np
import alphapept.performance
from numba.typed import List
from typing import Callable, Union
# def get_isotope_patterns(pre_isotope_patterns:list, hill_ptrs:np.ndarray, hill_data:np.ndarray, int_data:np.ndarray, scan_idx:np.ndarray, stats:np.ndarray, sortindex_:np.ndarray,  averagine_aa:Dict, isotopes:Dict, iso_charge_min:int = 1, iso_charge_max:int = 6, iso_mass_range:float = 5, iso_n_seeds:int = 100, cc_cutoff:float=0.6, iso_split_level:float = 1.3, callback:Union[Callable, None]=None) -> (np.ndarray, np.ndarray, np.ndarray):
#     """Wrapper function to iterate over pre_isotope_patterns.
#
#     Args:
#         pre_isotope_patterns (list): List of pre-isotope patterns.
#         hill_ptrs (np.ndarray): Array containing the bounds to the hill_data.
#         hill_data (np.ndarray): Array containing the indices to hills.
#         int_data (np.ndarray): Array containing the intensity to each centroid.
#         scan_idx (np.ndarray): Array containing the scan index for a centroid.
#         stats (np.ndarray): Stats array that contains summary statistics of hills.
#         sortindex_ (np.ndarray): Sortindex to access the hills from stats.
#         averagine_aa (Dict): Dict containing averagine masses.
#         isotopes (Dict): Dict containing isotopes.
#         iso_charge_min (int, optional): Minimum isotope charge. Defaults to 1.
#         iso_charge_max (int, optional): Maximum isotope charge. Defaults to 6.
#         iso_mass_range (float, optional): Mass search range. Defaults to 5.
#         iso_n_seeds (int, optional): Number of isotope seeds. Defaults to 100.
#         cc_cutoff (float, optional): Cuttoff for correlation.. Defaults to 0.6.
#         iso_split_level (float, optional): Isotope split level.. Defaults to 1.3.
#         callback (Union[Callable, None], optional): Callback function for progress. Defaults to None.
#     Returns:
#         list: List of isotope patterns.
#         np.ndarray: Iso idx.
#         np.ndarray: Array containing isotope charges.
#     """
#
#     isotope_patterns = []
#     isotope_charges = []
#
#     charge_range = List()
#
#     for i in range(iso_charge_min, iso_charge_max + 1):
#         charge_range.append(i)
#
#     isotope_patterns = []
#     isotope_charges = []
#
#     for idx, pre_pattern in enumerate(pre_isotope_patterns):
#         extract = True
#         while extract:
#             isotope_pattern, isotope_charge = isolate_isotope_pattern(np.array(pre_pattern), hill_ptrs, hill_data, int_data, scan_idx, stats, sortindex_, iso_mass_range, charge_range, averagine_aa, isotopes, iso_n_seeds, cc_cutoff, iso_split_level)
#             if isotope_pattern is None:
#                 length = 0
#             else:
#                 length = len(isotope_pattern)
#
#             if length > 1:
#                 isotope_charges.append(isotope_charge)
#                 isotope_patterns.append(isotope_pattern)
#
#                 pre_pattern = [_ for _ in pre_pattern if _ not in isotope_pattern]
#
#                 if len(pre_pattern) <= 1:
#                     extract = False
#             else:
#                 extract = False
#
#
#         if callback:
#             callback((idx+1)/len(pre_isotope_patterns))
#
#
#     iso_patterns = np.zeros(sum([len(_) for _ in isotope_patterns]), dtype=np.int64)
#
#     iso_idx = np.zeros(len(isotope_patterns)+1, dtype='int')
#
#
#     start = 0
#     for idx, _ in enumerate(isotope_patterns):
#         iso_patterns[start:start+len(_)] = _
#         start += len(_)
#         iso_idx[idx+1] = start
#
#
#
#     return iso_patterns, iso_idx, np.array(isotope_charges)
#
#
#
# @alphapept.performance.compile_function(compilation_mode="numba")
# def isolate_isotope_pattern(pre_pattern:np.ndarray, hill_ptrs:np.ndarray, hill_data:np.ndarray, int_data:np.ndarray, scan_idx:np.ndarray, stats:np.ndarray, sortindex_:np.ndarray, iso_mass_range:float, charge_range:List, averagine_aa:Dict, isotopes:Dict, iso_n_seeds:int, cc_cutoff:float, iso_split_level:float)->(np.ndarray, int):
#     """Isolate isotope patterns.
#
#     Args:
#         pre_pattern (np.ndarray): Pre isotope pattern.
#         hill_ptrs (np.ndarray): Array containing the bounds to the hill_data.
#         hill_data (np.ndarray): Array containing the indices to hills.
#         int_data (np.ndarray): Array containing the intensity to each centroid.
#         scan_idx (np.ndarray): Array containing the scan index for a centroid.
#         stats (np.ndarray): Stats array that contains summary statistics of hills.
#         sortindex_ (np.ndarray): Sortindex to access the hills from stats.
#         iso_mass_range (float): Mass range for checking isotope patterns.
#         charge_range (List): Charge range.
#         averagine_aa (Dict): Dict containing averagine masses.
#         isotopes (Dict): Dict containing isotopes.
#         iso_n_seeds (int): Number of seeds.
#         cc_cutoff (float): Cutoff value for what is considered correlating.
#         iso_split_level (float): Split level when isotopes are split.
#
#     Returns:
#         np.ndarray: Array with the best pattern.
#         int: Charge of the best pattern.
#     """
#     longest_trace = 0
#     champion_trace = None
#     champion_charge = 0
#     champion_intensity = 0
#
#     # Sort patterns by mass
#
#     sortindex = np.argsort(stats[pre_pattern][:,0]) #intensity
#     sorted_pattern = pre_pattern[sortindex]
#     massindex = np.argsort(stats[sorted_pattern][:,2])[::-1][:iso_n_seeds]
#
#     # Use all the elements in the pre_pattern as seed
#
#     for seed in massindex:  # Loop through all seeds
#         seed_global = sorted_pattern[seed]
#
#         trails = get_trails(seed, sorted_pattern, stats, charge_range, iso_mass_range, sortindex_, hill_ptrs, hill_data, int_data, scan_idx, cc_cutoff)
#
#         for index, trail in enumerate(trails):
#             if len(trail) >= longest_trace:  # Needs to be longer than the current champion
#                 arr = int_list_to_array(trail)
#                 intensity_profile = stats[arr][:,2]
#                 seedpos = np.nonzero(arr==seed_global)[0][0]
#
#                 # truncate around the seed...
#                 arr = truncate(arr, intensity_profile, seedpos, iso_split_level)
#                 intensity_profile = stats[arr][:,2]
#
#                 # Remove lower masses:
#                 # Take the index of the maximum and remove all masses on the left side
#                 if charge_range[index] * stats[seed_global, 0] < 1000:
#                     maxpos = np.argmax(intensity_profile)
#                     arr = arr[maxpos:]
#                     intensity_profile = stats[arr][:,2]
#
#                 if (len(arr) > longest_trace) | ((len(arr) == longest_trace) & (intensity_profile.sum() > champion_intensity)):
#                     # Averagine check
#                     cc = 1
#                     if cc > 0.6:
#                         # Update the champion
#                         champion_trace = arr
#                         champion_charge = charge_range[index]
#                         longest_trace = len(arr)
#                         champion_intensity = intensity_profile.sum()
#
#     return champion_trace, champion_charge
#
# @alphapept.performance.compile_function(compilation_mode="numba")
# def get_trails(seed:int, pattern:np.ndarray, stats:np.ndarray, charge_range:List, iso_mass_range:float, sortindex_:np.ndarray, hill_ptrs:np.ndarray, hill_data:np.ndarray, int_data:np.ndarray, scan_idx:np.ndarray, cc_cutoff:float)->List:
#     """Wrapper to extract trails for a given charge range.
#
#     Args:
#         seed (int): Seed index.
#         pattern (np.ndarray): Pre isotope pattern.
#         stats (np.ndarray): Stats array that contains summary statistics of hills.
#         charge_range (List): Charge range.
#         iso_mass_range (float): Mass range for checking isotope patterns.
#         sortindex_ (np.ndarray): Sortindex to access the hills from stats.
#         hill_ptrs (np.ndarray): Array containing the bounds to the hill_data.
#         hill_data (np.ndarray): Array containing the indices to hills.
#         int_data (np.ndarray): Array containing the intensity to each centroid.
#         scan_idx (np.ndarray): Array containing the scan index for a centroid.
#         cc_cutoff (float): Cutoff value for what is considered correlating.
#
#     Returns:
#         List: Trail of consistent hills.
#     """
#     trails = []
#     for charge in charge_range:
#         trail = grow_trail(seed, pattern, stats, charge, iso_mass_range, sortindex_, hill_ptrs, hill_data, int_data, scan_idx, cc_cutoff)
#
#         trails.append(trail)
#
#     return trails
#
# @alphapept.performance.compile_function(compilation_mode="numba")
# def grow_trail(seed:int, pattern:np.ndarray, stats:np.ndarray, charge:int, iso_mass_range:float, sortindex_:np.ndarray, hill_ptrs:np.ndarray, hill_data:np.ndarray, int_data:np.ndarray, scan_idx:np.ndarray, cc_cutoff:float)->List:
#     """Wrapper to grow an isotope pattern to the left and right side.
#
#     Args:
#         seed (int): Seed position.
#         pattern (np.ndarray): Isotope pattern.
#         stats (np.ndarray): Stats array that contains summary statistics of hills.
#         charge (int): Charge.
#         iso_mass_range (float): Mass range for checking isotope patterns.
#         sortindex_ (np.ndarray): Sortindex to access the hills from stats.
#         hill_ptrs (np.ndarray): Array containing the bounds to the hill_data.
#         hill_data (np.ndarray): Array containing the indices to hills.
#         int_data (np.ndarray): Array containing the intensity to each centroid.
#         scan_idx (np.ndarray): Array containing the scan index for a centroid.
#         cc_cutoff (float): Cutoff value for what is considered correlating.
#
#     Returns:
#         List: Isotope pattern.
#     """
#     x = pattern[seed]
#     trail = List()
#     trail.append(x)
#     trail = grow(trail, seed, -1, -1, 1, stats, pattern, charge, iso_mass_range, sortindex_, hill_ptrs, hill_data, int_data, scan_idx, cc_cutoff)
#     trail = grow(trail, seed, 1, 1, 1, stats, pattern, charge, iso_mass_range, sortindex_, hill_ptrs, hill_data, int_data, scan_idx, cc_cutoff)
#
#     return trail
#
# @alphapept.performance.compile_function(compilation_mode="numba")
# def grow(trail:List, seed:int, direction:int, relative_pos:int, index:int, stats:np.ndarray, pattern:np.ndarray, charge:int, iso_mass_range:float, sortindex_:np.ndarray, hill_ptrs:np.ndarray, hill_data:np.ndarray, int_data:np.ndarray, scan_idx:np.ndarray, cc_cutoff:float)->List:
#     """Grows isotope pattern based on a seed and direction.
#
#     Args:
#         trail (List): List of hills belonging to a pattern.
#         seed (int): Seed position.
#         direction (int): Direction in which to grow the trail
#         relative_pos (int): Relative position.
#         index (int): Index.
#         stats (np.ndarray): Stats array that contains summary statistics of hills.
#         pattern (np.ndarray): Isotope pattern.
#         charge (int): Charge.
#         iso_mass_range (float): Mass range for checking isotope patterns.
#         sortindex_ (np.ndarray): Sortindex to access the hills from stats.
#         hill_ptrs (np.ndarray): Array containing the bounds to the hill_data.
#         hill_data (np.ndarray): Array containing the indices to hills.
#         int_data (np.ndarray): Array containing the intensity to each centroid.
#         scan_idx (np.ndarray): Array containing the scan index for a centroid.
#         cc_cutoff (float): Cutoff value for what is considered correlating.
#
#     Returns:
#         List: List of hills belonging to a pattern.
#     """
#     x = pattern[seed]  # This is the seed
#     mass1 = stats[x,0]
#     delta_mass1 = stats[x,1]
#
#     k = sortindex_[x]
#     start = hill_ptrs[k]
#     end = hill_ptrs[k + 1]
#     idx_ = hill_data[start:end]
#     int_ = int_data[idx_]
#     scans_ = scan_idx[idx_]
#
#     growing = True
#
#     while growing:
#         if direction == 1:
#             if seed + relative_pos == len(pattern):
#                 growing = False
#                 break
#         else:
#             if seed + relative_pos < 0:
#                 growing = False
#                 break
#
#         y = pattern[seed + relative_pos]  # This is a reference peak
#
#         l = sortindex_[y]
#
#         mass2 = stats[y,0]
#         delta_mass2 = stats[y,1]
#
#         start = hill_ptrs[l]
#         end = hill_ptrs[l + 1]
#         idx_ = hill_data[start:end]
#         int_2 = int_data[idx_]
#         scans_2 = scan_idx[idx_]
#
#         if correlate(scans_, scans_2, int_, int_2) > cc_cutoff:
#             if check_isotope_pattern_directed(mass1, mass2, delta_mass1, delta_mass2, charge, -direction * index, iso_mass_range):
#                 if direction == 1:
#                     trail.append(y)
#                 else:
#                     trail.insert(0, y)
#                 index += (
#                     1
#                 )  # Greedy matching: Only one edge for a specific distance, will not affect the following matches
#
#         delta_mass = np.abs(mass1 - mass2)
#
#         if (delta_mass > (maximum_offset) * index):  # the pattern is sorted so there is a maximum to look back
#             break
#
#         relative_pos += direction
#
#     return trail
#
#
# @alphapept.performance.compile_function(compilation_mode="numba")
# def truncate(array:np.ndarray, intensity_profile:np.ndarray, seedpos:int, iso_split_level:float)->np.ndarray:
#     """Function to truncate an intensity profile around its seedposition.
#
#     Args:
#         array (np.ndarray):  Input array.
#         intensity_profile (np.ndarray): Intensities for the input array.
#         seedpos (int): Seedposition.
#         iso_split_level (float): Split level.
#
#     Returns:
#         np.ndarray: Truncated array.
#     """
#     minima = int_list_to_array(get_minpos(intensity_profile, iso_split_level))
#
#     if len(minima) > 0:
#         left_minima = minima[minima < seedpos]
#         right_minima = minima[minima > seedpos]
#
#         # If the minimum is smaller than the seed
#         if len(left_minima) > 0:
#             minpos = left_minima[-1]
#         else:
#             minpos = 0
#
#         if len(right_minima) > 0:
#             maxpos = right_minima[0]
#         else:
#             maxpos = len(array)
#
#         array = array[minpos:maxpos+1]
#
#     return array