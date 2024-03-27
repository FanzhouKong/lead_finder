import time

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from toolsets.file_io import get_file_list
import toolsets.T_rex as trx
import toolsets.raw_data_scaffold as rds
# mzml_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset/mzml'
# pl_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset/pl_new'
# eic_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset/eic'
master_parent = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/Alkaloids_lib'
master_child_qe = ['QE_neg', 'QE_pos', ]
master_child_ttof = ['TTOF_neg', 'TTOF_pos']
mzml_dir = 'mzml'
pl_dir = 'pl'
# from toolsets.search import search_feature, quick_search_values
# for f in [mzml_dir, pl_dir]:
#     for c in master_child:
#         if os.path.exists(os.path.join(master_parent, c, f))==False:
#             os.makedirs(os.path.join(master_parent, c, f))

from tqdm import tqdm
# for child in master_child_qe:
#     bad_files = []
#     mzml = os.path.join(master_parent, child, mzml_dir)
#     output_dir = os.path.join(master_parent, child, pl_dir)
#     file_list = get_file_list(mzml, '.mzML')
#     for f in tqdm(file_list):
#         ms1, ms2 = rds.read_mzml(f, mzml)
#         mass_sorted, intensity_sorted, index_sorted, rt_list = trx.build_index(ms1)
#         try:
#             df = trx.get_features(mass_sorted, intensity_sorted, index_sorted, rt_list,
#                              intensity_threshold=30000, n_neighbor = 2)
#             df.to_csv(os.path.join(output_dir, f+'.csv'), index = False)
#         except:
#             bad_files.append(f)
#     if len(bad_files)>0:
#         bad_files = pd.DataFrame(bad_files)
#         bad_files.to_csv(os.path.join(output_dir,'bad_files.csv'), index = False)

        # df.to_csv(os.path.join(pl_dir, f+'.csv'), index = False)
print('start processing qttof data')
for child in master_child_ttof:
    bad_files = []
    mzml = os.path.join(master_parent, child, mzml_dir)
    output_dir = os.path.join(master_parent, child, pl_dir)
    file_list = get_file_list(mzml, '.mzML')
    for f in tqdm(file_list):
        try:
            df = trx.find_feature(f, mzml, intensity_threshold=1000, n_neighbor=1)

            df[0].to_csv(os.path.join(output_dir, f+'.csv'), index = False)
        except:
            bad_files.append(f)
    if len(bad_files)>0:
        bad_files = pd.DataFrame(bad_files)
        bad_files.to_csv(os.path.join(output_dir,'bad_files.csv'), index = False)