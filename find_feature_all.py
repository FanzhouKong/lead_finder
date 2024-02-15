import time

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from toolsets.file_io import get_file_list
import toolsets.T_rex as trx
import toolsets.raw_data_scaffold as rds
mzml_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset/mzml'
pl_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset/pl_new'
eic_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset/eic'
master_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/benchmarking_dataset'
from toolsets.search import search_feature, quick_search_values
for f in [mzml_dir, pl_dir, eic_dir]:
    if os.path.exists(f)==False:
        os.makedirs(f)
file_list = get_file_list(mzml_dir, '.mzML', with_tail=False)
from tqdm import tqdm
for f in tqdm(file_list[1:]):
    ms1, ms2 = rds.read_mzml(f, mzml_dir)
    mass_sorted, intensity_sorted, index_sorted, rt_list = trx.build_index(ms1)
    try:
        df = trx.get_features(mass_sorted, intensity_sorted, index_sorted, rt_list,
                              base_name=ms1.iloc[0]['base_name'], intensity_threshold=30000, n_neighbor = 10)
        df.to_csv(os.path.join(pl_dir, f+'.csv'), index = False)
    except:
        print(f)