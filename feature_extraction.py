from tqdm import tqdm
import numpy as np
import os
from toolsets.file_io import get_file_list
import pandas as pd

import toolsets.raw_data_scaffold as rds
from toolsets.file_io import readin_peak_list
import toolsets.T_rex as trx
import pandas as pd
import sys
# start = int(sys.argv[1])
# end = int(sys.argv[2])
mzml_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/exposome_negBA'
pl_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/exposome_negBA_pl'
if __name__ == '__main__':
    if os.path.exists(pl_dir)==False:
        os.mkdir(pl_dir)
    file_list = get_file_list(mzml_dir, '.mzML', with_tail=False)
    bad_files = []
    # file_list = file_list[100:]
    for file in tqdm(file_list):
        try:
            featrue_temp = trx.find_features(file, mzml_dir)
            featrue_temp.to_csv(os.path.join(pl_dir, file+'.csv'), index = False)
        except:
            bad_files.append(file)
            # break
    if len(bad_files)>0:
        with open(os.path.join(pl_dir, 'bad_files.txt'), 'w') as f:
            for file in bad_files:
                f.write(file)
                f.write('\n')
        f.close()