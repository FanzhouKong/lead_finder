import os

from toolsets.file_io import get_file_list
from toolsets.ff_droup import process_mzml, feature_finding
from tqdm import tqdm
import pandas as pd
import sys
# start = int(sys.argv[1])
# end = int(sys.argv[2])
mzml_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/gut_microbiom/HILIC pos mode mzml'
pl_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/gut_microbiom/HILIC pos mode mzml_pl'
if __name__ == '__main__':
    if os.path.exists(pl_dir)==False:
        os.mkdir(pl_dir)
    file_list = get_file_list(mzml_dir, '.mzML', with_tail=False)
    bad_files = []
    # file_list = file_list[100:]
    for file in tqdm(file_list):
        try:
            ms1, ms2 = process_mzml(file+'.mzML', parent_dir=mzml_dir, rt_max=5)
            features = feature_finding(ms1, ms2)
            features = features[features['snr']>3]
            features.to_csv(os.path.join(pl_dir, file+'.csv'))
        except:
            bad_files.append(file)
            # break
    if len(bad_files)>0:
        with open(os.path.join(pl_dir, 'bad_files.txt'), 'w') as f:
            for file in bad_files:
                f.write(file)
                f.write('\n')
        f.close()