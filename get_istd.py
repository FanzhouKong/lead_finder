import pandas as pd
import os
mzml_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/gut_microbiom/HILIC pos mode mzml'
istd_info = pd.read_csv(os.path.join('/Users/fanzhoukong/Documents/GitHub/Libgen_data/gut_microbiom', 'CORE_posHILIC_mzrt_5minMtd_correctionSheet.csv'))
istd_info_rt = istd_info.copy()
istd_info_intensity = istd_info.copy()
from toolsets.ff_droup import get_istd_info_all
from toolsets.file_io import get_file_list
import time
if __name__ == '__main__':
    file_list = get_file_list(mzml_dir, '.mzML', with_tail=False)
    print('got file list')
    start = time.time()
    results = get_istd_info_all(istd_info, file_list, mzml_dir)
    end = time.time()
    for r in results:
        istd_info_rt[r[0]]=r[1]
        istd_info_intensity[r[0]]=r[2]
    istd_info_rt.to_csv(os.path.join('/Users/fanzhoukong/Documents/GitHub/Libgen_data/gut_microbiom', 'istd_summary_rt.csv'), index = False)
    istd_info_intensity.to_csv(os.path.join('/Users/fanzhoukong/Documents/GitHub/Libgen_data/gut_microbiom', 'istd_summary_intensity.csv'), index = False)
    print('i am done extrcting istd info, time used is: ')
    print(str((end-start)/60))

