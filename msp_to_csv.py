import os.path
from toolsets.file_io import read_msp_files
import pandas as pd
msp_dir = '/Volumes/scratch/Fanzhou/library'
msp_name = 'neg_orbi.msp'
csv_dir = '/Users/fanzhoukong/Documents/GitHub/Libgen_data/curated_library/csv'
csv_name = msp_name.split('.')[0]+'.csv'
msp = read_msp_files(os.path.join(msp_dir, msp_name))
print('i am done read in msp files')
msp.to_csv(os.path.join(csv_dir, csv_name))
print('i am done output into csv files')