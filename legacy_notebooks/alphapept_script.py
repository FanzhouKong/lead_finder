from feature_finding.feature_finding_alphapept import find_feature, find_feature_alphapept
import os
from tqdm import tqdm
import pandas as pd
mzml_dir = '/Users/fanzhou.kong/Dropbox (Brightseed)/Mac/Documents/GitHub/data_garage/brighseed/MZMLS/2068_AX_M6P1/pos'
feature_folder = '/Users/fanzhou.kong/Dropbox (Brightseed)/Mac/Documents/GitHub/data_garage/brighseed/alphapept_features/2068_AX_M6P1/pos'
if os.path.exists(feature_folder)== False:
    os.makedirs(feature_folder)
file_name_lists = []
for root, dirs, files in os.walk(mzml_dir):
    for file in tqdm(files, total = len(files)):
        if file.endswith('.mzML'):
            file_name_lists.append(file)
            # feature_table, lookup = find_feature_alphapept(os.path.join(mzml_dir,file))
for file_name in tqdm(file_name_lists):
    feature_table_temp, lookup = find_feature_alphapept(os.path.join(mzml_dir,file_name))
    base_name = file_name.split('.')[0]
    feature_table_temp.to_csv(os.path.join(feature_folder, base_name+'.csv'), index = False)