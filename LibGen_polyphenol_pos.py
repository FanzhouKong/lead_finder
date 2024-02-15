
import numpy as np
import pandas as pd
from toolsets.search import string_search, quick_search_values, quick_search_sorted
from tqdm import tqdm
import os
from toolsets.file_io import get_file_list
from rdkit import Chem
from toolsets.ff_droup import feature_finding, auto_EIC, EIC, get_feature
from toolsets.raw_data_scaffold import read_mzml
from toolsets.std_list_prep import neutrilize_salt, neutrilize_salt_df, std_list_cleanup
from toolsets.std_list_prep import complete_adducts
from toolsets.spectra_operations import entropy_identity
from toolsets.feature_std_matching import feature_matching
master_dir ='/Volumes/Samsung_T5/polyphenol_lib/'
mzml_dir = '/Volumes/Samsung_T5/polyphenol_lib/neg_high_ce/mzml'
adducts = ['[M]-','[M-H]-',  '[M+Cl]-', '[M-H2O-H]-']
pl_dir = '/Volumes/Samsung_T5/polyphenol_lib/neg_high_ce/pl'
if os.path.exists(pl_dir)==False:
    # print('i am in if')
    os.makedirs(pl_dir)
std_list_pos_adduct=pd.read_csv(os.path.join(master_dir, 'std_list_neg_adduct.csv'))
matched_all = pd.DataFrame()
p_files = []
for mix in tqdm(std_list_pos_adduct['mix'].unique()):
    # ms1, ms2 = read_mzml(mix, mzml_dir)
    # feature_mix = get_feature(ms1, ms2, filter=True)
    # std_list_mix = string_search(std_list_pos_adduct, 'mix', mix)
    # matched_mix = feature_matching(feature_mix, std_list_mix, adducts = adducts)

    try:
        ms1, ms2 = read_mzml(mix, mzml_dir)
        feature_mix = get_feature(ms1, ms2, filter=True)
        if len(feature_mix)>0:
            feature_mix.to_csv(os.path.join(pl_dir, mix+'_features.csv'), index = False)
        std_list_mix = string_search(std_list_pos_adduct, 'mix', mix)
        matched_mix = feature_matching(feature_mix, std_list_mix, adducts = adducts)

        matched_mix.to_csv(os.path.join(pl_dir, mix+'_matched.csv'), index = False)
    except:
        print(mix)


    # matched_all = pd.concat([matched_all, matched_mix], ignore_index=True)
    # try:
    #     matched_mix = feature_matching(feature_mix, std_list_mix, adducts = adducts)
    #     matched_all = pd.concat([matched_all, matched_mix], ignore_index=True)
    # except:
    #     print(mix)
    #     p_files.append(mix)
# matched_all.to_csv(os.path.join(master_dir, 'pos_matched.csv'), index=False)