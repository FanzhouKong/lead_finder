import subprocess
import os

import shutil
def get_dirs(workplace, mode = None):
    first_dirs = ['alignment_result', 'bioactivity', 'mzml', 'normalized_peak_list',
                  'peak_list', 'results', 'sirius_files']
    return_dirs = []
    if mode is not None:
        for item in first_dirs:
            if item != 'bioactivity':
                if os.path.exists(os.path.join(workplace, item, mode))== True:
                    return_dirs.append(os.path.join(workplace, item, mode))
            else:
                if os.path.exists(os.path.join(workplace, item))== True:
                    return_dirs.append(os.path.join(workplace, item))
    return return_dirs
def set_workspace(workplace, mode = None):
    if os.path.exists(workplace)==False:
        os.makedirs(workplace)
    else:
        print("the workspace exists!")
    # print('tt')
    first_dirs = ['alignment_result', 'bioactivity', 'mzml', 'normalized_peak_list',
                  'peak_list', 'results', 'sirius_files']
    second_dir = ['pos', 'neg']
    for item in first_dirs:
        if item != 'bioactivity':
            if os.path.exists(os.path.join(workplace, item))== False:
                os.makedirs(os.path.join(workplace, item))
                for item_sec in second_dir:
                    os.makedirs(os.path.join(workplace, item, item_sec))
                    if item_sec==mode:
                        # print('i am in if')
                        dir_to_return.append(os.path.join(workplace, item, item_sec))
        else:
            if os.path.exists(os.path.join(workplace, item))== False:
                os.makedirs(os.path.join(workplace, item))
                dir_to_return.append(os.path.join(workplace, item))
    print('set up complete')
def split_pos_neg(all_folder):
    pos_folder = os.path.join(all_folder, 'pos')
    neg_folder = os.path.join(all_folder, 'neg')
    for folder in [pos_folder, neg_folder]:
        if os.path.exists(folder)==False:
            os.makedirs(folder)

    # file_lists = []
    for root, dirs, files in os.walk(all_folder):
        for file in files:
            # print(file)

            if file.endswith('.mzML'):
                if len(file.split('.'))==2:
                    base_name = file.split('.')[0]
                    # print(base_name)
                    if base_name[-1]=='P':
                        shutil.move(os.path.join(all_folder, file), os.path.join(pos_folder, file))
                    elif base_name[-1]=='N':
                        shutil.move(os.path.join(all_folder, file), os.path.join(neg_folder, file))