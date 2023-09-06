import os
import shutil
import sys

all_folder = sys.argv[1]
# pos_neg_folder = '/Volumes/Brother_cow/brighseed/2068_AX_M6P1/2068_AX_M6P1'
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