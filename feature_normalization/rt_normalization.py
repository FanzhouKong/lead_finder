from feature_finding.feature_finding_alphapept import find_feature
def find_seeded(seeded, test_file):
    test_rt = []
    test_intensity = []
    # print('i am in new')
    for index, row in seeded.iterrows():
        istd_temp = find_feature(test_file, mz = row['pmz_protonated'], rt = row['RT'], mz_column= 'Precursor m/z', rt_column='RT_adjusted')
        if len(istd_temp)>0:
            istd_temp.sort_values(by = 'Height', ascending=False, inplace=True)
            test_rt.append(istd_temp.iloc[0]['RT_adjusted'])
            test_intensity.append(istd_temp.iloc[0]['Height'])
    return(test_rt, test_intensity)
    # pass



def find_istd(istd_info, sample_file, rt_column = 'RT (min)', return_adjusted_rt = False):
    seed_rt = []
    seed_intensity = []
    # print('i am in new')
    for index, row in istd_info.iterrows():
        istd_temp = find_feature(sample_file, mz = row['[M+H]+'], rt = row['RT_suggested'], mz_column= 'Precursor m/z', rt_column=rt_column)
        if len(istd_temp)>0:
            istd_temp.sort_values(by = 'Height', ascending=False, inplace=True)
            if return_adjusted_rt == False:
                seed_rt.append(istd_temp.iloc[0]['RT (min)'])
            else:
                seed_rt.append(istd_temp.iloc[0]['RT_adjusted'])
            seed_intensity.append(istd_temp.iloc[0]['Height'])
    return(seed_rt, seed_intensity)