import numpy as np
def build_index(mass_flatten, intensity_flatten, index_flatten):

    mass_sorted, intensity_sorted, index_sorted = zip(*sorted(zip(mass_flatten, intensity_flatten, index_flatten)))
    mass_sorted = np.array(mass_sorted)
    intensity_sorted= np.array(intensity_sorted)
    index_sorted = np.array(index_sorted)
    return(mass_sorted, intensity_sorted, index_sorted)

def flash_eic(pmz, mass_error, mass_sorted, intensity_sorted, index_sorted):
    index_start, index_end = mass_sorted.searchsorted([pmz-mass_error, pmz+mass_error])
    index_range = index_sorted[index_start:index_end]
    intensity_range = intensity_sorted[index_start:index_end]
    intensity_list = np.zeros(np.max(index_sorted)+1)
    for idx in range(0,len(index_range)):
        intensity_list[index_range[idx]]= intensity_list[index_range[idx]]+intensity_range[idx]
    return(intensity_list)