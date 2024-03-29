import numpy as np
import scipy.stats
import scipy.special
from typing import Union
from .tools import convert_spectrum_to_numpy_array, clean_spectrum, match_two_spectra
from .tools_spectral_entropy import unweighted_entropy_similarity as fast_entropy_similarity, \
    spectral_entropy as spectral_entropy_fast, intensity_entropy, apply_weight_to_intensity


def entropy_similarity(spectrum_1: Union[list, np.ndarray], spectrum_2: Union[list, np.ndarray],
                       ms2_ppm: float = None, ms2_da: float = 0.05,
                       clean_spectra: bool = True) -> float:
    """
    Calculate the entropy similarity between two spectra.

    :param spectrum_1: The first spectrum.
    :param spectrum_2: The second spectrum.
    :param ms2_ppm: The mass accuracy in ppm.
    :param ms2_da: The mass accuracy in Da.
    :param clean_spectra: Whether to clean the spectra before calculating the similarity or not.
                          It can only be set to False when the spectra are already cleaned by the "clean_spectrum" function.
    :return: The entropy similarity.
    """

    # Clean the spectra if needed.
    if clean_spectra:
        spectrum_1 = clean_spectrum(spectrum_1, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        spectrum_2 = clean_spectrum(spectrum_2, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    if spectrum_1.shape[0] == 0 or spectrum_2.shape[0] == 0:
        return 0.

    # Weight the spectra.
    spectrum_1 = apply_weight_to_intensity(spectrum_1)
    spectrum_2 = apply_weight_to_intensity(spectrum_2)

    # Calculate the similarity.
    # entropy_1 = np.sum(scipy.special.entr(spectrum_1[:, 1]))
    # entropy_2 = np.sum(scipy.special.entr(spectrum_2[:, 1]))
    # entropy_merged = merged_spectral_entropy_log2(spectrum_1, spectrum_2, ms2_ppm=ms2_ppm, ms2_da=ms2_da) * np.log(2)
    # return 1 - (2 * entropy_merged - entropy_1 - entropy_2) / np.log(4)
    similarity = fast_entropy_similarity(spectrum_1, spectrum_2, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    return similarity
    # entropy_1 = spectral_entropy_log2(spectrum_1)
    # entropy_2 = spectral_entropy_log2(spectrum_2)
    # return 1 - entropy_merged + 0.5 * (entropy_1 + entropy_2)


def reverse_entropy_similarity(spectrum_query: Union[list, np.array], spectrum_library: Union[list, np.array],
                               ms2_ppm: float = None, ms2_da: float = 0.05,
                               clean_spectra: bool = True) -> float:
    """
    Calculate the reverse similarity between a query spectrum and a library spectrum.
    :param spectrum_query: The query spectrum.
    :param spectrum_library: The library spectrum.
    :param ms2_ppm: The mass accuracy in ppm.
    :param ms2_da: The mass accuracy in Da.
    :param clean_spectra: Whether to clean the spectra before calculating the similarity or not.
    :return: The reverse similarity and the total uesd peaks' intensity in the query spectrum.
    """
    if clean_spectra:
        spectrum_query = clean_spectrum(spectrum_query, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        spectrum_library = clean_spectrum(spectrum_library, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    if spectrum_1.shape[0] == 0 or spectrum_2.shape[0] == 0:
        return 0.

    # Select peaks only exist in the library spectrum.
    match_result = match_two_spectra(spectrum_query, spectrum_library, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    spectrum_query = match_result[np.bitwise_and(match_result[:, 1] > 0, match_result[:, 2] > 0), 0:2]
    spectrum_query_sum = np.sum(spectrum_query[:, 1])
    if spectrum_query_sum > 0:
        spectrum_query[:, 1] = spectrum_query[:, 1] / spectrum_query_sum
        return entropy_similarity(spectrum_query, spectrum_library, ms2_ppm=ms2_ppm, ms2_da=ms2_da, clean_spectra=False), spectrum_query_sum
    else:
        return 0, 0


def unweighted_entropy_similarity(spectrum_1: np.ndarray, spectrum_2: np.ndarray,
                                  ms2_ppm: float = None, ms2_da: float = 0.05,
                                  clean_spectra: bool = True) -> float:
    """
    Calculate the unweighted entropy similarity between two spectra.
    :param spectrum_1: The first spectrum.
    :param spectrum_2: The second spectrum.
    :param ms2_ppm: The mass accuracy in ppm.
    :param ms2_da: The mass accuracy in Da.
    :param clean_spectra: Whether to clean the spectra before calculating the similarity or not.
    :return: The unweighted entropy similarity.
    """
    if clean_spectra:
        spectrum_1 = clean_spectrum(spectrum_1, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        spectrum_2 = clean_spectrum(spectrum_2, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    if spectrum_1.shape[0] == 0 or spectrum_2.shape[0] == 0:
        return 0.
    
    # Calculate the similarity.
    return fast_entropy_similarity(spectrum_1, spectrum_2, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    # entropy_1 = np.sum(scipy.special.entr(spectrum_1[:, 1]))
    # entropy_2 = np.sum(scipy.special.entr(spectrum_2[:, 1]))
    # entropy_merged = merged_spectral_entropy(spectrum_1, spectrum_2, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    # return 1 - (2 * entropy_merged - entropy_1 - entropy_2) / np.log(4)
    # entropy_1 = spectral_entropy_log2(spectrum_1)
    # entropy_2 = spectral_entropy_log2(spectrum_2)
    # entropy_merged = merged_spectral_entropy_log2(spectrum_1, spectrum_2, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    # return 1 - entropy_merged + 0.5 * (entropy_1 + entropy_2)


def spectral_entropy(spectrum: np.ndarray,
                     max_mz: float = None,
                     noise_threshold: float = 0.01,
                     ms2_ppm: float = None, ms2_da: float = 0.05,
                     clean_spectra: bool = True) -> float:
    """
    Calculate the spectral entropy of the spectrum.
    :param spectrum: The spectrum.
    :param max_mz: The maximum m/z value.
    :param noise_threshold: The noise threshold.
    :param ms2_ppm: The mass accuracy in ppm.
    :param ms2_da: The mass accuracy in Da.
    :param clean_spectra: Whether to clean the spectra before calculating the spectral entropy or not.
    :return: The spectral entropy.
    """

    if clean_spectra:
        spectrum = clean_spectrum(spectrum,
                                  max_mz=max_mz, noise_threshold=noise_threshold,
                                  ms2_ppm=ms2_ppm, ms2_da=ms2_da)
    else:
        spectrum = convert_spectrum_to_numpy_array(spectrum)

    if spectrum.shape[0] == 0:
        return 0.
    else:
        return spectral_entropy_fast(spectrum)


def weight_spectrum(intensity: np.ndarray) -> np.ndarray:
    """
    Weight the spectral intensity. To save time, the intensity should be normalized to sum to 1.
    :param intensity: The spectral intensity.
    :return: The weighted spectral intensity.
    """
    # Calculate spectral entropy.
    spectrum_entropy = intensity_entropy(intensity)
    if spectrum_entropy >= 3:
        return intensity
    else:
        # Weight the spectrum.
        weight = 0.25 + 0.25 * spectrum_entropy
        intensity_weighted = np.power(intensity, weight)
        intensity_sum = np.sum(intensity_weighted)

        if intensity_sum == 0:
            return intensity_weighted
        else:
            return intensity_weighted / intensity_sum
