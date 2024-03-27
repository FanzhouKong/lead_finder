import pandas as pd
from molmass import Formula
from functools import reduce
import os
import numpy as np
from fuzzywuzzy import fuzz

import re

def find_floats(input_string):
    # Regular expression pattern for matching floats
    # This pattern matches an optional '+' or '-' sign, followed by one or more digits,
    # an optional decimal point, and zero or more digits after the decimal point
    float_pattern = r'[-+]?\d*\.\d+'

    # Find all matches of the pattern in the input string
    # The findall function returns all non-overlapping matches of the pattern in the string, as a list of strings
    floats = re.findall(float_pattern, input_string)

    # Convert found string matches to float values
    float_values = [float(f) for f in floats]

    return float_values



def save_value_counts(data, column):
    vc = data[column].value_counts().rename_axis('unique_values').to_frame('counts')
    vc.index.name = 'unique_values'
    vc.reset_index(inplace=True)
    return(vc)
def check_missing_compound(lst1, lst2):
    set1 = set(lst1)
    set2 = set(lst2)

    # Find items that are in either set1 or set2 but not in both using symmetric difference
    unique_items = set1.symmetric_difference(set2)

    # Convert the set back to a list if needed
    unique_items_list = list(unique_items)
    return(unique_items_list)
def find_files(base, pattern):
    # print('i am in new new')
    score = []
    for filename in os.listdir(base):
        score.append(fuzz.ratio(filename,pattern))
    return(os.listdir(base)[np.argmax(score)])
def specify_column(keyword, column_name):
    score = []
    for name in column_name:
        score.append(fuzz.token_sort_ratio(keyword,name))
    return(column_name[np.argmax(score)])




