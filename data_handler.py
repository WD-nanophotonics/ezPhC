import numpy as np

'This file stores relatively irrelavent functions'

def get_data_structure(data):
    if isinstance(data, np.ndarray):
        return f'np.array:{data.shape}'
    elif isinstance(data, list):
        sub_structure = []
        for item in data:
            sub_structure.append(get_data_structure(item))

        return sub_structure

    return 'unknown'

def is_one_layer_pattern_list(data): # np arrays stored in a list, the standard data structure for lattice dots.
    if isinstance(data, list):
        return all(isinstance(item, np.ndarray) for item in data)
    return False

def is_two_layer_pattern_list(data):
    if isinstance(data, list):
        return all(is_one_layer_pattern_list(item) for item in data)
    return False


def convet_to_two_layer_pattern_list(data):
    result=[]
    if isinstance(data, list):
        if is_two_layer_pattern_list(data):
            return data
        elif is_one_layer_pattern_list(data):
            return [data]
        else:
            for item in data:
                result+=convet_to_two_layer_pattern_list(item)
    elif isinstance(data, np.ndarray):
        return [[data]]
    else:
        raise ValueError("Invalid data type provided.")