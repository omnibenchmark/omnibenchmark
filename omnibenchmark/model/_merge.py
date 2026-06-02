"""Dictionary merge helper for assembling stage/module outputs."""


def merge_dict_list(list_of_dicts):
    """Merge a list of dictionaries into a single dictionary."""
    merged_dict = {
        key: value for d in list_of_dicts if d is not None for key, value in d.items()
    }

    return merged_dict
