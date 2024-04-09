""" General utils functions"""
from typing import List, Union, Any


def as_list(input: Union[List, Any]):
    return input if isinstance(input, List) else [input]
