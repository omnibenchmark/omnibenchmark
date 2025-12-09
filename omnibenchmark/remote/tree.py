from collections import defaultdict


def tree_string_from_list(file_list):
    file_tree = _tree()

    for file_path in file_list:
        parts = file_path.parts
        current_level = file_tree
        for part in parts:
            current_level = current_level[part]

    def print_tree(current_level, indent=""):
        out = ""
        for key, subtree in current_level.items():
            out += f"{indent}{key}/\n"
            out += print_tree(subtree, indent + "    ")
        return out

    return print_tree(file_tree)


def _tree():
    """
    A recursive function that creates a tree structure from a list of file paths.
    """
    return defaultdict(_tree)
