# from: https://stackoverflow.com/a/1094933
# TODO(ben): use humanize
def sizeof_fmt(num: float, suffix: str = "B"):
    if abs(num) < 1024:
        return f"{num: 5d}{suffix}"
    for unit in ("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"):
        if abs(num) < 1024:
            return f"{num:3.1f}{unit}{suffix}"
        num /= 1024
    return f"{num:.1f}Yi{suffix}"
