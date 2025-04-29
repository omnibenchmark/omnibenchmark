import os
from pathlib import Path
from typing import Dict, List, Tuple, Union
import glob

from omnibenchmark.benchmark import Benchmark
from omnibenchmark.io.exception import MinIOStorageVersioningCorruptionException
from omnibenchmark.io.S3config import S3_DEFAULT_STORAGE_OPTIONS


def get_objects_to_tag(
    objdic: Dict,
    storage_options: S3_DEFAULT_STORAGE_OPTIONS = S3_DEFAULT_STORAGE_OPTIONS(),
) -> Tuple[List, List]:
    """
    Get a list of objects that need to be tagged with the current version.

    Args:
    - objdic: A dictionary of objects and associated metadata.
    - storage_options: The storage options.

    Returns:
    - A list of object names and a list of version ids.
    """
    object_names = sorted(list(objdic.keys()))
    object_names_to_tag = []
    versionid_of_objects_to_tag = []
    for object_name in object_names:
        # get newest version
        newest_version = sorted(
            objdic[object_name].items(),
            key=lambda it: it[1]["last_modified"],
            reverse=True,
        )[0][0]
        # get root directory
        object_name_red = object_name
        while not os.path.split(object_name_red)[0] == "":
            object_name_red = os.path.split(object_name_red)[0]
            if object_name_red in storage_options.tracked_directories:
                break
        # check if newest version exists and
        # if object is in tracked directories
        if (object_name_red in storage_options.tracked_directories) and (
            not objdic[object_name][newest_version]["is_delete_marker"]
        ):
            object_names_to_tag.append(object_name)
            versionid_of_objects_to_tag.append(newest_version)
    return object_names_to_tag, versionid_of_objects_to_tag


def get_expected_benchmark_output_files(
    benchmark: Benchmark,
    storage_options: S3_DEFAULT_STORAGE_OPTIONS = S3_DEFAULT_STORAGE_OPTIONS(),
) -> List:
    object_names_to_keep = benchmark.get_output_paths()
    if storage_options.extra_files_to_version_not_in_benchmark_yaml:
        for (
            glob_expression
        ) in storage_options.extra_files_to_version_not_in_benchmark_yaml:
            found_files = glob.glob(glob_expression, recursive=True)
            if isinstance(found_files, str):
                found_files = [found_files]
            for found_file in found_files:
                object_names_to_keep.add(found_file)
    return list(object_names_to_keep)


# filter objects based on workflow.
def filter_objects_to_tag(
    object_names_to_tag: List,
    versionid_of_objects_to_tag: List,
    benchmark: Union[None, Benchmark] = None,
    storage_options: S3_DEFAULT_STORAGE_OPTIONS = S3_DEFAULT_STORAGE_OPTIONS(),
) -> Tuple[List, List]:
    if benchmark is None:
        return object_names_to_tag, versionid_of_objects_to_tag
    else:
        object_names_to_keep = get_expected_benchmark_output_files(
            benchmark, storage_options
        )

        rootdirs = [Path(obj).parents[-2].name for obj in object_names_to_tag]
        overwrite_dirs_to_keep = [
            t
            for t in storage_options.tracked_directories
            if t not in storage_options.results_directories
        ]
        return zip(
            *[
                (obj, versionid_of_objects_to_tag[i])
                for i, obj in enumerate(object_names_to_tag)
                if obj in object_names_to_keep or rootdirs[i] in overwrite_dirs_to_keep
            ]
        )


def get_single_remoteversion_from_bmversion(
    di: Dict, object_name: str, query_version: str
) -> Union[str, None]:
    """
    Get the object version where the tag matches the benchmark version.

    Args:
    - di: A dictionary of objects and associated metadata.
    - object_name: The name of the object.
    - query_version: The benchmark version to query.

    Returns:
    - The object version or None where the tag matches the benchmark version.
    """
    version_ls = list(
        filter(
            lambda v: query_version in di[object_name][v]["tags"].keys(),
            di[object_name].keys(),
        )
    )
    if len(version_ls) > 1:
        raise MinIOStorageVersioningCorruptionException(
            f"Multiple versions found for object {object_name}"
        )
    return version_ls[0] if version_ls else None


def get_remoteversion_from_bmversion(di: Dict, query_version: str) -> List:
    """
    Get the versions of all objects where the tag matches the benchmark version.

    Args:
    - di: A dictionary of objects and associated metadata.
    - query_version: The benchmark version to query.

    Returns:
    - A list of objects and associated versions where the tag matches the benchmark
    """
    summary_ls = []
    for object_name in di.keys():
        version_id = get_single_remoteversion_from_bmversion(
            di, object_name, query_version
        )
        if version_id is not None:
            summary_ls.append(
                [
                    object_name,
                    version_id,
                    di[object_name][version_id]["last_modified"],
                    di[object_name][version_id]["size"],
                    di[object_name][version_id]["etag"],
                ]
            )
    return summary_ls


def prepare_csv_remoteversion_from_bmversion(summary_ls: List) -> str:
    """
    Prepare a CSV string from a list of objects and associated versions.
    """
    outstr = "name,version_id,last_modified,size,etag\n"
    for element in summary_ls:
        outstr += f"{element[0]},{element[1]},{element[2]},{element[3]},{element[4]}\n"
    return outstr
