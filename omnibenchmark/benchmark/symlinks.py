import os
from pathlib import Path
from filelock import FileLock

from .params import Params

UNSAFE_CHARS = {
    "/": "_",
    "\\": "_",
    ":": "_",
    "*": "_",
    "?": "_",
    '"': "_",
    "<": "_",
    ">": "_",
    "|": "_",
    " ": "_",
}


class SymlinkManager:
    """
    Manages a two-level storage system for parameters with human-readable access.

    This class implements a storage pattern where parameters are stored in a content-addressable
    storage (using their hash as identifier) while maintaining human-readable symbolic links
    for easier browsing and access.

    The storage structure is:
        {base_dir}/               # Human-readable symlinks
            a_1-b_2-c_3/          # -> {.hash123}
            x_1-y_2/              # -> {.hash456}
            {.hash123}/           # Actual storage
                parameters.json
            {.hash456}/
                parameters.json

    Args:
        base_dir (str): Base directory for storage

    Example:
        >>> manager = SymlinkManager("./data")
        >>> info = manager.store(params)
        >>> print(f"Stored in {info['human']}, hash: {info['folder']}")
        >>> retrieved_params = manager.get_params("param_a1_b2_c3")
    """

    def __init__(self, base_dir: str | Path = "."):
        self.base_dir: Path = Path(base_dir)
        self.lock = FileLock(self.base_dir / "storage.lock")

    def _sanitize_name(self, name):
        """Basic character sanitization for filesystem safety."""
        for unsafe, safe in UNSAFE_CHARS.items():
            name = name.replace(unsafe, safe)
        return name

    def _make_human_name(self, params, max_len=255):
        """
        Create human-readable name, falling back appending short hash if too long.

        Args:
            params: Params object to create name from
            max_len (int): Maximum length of the folder name.

        Returns:
            str: Human-readable name within filesystem constraints
        """
        parts = [f"{k}-{v}" for k, v in params.items()]
        name = "_".join(parts)
        name = self._sanitize_name(name)

        # If too long, append a short 8 characters hash
        if len(name) > max_len:
            short_hash = params.hash()[:8]
            name = name[: max_len - 9] + "_" + short_hash
        return name

    def _write_params(self, params, hash_path, human_path):
        params_hash = params.hash()
        params_file = hash_path / "parameters.json"
        with open(params_file, "w") as f:
            f.write(params.serialize())

        params_dict_tsv = hash_path.parent / "parameters_dict.tsv"
        if os.path.exists(params_dict_tsv):
            with open(params_dict_tsv, "r") as file:
                existing_lines = file.readlines()

            # If the params_hash exists, skip writing
            if any(params_hash in line for line in existing_lines):
                return
        else:
            with open(params_dict_tsv, "w") as file:
                file.write("id\tbase_path\talias_path\tparameters\n")

        # If not found, write the new entry to parameters_dict.tsv
        with open(params_dict_tsv, "a") as file:
            file.write(
                f"{params_hash}\t{hash_path}\t{human_path}\t{params.serialize()}\n"
            )

    def store(self, params):
        """
        Store parameters and create human-readable symlink.

        Args:
            params: Params object to store

        Returns:
            dict: Information about stored params including:
                - folder: Hash-based folder name
                - params: Original params object
                - human: Human-readable name
        """
        if params:
            with self.lock:
                hash_id = params.hash()
                hash_folder = f".{hash_id}"
                human_name = self._make_human_name(params)

                # Create hash directory and store params
                hash_path = self.base_dir / hash_folder
                hash_path.mkdir(exist_ok=True)

                # Create/update symlink
                symlink_path = self.base_dir / human_name
                if symlink_path.exists() or symlink_path.is_symlink():
                    symlink_path.unlink()

                # Create relative path from visible_dir to hidden_dir
                real_target = self.base_dir / hash_folder
                relative_path = os.path.relpath(real_target, symlink_path.parent)
                symlink_path.symlink_to(relative_path)

                self._write_params(params, hash_path, symlink_path)

                return {"folder": hash_folder, "params": params, "human": human_name}

    def get_params(self, human_name):
        """
        Retrieve parameters from human-readable name.
        """
        symlink_path = self.base_dir / human_name
        if not symlink_path.exists():
            raise FileNotFoundError(f"No parameters found for name: {human_name}")

        # Follow symlink to get the actual parameters.json file
        target = symlink_path.readlink()
        params_file = symlink_path.parent / target / "parameters.json"

        if not params_file.exists():
            raise ValueError(
                f"Parameters file missing for: {human_name} (looking in {params_file})"
            )

        with open(params_file) as f:
            return Params.deserialize(f.read())
