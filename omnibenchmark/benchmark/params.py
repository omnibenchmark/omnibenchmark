from collections import OrderedDict
from collections.abc import Mapping
import json
import hashlib


class Params:
    """
    Params is a data structure for workflow steps that provides ordered storage,
    serialization, and stable hashing capabilities.
    """

    def __init__(self, params=None):
        """Initialize with optional parameters dict."""
        self._params = OrderedDict()
        if params:
            # Convert input to ordered dict, sorting by keys
            if isinstance(params, (dict, OrderedDict)):
                for k in sorted(params.keys()):
                    self._params[k] = params[k]
            elif isinstance(params, Params):
                self._params = params._params.copy()

    def __eq__(self, other):
        """
        Compare for equality with another object.
        Handles comparison with other Params objects and dict-like objects.
        """
        if isinstance(other, Params):
            return self.serialize() == other.serialize()
        elif isinstance(other, Mapping):
            # Convert dict-like object to Params and compare
            return self == Params(other)
        return NotImplemented

    def __getitem__(self, key):
        """Allow dict-like access."""
        return self._params[key]

    def __setitem__(self, key, value):
        """Allow dict-like assignment."""
        self._params[key] = value

    def __str__(self):
        """Return canonical string representation."""
        return self.serialize()

    def __repr__(self):
        """Return canonical string representation."""
        return self.__str__()

    def hash(self):
        """Generate stable hash of canonical form."""
        canonical = self.serialize()
        return hashlib.sha256(canonical.encode()).hexdigest()

    def items(self):
        """Return items iterator."""
        return self._params.items()

    def update(self, other):
        """Update params from dict or another Params object."""
        if isinstance(other, Params):
            self._params.update(other._params)
        else:
            self._params.update(other)

    def serialize(self):
        """Convert to canonical JSON string."""
        # Sort keys to ensure canonical form
        return json.dumps(self._params, sort_keys=True)

    @classmethod
    def deserialize(cls, serialized):
        """Create new Params instance from serialized string."""
        params = json.loads(serialized)
        return cls(params)

    def to_cli_args(self, style="gnu"):
        """
        Convert parameters to command line arguments.

        Args:
            style (str): Either 'gnu' (--key value) or
                        'equals' (--key=value)

        Returns:
            list: List of command line argument strings
        """
        if style not in ["gnu", "equals"]:
            raise ValueError("style must be either 'gnu' or 'equals'")

        args = []
        for key, value in self._params.items():
            if isinstance(value, bool) and value:
                args.append(f"--{key}")
            else:
                if style == "gnu":
                    args.extend([f"--{key}", str(value)])
                else:  # equals style
                    args.append(f"--{key}={value}")
        return args

    @classmethod
    def from_cli_args(cls, args):
        """
        Create Params from command line arguments.
        Supports:
        - Long options: --key value and --key=value
        - Short options: -k value and -k=value
        - Implicit flags: --key and -k (True)

        Args:
            args (list): List of command line argument strings

        Returns:
            Params: New Params instance
        """
        params = {}
        i = 0
        while i < len(args):
            arg = args[i]

            # If arg contains a space inside itself ("-a 0")
            if " " in arg:
                first_part, second_part = arg.split(maxsplit=1)
                if first_part.startswith("--") or first_part.startswith("-"):
                    prefix_len = 2 if first_part.startswith("--") else 1
                    key = first_part[prefix_len:]
                    value = second_part
                    params[key] = value
                i += 1
                continue

            if arg.startswith("--") or arg.startswith("-"):
                # Remove -- or - prefix
                prefix_len = 2 if arg.startswith("--") else 1
                key = arg[prefix_len:]

                # Handle --key=value or -k=value format
                if "=" in key:
                    key, value = key.split("=", 1)
                    params[key] = value
                    i += 1
                    continue

                # Handle --key value or -k value format
                if i + 1 < len(args) and not args[i + 1].startswith("-"):
                    params[key] = args[i + 1]
                    i += 2
                else:
                    # Implicit flags (no value means True)
                    params[key] = True
                    i += 1
            else:
                # Skip non-option arguments
                i += 1

        return cls(params)
