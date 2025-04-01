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
            if style == "gnu":
                args.extend([f"--{key}", str(value)])
            else:  # equals style
                args.append(f"--{key}={value}")
        return args

    @classmethod
    def from_cli_args(cls, args):
        """
        Create Params from command line arguments.
        Accepts both --key value and --key=value formats.

        Args:
            args (list): List of command line argument strings

        Returns:
            Params: New Params instance
        """
        params = {}
        i = 0
        while i < len(args):
            arg = args[i]
            if not arg.startswith("--"):
                i += 1
                continue

            key = arg[2:]  # Remove '--' prefix

            # Handle --key=value format
            if "=" in key:
                key, value = key.split("=", 1)
                params[key] = value
                i += 1
                continue

            # Handle --key value format
            if i + 1 < len(args) and not args[i + 1].startswith("--"):
                params[key] = args[i + 1]
                i += 2
            else:
                # Handle flag-like parameters (no value)
                params[key] = True
                i += 1

        return cls(params)
