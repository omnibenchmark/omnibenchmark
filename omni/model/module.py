from typing import List, Optional

import yaml
from pydantic import BaseModel, HttpUrl
from sdpx_license_list import LICENSES


class DerivedSoftware(BaseModel):
    """Model for sources that this work derives from"""

    description: str
    author: str
    url: HttpUrl
    doi: Optional[str] = None

    # TODO: should enforce URI-like, either file:// or https://
    easyconfig: Optional[str]
    conda: Optional[str]


class ModuleMetadata(BaseModel):
    """Metadata model for an Omnibenchmark Module"""

    name: str
    description: str
    author: str
    derives_from: List[DerivedSoftware]
    license: str

    def validate_license(self):
        """Validate that license is a recognized open source license"""
        if self.license not in LICENSES.keys():
            raise ValueError(
                f"License {self.license} not in recognized open source licenses"
            )

    @classmethod
    def from_yaml(cls, yaml_str: str):
        """Alternative constructor that loads from YAML string"""
        data = yaml.safe_load(yaml_str)
        return cls(**data)
