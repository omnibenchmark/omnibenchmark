from typing import List, Optional

import yaml
from pydantic import BaseModel, HttpUrl
import spdx_license_list


class DerivedSoftware(BaseModel):
    """Model for sources that this work derives from"""

    description: str
    author: str
    url: HttpUrl
    doi: Optional[str] = None

    # TODO: should enforce URI-like, either file:// or https://
    easyconfig: Optional[str] = None
    conda: Optional[str] = None


class ModuleMetadata(BaseModel):
    """Metadata model for an Omnibenchmark Module"""

    name: str
    description: str
    author: str
    derives_from: List[DerivedSoftware]
    license: str

    def validate_license(self):
        """Validate that license is a recognized open source license"""
        if self.license not in spdx_license_list.LICENSES.keys():
            raise ValueError(
                f"License {self.license} not in recognized open source licenses"
            )

    @classmethod
    def from_yaml(cls, yaml_str: str):
        """Alternative constructor that loads from YAML string"""
        data = yaml.safe_load(yaml_str)
        return cls(**data)
