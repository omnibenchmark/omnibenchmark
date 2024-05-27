#!/bin/bash
##
## Installs easybuild
##
## Izaskun Mallona
## 27th May 2024

## is this needed? rather add to the omni-py pyenv
pip install easybuild
eb --confighelp > $HOME/.config/easybuild/config.cfg
