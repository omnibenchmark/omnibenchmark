## updated from Snakemake's test suite
##
##https://raw.githubusercontent.com/snakemake/snakemake/main/tests/tests.py
## Derivative of (c) 2012-2022 Johannes Köster johannes.koester@uni-due.com
## From Koester's license above:
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# import os
import sys
import os.path as op
# import subprocess as sp
# from pathlib import Path
from utils.run import run
import pytest

# from utils.resources import DefaultResources, GroupResources
# from utils.enums import RerunTrigger

# from snakemake.shell import shell

sys.path.insert(0, op.dirname(__file__))

# from snakemake_interface_executor_plugins.settings import DeploymentMethod

def test_conda():
    run(Snakefile = op.join("00_conda", "Snakefile"),
        produced = op.join('test0.out'),
        expected = op.join('00_conda', 'expected_results', 'test0.out'),
        method= 'conda')

def test_singularity():
    run(Snakefile = op.join("01_singularity", "Snakefile"),
        produced = op.join('test0.out'),
        expected = op.join('01_singularity', 'expected_results', 'test0.out'),
        method= 'apptainer')

def test_singularity_nonexistent():
    with pytest.raises(Exception):
        run(Snakefile = op.join("02_singularity_nonexistent", "Snakefile"),
        produced = op.join('test0.out'),
        expected = op.join('02_singularity_nonexistent', 'expected_results', 'test0.out'),
        method= 'apptainer')
