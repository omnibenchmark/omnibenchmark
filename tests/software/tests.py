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

import sys
import os
import os.path as op
import pytest

from tests.software.utils.run import *

sys.path.insert(0, op.dirname(__file__))
WD = op.dirname(__file__)

def test_conda():
    run(Snakefile = op.join(WD, "00_conda", "Snakefile"),
        produced = op.join('test0.out'),
        expected = op.join(WD, '00_conda', 'expected_results', 'test0.out'),
        method= 'conda')

def test_singularity():
    run(Snakefile = op.join(WD, "01_singularity", "Snakefile"),
        produced = op.join( 'test0.out'),
        expected = op.join(WD, '01_singularity', 'expected_results', 'test0.out'),
        method= 'apptainer')

def test_singularity_nonexistent():
    with pytest.raises(Exception):
        run(Snakefile = op.join(WD, "02_singularity_nonexistent", "Snakefile"),
        produced = op.join('test0.out'),
        expected = op.join(WD, '02_singularity_nonexistent', 'expected_results', 'test0.out'),
        method= 'apptainer')

def test_easybuild_cmd():
    run_subprocess('eb --version', 'This is EasyBuild')

def test_easybuild_build():
    run(Snakefile = op.join(WD, "03_easybuild_build", "Snakefile"),
        produced = op.join('test_03', 'binutils-2.35.eb_hello.txt'),
        expected = op.join(WD, '03_easybuild_build', 'expected_results', 'binutils-2.35.eb_hello.txt'),
        method= 'apptainer')
