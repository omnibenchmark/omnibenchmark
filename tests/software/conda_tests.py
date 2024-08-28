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

from utils.run import run, check_cmd_zero_exit

sys.path.insert(0, op.dirname(__file__))
# WD = op.dirname(__file__)


def test_conda():
    run(
        Snakefile=op.join("00_conda", "Snakefile"),
        produced=op.join("test0.out"),
        expected=op.join("00_conda", "expected_results", "test0.out"),
        method="conda",
    )


def test_omni_python_import():
    import omni


def test_omni_easybuild_import():
    from omni.software import easybuild_backend as easy
