## updated from Snakemake's test suite
##
## https://raw.githubusercontent.com/snakemake/snakemake/main/tests/tests.py
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
import os.path as op
import pytest

from utils.run import run, check_cmd_zero_exit

sys.path.insert(0, op.dirname(__file__))


def test_singularity():
    run(
        Snakefile=op.join("01_singularity", "Snakefile"),
        produced=op.join("test0.out"),
        expected=op.join("01_singularity", "expected_results", "test0.out"),
        method="apptainer",
    )


def test_singularity_nonexistent():
    with pytest.raises(Exception):
        run(
            Snakefile=op.join("02_singularity_nonexistent", "Snakefile"),
            produced=op.join("test0.out"),
            expected=op.join(
                "02_singularity_nonexistent", "expected_results", "test0.out"
            ),
            method="apptainer",
        )


def test_easybuild_cmd():
    check_cmd_zero_exit("eb --version")


def test_singularity_cmd():
    check_cmd_zero_exit("singularity --version")


def test_bash_cmd():
    check_cmd_zero_exit("bash --version")


@pytest.mark.easybuild_toolchain
def test_easybuild_sys_toolchain_build():
    run(
        Snakefile=op.join("03_easybuild_build", "Snakefile"),
        produced=op.join("binutils-2.35.eb_ld.txt"),
        expected=op.join(
            "03_easybuild_build", "expected_results", "binutils-2.35.eb_ld.txt"
        ),
        method="apptainer",
    )
