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

# Get the directory containing this test file
TEST_DIR = op.dirname(__file__)


# As they are, these tests are just exercising snakemake's apptainer integration, they don't belong here.
@pytest.mark.skip(reason="snakemake-apptainer tests disabled")
def test_apptainer(tmp_path):
    run(
        snakefile_path=op.join(TEST_DIR, "01_apptainer", "Snakefile"),
        expected_path=op.join(
            TEST_DIR, "01_apptainer", "expected_results", "test0.out"
        ),
        method="apptainer",
        tmp_path=tmp_path,
    )


# As they are, these tests are just exercising snakemake's apptainer integration, they don't belong here.
@pytest.mark.skip(reason="snakemake-apptainer tests disabled")
def test_apptainer_nonexistent(tmp_path):
    with pytest.raises(Exception):
        run(
            snakefile_path=op.join(TEST_DIR, "02_apptainer_nonexistent", "Snakefile"),
            expected_path=op.join(
                TEST_DIR, "02_apptainer_nonexistent", "expected_results", "test0.out"
            ),
            method="apptainer",
            tmp_path=tmp_path,
        )


@pytest.mark.skip(reason="Apptainer tests disabled")
def test_easybuild_cmd():
    check_cmd_zero_exit("eb --version")


@pytest.mark.skip(reason="Apptainer tests disabled")
def test_apptainer_cmd():
    check_cmd_zero_exit("apptainer --version")


@pytest.mark.skip(reason="Apptainer tests disabled")
def test_bash_cmd():
    check_cmd_zero_exit("bash --version")
