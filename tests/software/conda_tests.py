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
import os.path as op

from utils.run import run

import pytest

sys.path.insert(0, op.dirname(__file__))

# Get the directory containing this test file
TEST_DIR = op.dirname(__file__)


@pytest.mark.skip(reason="Conda tests disabled")
def test_conda(tmp_path):
    run(
        snakefile_path=op.join(TEST_DIR, "00_conda", "Snakefile"),
        expected_path=op.join(TEST_DIR, "00_conda", "expected_results", "test0.out"),
        method="conda",
        tmp_path=tmp_path,
    )
