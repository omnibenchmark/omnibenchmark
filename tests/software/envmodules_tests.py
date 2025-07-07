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
import os
import os.path as op
import shutil

import pytest

from utils.run import run

# Get the absolute path to the test directory
TEST_DIR = op.dirname(op.abspath(__file__))
sys.path.insert(0, TEST_DIR)


## this won't work on Mac
## Even in linux, building binutils takes too long to complete in a reasonable time (does not finish under 5 minutes)
## Perhaps we should consider testing with a shorter compilation. To be fair, I have no idea if this test
## was running before, or why did it started failing randomly after an unrelated change.
## I suspect we might have unadvertedly relying on cached results, which might have led to false positives.
@pytest.mark.skip(reason="Test takes too long to complete")
@pytest.mark.easybuild
@pytest.mark.timeout(300)  # 5 minutes timeout
def test_easybuild_sys_toolchain_build(tmp_path):
    # Copy the Snakefile to the temporary directory
    source_snakefile = op.join(TEST_DIR, "04_easybuild_build_envmodules", "Snakefile")
    dest_snakefile = op.join(tmp_path, "Snakefile")
    shutil.copy2(source_snakefile, dest_snakefile)

    # Copy the expected_results folder to the temporary directory
    source_expected_dir = op.join(
        TEST_DIR, "04_easybuild_build_envmodules", "expected_results"
    )
    dest_expected_dir = op.join(tmp_path, "expected_results")
    os.makedirs(dest_expected_dir, exist_ok=True)

    # Copy all files from the expected_results directory
    for item in os.listdir(source_expected_dir):
        source_item = op.join(source_expected_dir, item)
        dest_item = op.join(dest_expected_dir, item)
        if os.path.isfile(source_item):
            shutil.copy2(source_item, dest_item)

    # Set up paths
    produced_file = op.join(tmp_path, "binutils-2.35.eb_ld.txt")
    expected_file = op.join(
        tmp_path,
        "expected_results",
        "binutils-2.35.eb_ld.txt",
    )

    # Change to the temporary directory and run the test
    os.chdir(tmp_path)
    run(
        Snakefile=dest_snakefile,
        produced=produced_file,
        expected=expected_file,
        method="envmodules",
    )
