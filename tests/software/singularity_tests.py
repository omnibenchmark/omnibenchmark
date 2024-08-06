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

import sys, os  # , subprocess
import os.path as op
import pytest
from easybuild.tools.modules import modules_tool
from utils.run import run, check_cmd_zero_exit
from omni.software import easybuild_backend as eb

sys.path.insert(0, op.dirname(__file__))
# WD = op.dirname(__file__)


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


# def test_module_cmd():
#     from omni.software import easybuild_backend as easy
#     # assert easy.check_available_modules() == ['lmod', 'settarg']
#     mod_tool = modules_tool()
#     # assert len(mod_tool.available(envmodule)) != 0


def test_omni_python_import():
    import omni


def test_omni_easybuild_import():
    from omni.software import easybuild_backend as easy


# def test_env_export():
#     eb.export_lmod_env_vars()


def test_easybuild_sys_toolchain_build():
    # try:
    run(
        Snakefile=op.join("03_easybuild_build", "Snakefile"),
        produced=op.join("binutils-2.35.eb_ld.txt"),
        expected=op.join(
            "03_easybuild_build", "expected_results", "binutils-2.35.eb_ld.txt"
        ),
        method="apptainer",
    )


# def test_easybuild_build_toolchain():
#     run(Snakefile = op.join("04_easybuild_build_toolchain", "Snakefile"),
#         produced = op.join('testing_04', 'placeholder'),
#         expected = op.join('04_easybuild_build', 'expected_results', 'placeholder'),
#         method= 'apptainer')

# def test_easybuild_cli():
