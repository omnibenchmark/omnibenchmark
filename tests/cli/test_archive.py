from typing import Optional


from .fixtures import minio_storage, _minio_container  # noqa: F401


def do_first_run(clisetup, file: str, cwd: Optional[str] = None):
    run1 = clisetup.call(
        [
            "run",
            "benchmark",
            "--benchmark",
            file,
        ],
        cwd=cwd,
    )
    assert run1.returncode == 0


# TODO: Uncomment when run_module functionality is restored
# def test_archive_config(minio_storage):  # noqa: F811
#     """Test archiving a benchmark configuration."""
#     expected_output = "Created archive:"

#     with OmniCLISetup() as omni:
#         # First run the benchmark to generate data
#         do_first_run(omni, str(minio_storage.benchmark_file), cwd=minio_storage.out_dir)

#         # Create a version which is needed before archiving
#         run2 = omni.call(
#             [
#                 "storage",
#                 "create-version",
#                 "--benchmark",
#                 str(minio_storage.benchmark_file),
#             ],
#             cwd=minio_storage.out_dir,
#         )
#         assert run2.returncode == 0

#         # Now run the archive command
#         run3 = omni.call(
#             [
#                 "storage",
#                 "archive",
#                 "--benchmark",
#                 str(minio_storage.benchmark_file),
#             ],
#             cwd=minio_storage.out_dir,
#         )

#         # Check the results
#         assert run3.returncode == 0
#         assert run3.stdout.startswith(expected_output)

#         # Get the benchmark name and version to check the archive file
#         benchmark = Benchmark(Path(minio_storage.benchmark_file))
#         outfile = f"{benchmark.get
