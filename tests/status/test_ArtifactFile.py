from omnibenchmark.benchmark.execution_path import (
    n2f,
    LocalArtifactFile,
    ArtifactFileConfig,
)

import pytest


@pytest.mark.short
def test_n2f():
    # Test None converts to False
    assert n2f(None) is False

    # Test True stays True
    assert n2f(True) is True

    # Test False converts to False
    assert n2f(False) is False


@pytest.mark.short
def test_local_artifact_file_exists(tmp_path):
    # Create a temporary file to test LocalArtifactFile
    test_file = tmp_path / "test_file.txt"
    test_file.write_text("This is a test file.")

    artifact_file = LocalArtifactFile(file_path=test_file)

    # Test that the file exists
    assert artifact_file.exists() is True


@pytest.mark.short
def test_local_artifact_file_is_empty(tmp_path):
    # Create a temporary empty file to test LocalArtifactFile
    empty_file = tmp_path / "empty_file.txt"
    empty_file.touch()

    artifact_file = LocalArtifactFile(file_path=empty_file)

    # Test that the file is empty
    assert artifact_file.is_empty() is True


@pytest.mark.short
def test_local_artifact_file_get_timestamp(tmp_path):
    # Create a temporary file to test LocalArtifactFile
    test_file = tmp_path / "timestamp_file.txt"
    test_file.write_text("This is a test file.")

    artifact_file = LocalArtifactFile(file_path=test_file)

    # Test that the timestamp is not None
    timestamp = artifact_file.get_timestamp()
    assert timestamp is not None
    assert isinstance(timestamp, float)


@pytest.mark.short
def test_artifact_file_config_default():
    # Test that default implementation is LocalArtifactFile
    ArtifactFileConfig._implementation = None  # Reset to default
    implementation = ArtifactFileConfig.get_implementation()
    assert implementation == LocalArtifactFile


@pytest.mark.short
def test_artifact_file_config_set_implementation():
    # Create a mock implementation class
    class MockArtifactFile:
        def __init__(self, file_path):
            self.file_path = file_path

    # Set the implementation
    ArtifactFileConfig.set_implementation(MockArtifactFile)

    # Verify it was set correctly
    implementation = ArtifactFileConfig.get_implementation()
    assert implementation == MockArtifactFile

    # Reset to default
    ArtifactFileConfig._implementation = None
    assert ArtifactFileConfig.get_implementation() == LocalArtifactFile


@pytest.mark.short
def test_create_artifact_file_with_default(tmp_path):
    # Ensure default implementation is set
    ArtifactFileConfig._implementation = None

    test_file = tmp_path / "test_file.txt"
    test_file.write_text("Test content")

    # Create artifact file using factory function
    artifact = ArtifactFileConfig.get_implementation()(test_file)

    # Verify it's a LocalArtifactFile instance
    assert isinstance(artifact, LocalArtifactFile)
    assert artifact.exists() is True


@pytest.mark.short
def test_create_artifact_file_with_custom_implementation(tmp_path):
    # Create a custom implementation
    class CustomArtifactFile:
        def __init__(self, file_path):
            self.file_path = file_path
            self.custom_flag = True

    # Set custom implementation
    ArtifactFileConfig.set_implementation(CustomArtifactFile)

    test_file = tmp_path / "test_file.txt"

    # Create artifact file using factory function
    artifact = ArtifactFileConfig.get_implementation()(test_file)

    # Verify it's a CustomArtifactFile instance
    assert isinstance(artifact, CustomArtifactFile)
    assert hasattr(artifact, "custom_flag")
    assert artifact.custom_flag is True

    # Reset to default
    ArtifactFileConfig._implementation = None
