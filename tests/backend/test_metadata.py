"""Unit tests for omnibenchmark/backend/_metadata.py."""

from omnibenchmark.backend._metadata import save_metadata

from tests.backend.test_snakemake_gen import _make_node


class TestSaveMetadata:
    def test_creates_metadata_dir(self, tmp_path):
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        save_metadata(benchmark_yaml, out_dir, [_make_node()])
        assert (out_dir / ".metadata").is_dir()
        assert (out_dir / ".metadata" / "benchmark.yaml").exists()
        assert (out_dir / ".metadata" / "modules.txt").exists()

    def test_modules_txt_contains_module_info(self, tmp_path):
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        save_metadata(benchmark_yaml, out_dir, [_make_node()])
        modules_txt = (out_dir / ".metadata" / "modules.txt").read_text()
        assert "https://github.com/example/repo" in modules_txt
        assert "abc1234" in modules_txt

    def test_duplicate_modules_deduplicated(self, tmp_path):
        benchmark_yaml = tmp_path / "benchmark.yaml"
        benchmark_yaml.write_text("name: test")
        out_dir = tmp_path / "out"
        out_dir.mkdir()
        node1 = _make_node(node_id="stage1-mod1-p1", param_id="p1")
        node2 = _make_node(node_id="stage1-mod1-p2", param_id="p2")
        save_metadata(benchmark_yaml, out_dir, [node1, node2])
        modules_txt = (out_dir / ".metadata" / "modules.txt").read_text()
        assert modules_txt.count("https://github.com/example/repo") == 1
