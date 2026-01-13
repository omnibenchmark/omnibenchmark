#!/usr/bin/env python3
"""
Benchmark script to test omnibenchmark parallelism across different core counts.

Usage:
    python benchmark_parallelism.py Clustering_oras.yml
    python benchmark_parallelism.py Clustering_oras.yml --max-cores 32
"""

import argparse
import subprocess
import time
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime


def run_benchmark(
    yaml_file: str, cores: int, output_dir: Path, just_snakemake: bool = False
) -> dict:
    """Run ob benchmark with specified number of cores and measure time."""
    print(f"\n{'='*60}")
    print(f"Running benchmark with {cores} cores...")
    print(f"{'='*60}")

    if just_snakemake:
        # Clean only data and plotting directories (we're already in out/)
        print("  Cleaning 'data' and 'plotting' directories...")
        subprocess.run(["rm", "-rf", "data"], check=False)
        subprocess.run(["rm", "-rf", "plotting"], check=False)
    else:
        # Clean previous outputs before each run
        print("  Cleaning 'out' directory...")
        subprocess.run(["rm", "-rf", "out"], check=False)

    start_time = time.time()

    # Prepare log file for this run
    log_file = output_dir / f"run_cores_{cores}.log"

    print(f"  Logging output to: {log_file}")

    # Run without check=True to allow non-zero exit codes (e.g., metric collector failures)
    # Use -k flag to continue on errors
    # Redirect stdin to /dev/null to prevent interactive prompts from blocking
    import sys

    if just_snakemake:
        # Run bare snakemake (already in out/ directory)
        cmd = ["snakemake", "--use-conda", "--cores", str(cores), "-k"]
        cwd = None
    else:
        # Run ob run command
        cmd = ["ob", "run", yaml_file, "--cores", str(cores), "-k", "--yes"]
        cwd = None

    with open(log_file, "w", buffering=1) as log_f, open("/dev/null", "r") as devnull:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # Merge stderr into stdout
            stdin=devnull,  # Prevent reading from stdin
            text=True,
            bufsize=1,  # Line buffered
            cwd=cwd,
        )

        # Read line by line and output to both console and log
        for line in iter(process.stdout.readline, ""):
            if line:
                sys.stdout.write(line)
                sys.stdout.flush()
                log_f.write(line)
            else:
                break

        process.wait()
        exit_code = process.returncode

    elapsed_time = time.time() - start_time

    # Always consider it a success for timing purposes
    # (we're measuring parallelism, not correctness)
    success = True

    if exit_code == 0:
        print(f"\n✓ Completed in {elapsed_time:.2f} seconds (exit code: 0)")
    else:
        print(
            f"\n⚠ Completed in {elapsed_time:.2f} seconds (exit code: {exit_code}, but continuing)"
        )

    return {
        "cores": cores,
        "time": elapsed_time,
        "success": success,
        "exit_code": exit_code,
        "log_file": str(log_file),
        "timestamp": datetime.now().isoformat(),
    }


def fit_polynomial_and_analyze(cores_list, times_list, degree=2):
    """Fit Amdahl's law or power law to timing data."""
    from scipy.optimize import curve_fit

    cores_array = np.array(cores_list, dtype=float)
    times_array = np.array(times_list, dtype=float)

    # Find baseline (1-core time)
    baseline_idx = cores_list.index(min(cores_list))
    t1 = times_array[baseline_idx]

    # Amdahl's law: T(n) = T(1) * (s + (1-s)/n)
    def amdahl(n, s):
        return t1 * (s + (1 - s) / n)

    # Power law: T(n) = a * n^b (common for sub-linear speedup)
    def power_law(n, a, b):
        return a * np.power(n, b)

    try:
        # Fit Amdahl's law (s should be between 0 and 1)
        popt_amdahl, _ = curve_fit(
            amdahl, cores_array, times_array, bounds=(0, 1), p0=[0.1]
        )
        fitted_amdahl = amdahl(cores_array, *popt_amdahl)
        ss_res_amdahl = np.sum((times_array - fitted_amdahl) ** 2)

        serial_fraction = popt_amdahl[0]
    except (RuntimeError, ValueError, TypeError):
        serial_fraction = None
        ss_res_amdahl = float("inf")

    try:
        # Fit power law
        popt_power, _ = curve_fit(
            power_law, cores_array, times_array, p0=[t1, -0.5], maxfev=10000
        )
        fitted_power = power_law(cores_array, *popt_power)
        ss_res_power = np.sum((times_array - fitted_power) ** 2)
    except (RuntimeError, ValueError, TypeError):
        popt_power = None
        ss_res_power = float("inf")

    # Choose best fit
    if ss_res_amdahl < ss_res_power:
        fitted_times = fitted_amdahl
        model_type = "Amdahl"
        model_params = {"serial_fraction": serial_fraction}

        def fit_func(n):
            return amdahl(n, serial_fraction)
    else:
        fitted_times = fitted_power
        model_type = "Power"
        model_params = {"a": popt_power[0], "b": popt_power[1]}

        def fit_func(n):
            return power_law(n, *popt_power)

    # Calculate R-squared
    residuals = times_array - fitted_times
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((times_array - np.mean(times_array)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

    # Calculate RMSE
    rmse = np.sqrt(np.mean(residuals**2))

    # Calculate relative dispersion
    mean_time = np.mean(times_array)
    dispersion = rmse / mean_time if mean_time != 0 else float("inf")

    return {
        "model_type": model_type,
        "model_params": model_params,
        "fit_function": fit_func,
        "polynomial": fit_func,  # Keep for backward compat
        "fitted_times": fitted_times,
        "r_squared": r_squared,
        "rmse": rmse,
        "dispersion": dispersion,
        "residuals": residuals,
    }


def calculate_speedup_efficiency(cores_list, times_list):
    """Calculate speedup and parallel efficiency metrics relative to 1 core."""
    # Find the 1-core baseline
    baseline_idx = cores_list.index(min(cores_list))
    baseline_time = times_list[baseline_idx]
    baseline_cores = cores_list[baseline_idx]

    speedups = []
    efficiencies = []

    for cores, time_val in zip(cores_list, times_list):
        # Speedup: T(baseline) / T(n)
        speedup = baseline_time / time_val
        speedups.append(speedup)

        # Efficiency: Speedup / cores * 100%
        ideal_speedup = cores / baseline_cores
        efficiency = (speedup / ideal_speedup) * 100 if ideal_speedup > 0 else 0
        efficiencies.append(efficiency)

    return speedups, efficiencies


def plot_results(results, output_dir: Path, fit_results):
    """Create comprehensive plots of benchmark results."""
    cores_list = [r["cores"] for r in results if r["success"]]
    times_list = [r["time"] for r in results if r["success"]]

    if len(cores_list) < 2:
        print("Not enough successful runs to plot")
        return

    speedups, efficiencies = calculate_speedup_efficiency(cores_list, times_list)

    # Create figure with 2 subplots (1 row, 2 columns)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Omnibenchmark Parallelism Analysis", fontsize=16, fontweight="bold")

    # Plot 1: Execution time vs cores
    ax1 = axes[0]
    ax1.plot(cores_list, times_list, "o-", label="Measured", markersize=8, linewidth=2)

    # Plot model fit
    cores_smooth = np.linspace(min(cores_list), max(cores_list), 100)
    fitted_smooth = fit_results["polynomial"](cores_smooth)
    model_label = f'{fit_results.get("model_type", "Model")} fit (R²={fit_results["r_squared"]:.3f})'
    ax1.plot(
        cores_smooth, fitted_smooth, "--", label=model_label, linewidth=2, alpha=0.7
    )

    ax1.set_xlabel("Number of Cores", fontsize=12)
    ax1.set_ylabel("Execution Time (seconds)", fontsize=12)
    ax1.set_title("Execution Time vs Cores", fontsize=13, fontweight="bold")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Speedup vs cores (relative to 1 core)
    ax2 = axes[1]
    ideal_speedup = [c for c in cores_list]  # Ideal is linear with cores
    ax2.plot(
        cores_list,
        speedups,
        "o-",
        label="Actual Speedup",
        markersize=8,
        linewidth=2,
        color="green",
    )
    ax2.plot(
        cores_list,
        ideal_speedup,
        "--",
        label="Ideal (Linear)",
        linewidth=2,
        alpha=0.7,
        color="gray",
    )
    ax2.set_xlabel("Number of Cores", fontsize=12)
    ax2.set_ylabel("Speedup (relative to 1 core)", fontsize=12)
    ax2.set_title("Speedup vs Cores", fontsize=13, fontweight="bold")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()

    # Save plot
    plot_file = output_dir / "parallelism_benchmark.png"
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    print(f"\n📊 Plot saved to: {plot_file}")

    plt.show()


def print_summary(results, fit_results):
    """Print summary statistics."""
    cores_list = [r["cores"] for r in results if r["success"]]
    times_list = [r["time"] for r in results if r["success"]]

    if len(cores_list) < 2:
        print("\n⚠️  Not enough successful runs for analysis")
        return

    speedups, efficiencies = calculate_speedup_efficiency(cores_list, times_list)

    print("\n" + "=" * 70)
    print("BENCHMARK SUMMARY")
    print("=" * 70)

    print("\nExecution Times:")
    print(f"{'Cores':<10} {'Time (s)':<15} {'Speedup':<12} {'Efficiency (%)':<15}")
    print("-" * 70)
    for cores, time_val, speedup, eff in zip(
        cores_list, times_list, speedups, efficiencies
    ):
        print(f"{cores:<10} {time_val:<15.2f} {speedup:<12.2f} {eff:<15.1f}")

    print("\n" + "-" * 70)
    print(f"{fit_results.get('model_type', 'Model')} Fit Analysis:")
    print("-" * 70)
    print(f"Model type:          {fit_results.get('model_type', 'Unknown')}")
    print(f"R² (goodness of fit): {fit_results['r_squared']:.4f}")
    print(f"RMSE:                {fit_results['rmse']:.2f} seconds")
    print(f"Relative Dispersion: {fit_results['dispersion']:.2%}")

    if fit_results.get("model_type") == "Amdahl":
        s = fit_results["model_params"]["serial_fraction"]
        print("\nAmdahl's Law: T(n) = T(1) × (s + (1-s)/n)")
        print(f"  Serial fraction (s):   {s:.4f} ({s*100:.2f}%)")
        print(f"  Parallel fraction:     {1-s:.4f} ({(1-s)*100:.2f}%)")
        max_speedup = 1 / s if s > 0 else float("inf")
        print(f"  Max theoretical speedup: {max_speedup:.1f}x")
    elif fit_results.get("model_type") == "Power":
        a = fit_results["model_params"]["a"]
        b = fit_results["model_params"]["b"]
        print("\nPower Law: T(n) = a × n^b")
        print(f"  Coefficient a: {a:.2f}")
        print(f"  Exponent b:    {b:.4f}")
        if b < 0:
            print("  (Negative exponent indicates speedup with more cores)")

    print("\n" + "-" * 70)
    print("Performance Indicators:")
    print("-" * 70)

    # Best efficiency
    max_eff_idx = efficiencies.index(max(efficiencies))
    print(
        f"Best efficiency:     {efficiencies[max_eff_idx]:.1f}% at {cores_list[max_eff_idx]} cores"
    )

    # Average efficiency
    avg_efficiency = np.mean(efficiencies)
    print(f"Average efficiency:  {avg_efficiency:.1f}%")

    # Time saved at max cores (compare min cores vs max cores)
    min_cores_idx = cores_list.index(min(cores_list))
    max_cores_idx = cores_list.index(max(cores_list))
    time_1core = times_list[min_cores_idx]
    time_maxcore = times_list[max_cores_idx]
    time_saved = time_1core - time_maxcore
    time_saved_pct = (time_saved / time_1core) * 100
    print(
        f"\nTime saved (1→{max(cores_list)} cores): {time_saved:.1f}s ({time_saved_pct:.1f}%)"
    )

    print("=" * 70 + "\n")


def plot_only(results_file: Path, output_dir: Path):
    """Load existing results and generate plots without re-running benchmarks."""
    print(f"📊 Loading results from: {results_file}")

    with open(results_file, "r") as f:
        results = json.load(f)

    # Filter successful results
    successful_results = [r for r in results if r.get("success", True)]

    if len(successful_results) < 2:
        print("⚠️  Need at least 2 successful runs for analysis")
        return 1

    # Analyze results
    cores_list = [r["cores"] for r in successful_results]
    times_list = [r["time"] for r in successful_results]

    fit_results = fit_polynomial_and_analyze(cores_list, times_list)

    # Print summary
    print_summary(successful_results, fit_results)

    # Plot results
    plot_results(successful_results, output_dir, fit_results)

    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark omnibenchmark parallelism across different core counts"
    )
    parser.add_argument("yaml_file", nargs="?", help="Benchmark YAML file to run")
    parser.add_argument(
        "--cores",
        type=int,
        nargs="+",
        default=[1, 2, 4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 120],
        help="List of core counts to test (default: 1 2 4 8 12 16 20 24 28 32 40 48 56 64 80 96 120)",
    )
    parser.add_argument(
        "--max-cores",
        type=int,
        help="Maximum cores to test (overrides --cores with automatic sequence)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory for results (default: current directory)",
    )
    parser.add_argument(
        "--skip-failed", action="store_true", help="Continue even if a run fails"
    )
    parser.add_argument(
        "--plot-only",
        type=Path,
        metavar="RESULTS_JSON",
        help="Plot results from existing benchmark_results.json without re-running",
    )
    parser.add_argument(
        "--just-snakemake",
        action="store_true",
        help="Run bare snakemake from out/ directory instead of 'ob run' (requires pre-generated Snakefile)",
    )

    args = parser.parse_args()

    # Handle plot-only mode
    if args.plot_only:
        return plot_only(args.plot_only, args.output_dir)

    if not args.yaml_file and not args.just_snakemake:
        parser.error(
            "yaml_file is required unless using --plot-only or --just-snakemake"
        )

    # Determine core counts to test
    if args.max_cores:
        cores_to_test = [
            1,
            2,
            4,
            8,
            12,
            16,
            20,
            24,
            28,
            32,
            40,
            48,
            56,
            64,
            80,
            96,
            120,
        ]
        cores_to_test = [c for c in cores_to_test if c <= args.max_cores]
    else:
        cores_to_test = sorted(args.cores)

    # Reverse the order to start with highest parallelism and go down
    cores_to_test = list(reversed(cores_to_test))

    print("\n🚀 Starting parallelism benchmark")
    if args.just_snakemake:
        print("Mode: bare snakemake")
    else:
        print(f"Benchmark file: {args.yaml_file}")
    print(f"Testing with cores: {cores_to_test}")
    print(f"Output directory: {args.output_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Run benchmarks
    results = []
    for cores in cores_to_test:
        result = run_benchmark(
            args.yaml_file, cores, args.output_dir, args.just_snakemake
        )
        results.append(result)

        if not result["success"] and not args.skip_failed:
            print(f"\n❌ Benchmark failed with {cores} cores. Stopping.")
            break

    # Save raw results
    results_file = args.output_dir / "benchmark_results.json"
    with open(results_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n💾 Results saved to: {results_file}")

    # Filter successful results for analysis
    successful_results = [r for r in results if r["success"]]

    if len(successful_results) < 2:
        print("\n⚠️  Need at least 2 successful runs for analysis")
        return 1

    # Analyze results
    cores_list = [r["cores"] for r in successful_results]
    times_list = [r["time"] for r in successful_results]

    fit_results = fit_polynomial_and_analyze(cores_list, times_list, degree=2)

    # Print summary
    print_summary(successful_results, fit_results)

    # Plot results
    plot_results(successful_results, args.output_dir, fit_results)

    return 0


if __name__ == "__main__":
    exit(main())
