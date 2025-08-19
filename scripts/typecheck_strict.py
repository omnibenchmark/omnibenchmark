#!/usr/bin/env python3
"""Run strict type checking on specific modules."""

import subprocess
import sys
import shutil
from pathlib import Path
from typing import List, Tuple

# Modules that should have strict type checking
STRICT_MODULES = [
    "omnibenchmark/dag",
    "omnibenchmark/model",
    "omnibenchmark/benchmark",
    "tests/dag",
]


def check_tool_available(tool: str) -> bool:
    """Check if a type checking tool is available."""
    if tool == "ty":
        # Check if ty is available via uvx
        if shutil.which("uvx"):
            try:
                result = subprocess.run(
                    ["uvx", "ty", "--version"],
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=5,
                )
                return result.returncode == 0
            except (subprocess.TimeoutExpired, Exception):
                return False
        # Fall back to direct ty installation
        return shutil.which(tool) is not None
    return shutil.which(tool) is not None


def run_ty_on_module(module_path: str) -> Tuple[bool, str]:
    """Run ty on a specific module.

    Args:
        module_path: Path to the module to check

    Returns:
        Tuple of (success, output)
    """
    # Use uvx if available, otherwise direct ty
    if shutil.which("uvx"):
        cmd = ["uvx", "ty", "check", module_path]
    else:
        cmd = ["ty", "check", module_path]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=Path(__file__).parent.parent,
            timeout=30,  # Add timeout for ty
        )
        success = result.returncode == 0
        output = result.stdout if result.stdout else result.stderr

        # Check for ty pre-release warnings/errors
        if "pre-release software" in output or "not ready for production" in output:
            return False, f"ty pre-release error: {output}"

        # ty output is usually more concise
        if success and output.strip():
            lines = output.strip().split("\n")
            if lines and ("Success" in lines[-1] or "✓" in output):
                output = f"ty: {lines[-1]}"

        return success, output
    except subprocess.TimeoutExpired:
        return False, "ty timed out (ty is pre-release software)"
    except Exception as e:
        return False, f"ty failed (pre-release): {e}"


def run_pyright_on_module(module_path: str) -> Tuple[bool, str]:
    """Run pyright on a specific module with strict checking.

    Args:
        module_path: Path to the module to check

    Returns:
        Tuple of (success, output)
    """
    config_file = Path(__file__).parent.parent / "pyrightconfig.strict.json"
    cmd = [sys.executable, "-m", "pyright", "--project", str(config_file), module_path]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
            cwd=Path(__file__).parent.parent,
        )
        success = result.returncode == 0
        output = result.stdout if result.stdout else result.stderr

        return success, output
    except Exception as e:
        return False, f"Error running pyright: {e}"


def main() -> int:
    """Run strict type checking on all specified modules."""
    # Determine which tool to use - try ty first for speed
    use_ty = check_tool_available("ty")
    has_pyright = check_tool_available("pyright")

    if not use_ty and not has_pyright:
        print(
            "❌ Neither ty nor pyright is available. Please install pyright: pip install pyright"
        )
        return 1

    # Start with ty but fallback to pyright if ty fails
    tool_name = "ty (fast, with pyright fallback)" if use_ty else "pyright"

    print(f"🔍 Running strict type checking on core modules using {tool_name}...")
    print("=" * 70)

    if use_ty:
        print("\n⚠️  ty is pre-release software - will fallback to pyright on errors")
        print("Using ty.toml configuration for fast type checking.")
    else:
        print("\nUsing pyrightconfig.strict.json configuration.")

    print("⚠️  During gradual migration, failures are expected and won't block CI.\n")

    all_passed = True
    results: List[Tuple[str, bool, str]] = []
    ty_failures = 0

    for module in STRICT_MODULES:
        module_path = Path(module)
        if not module_path.exists():
            print(f"⚠️  Skipping {module} (path does not exist)")
            continue

        print(f"Checking {module}...", end=" ", flush=True)

        success = False
        output = ""

        # Try ty first if available
        if use_ty:
            success, output = run_ty_on_module(module)
            if not success and ("pre-release" in output or "ty failed" in output):
                ty_failures += 1
                print("ty failed, trying pyright...", end=" ", flush=True)
                success, output = run_pyright_on_module(module)
        else:
            success, output = run_pyright_on_module(module)

        results.append((module, success, output))

        if success:
            # Extract summary from output
            lines = output.strip().split("\n")
            summary = lines[-1] if lines else "No output"

            # Check for common success patterns
            if any(pattern in summary for pattern in ["0 errors", "Success", "✓"]):
                print(f"✅ PASS ({summary})")
            else:
                print("✅ PASS")
        else:
            print("❌ FAIL (expected during migration)")
            all_passed = False

    # Print detailed results only for failures
    failures = [(m, o) for m, s, o in results if not s]
    if failures:
        print("\n" + "=" * 70)
        print("FAILED MODULES - DETAILED OUTPUT")
        print("=" * 70)

        for module, output in failures:
            print(f"\n--- {module} ---")
            print(output)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    passed = sum(1 for _, success, _ in results if success)
    total = len(results)

    print(f"\nModules checked: {total}")
    print(f"Passed: {passed}")
    print(f"Failed: {total - passed}")

    if ty_failures > 0:
        print(
            f"\n⚠️  ty failed {ty_failures} times (pre-release software) - used pyright fallback"
        )

    if all_passed:
        print("\n🎉 All modules passed strict type checking!")
        print(
            "\nThese modules are configured with strict type checking using pyrightconfig.strict.json."
        )
        print("This ensures:")
        print("- All functions have type annotations")
        print("- All variables have clear types")
        print("- No implicit Any types")
        print("- Strict null/None checking")
        return 0
    else:
        print("\n⚠️  Some modules have type errors (expected during gradual migration).")
        print("\n📝 To gradually fix type errors:")
        print("1. Add type annotations to all function parameters and return values")
        print("2. Ensure all variables have explicit types or can be inferred")
        print("3. Use 'typing' module for complex types (List, Dict, Optional, etc.)")
        print("4. Avoid using Any type unless absolutely necessary")
        print("5. Handle None values explicitly with Optional[T] or Union[T, None]")
        print("\n💡 This check will become required once migration is complete.")
        print("\n🔧 Note: ty is pre-release software - use pyright for stable checking")
        # Return 0 during migration to not block CI
        return 0


if __name__ == "__main__":
    sys.exit(main())
