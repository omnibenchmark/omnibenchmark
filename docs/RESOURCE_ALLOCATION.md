# Resource Allocation for Parallel Execution

This document describes how to configure resource requirements in OmniBenchmark to optimize parallel execution with Snakemake.

## Overview

OmniBenchmark allows you to specify resource requirements (cores, memory, disk, runtime) at both the **stage level** and **module level**. These hints are translated into Snakemake `resources:` directives, enabling efficient parallel scheduling across your available compute resources.

### Why Resource Allocation Matters

When running benchmarks on high-performance systems with many cores (e.g., 100+ cores), you want to:

1. **Maximize throughput**: Run as many tasks in parallel as possible
2. **Prevent resource contention**: Ensure tasks don't compete for limited resources
3. **Optimize scheduling**: Give Snakemake accurate hints about task requirements

Without resource hints, Snakemake treats all tasks equally and may:
- Over-allocate resources (running too few tasks in parallel)
- Under-allocate resources (causing tasks to compete and slow down)

## Resource Hierarchy

Resources are resolved using a **most-specific-wins** hierarchy:

1. **Module-level resources** (highest priority)
   - Defined on individual modules
   - Overrides stage-level settings
   - Use for modules with special requirements

2. **Stage-level resources** (medium priority)
   - Defined on stages
   - Applies to all modules in that stage
   - Use as sensible defaults for a stage

3. **Global default** (lowest priority)
   - Default: `cores=2`
   - Applied when no resources are specified
   - Conservative default for safety

## Syntax

### Stage-Level Resources

Apply to all modules in a stage:

```yaml
stages:
  - id: clustering
    resources:
      cores: 20        # Cores per task
      mem_mb: 16000    # Memory in MB
      disk_mb: 5000    # Disk space in MB (optional)
      runtime: 120     # Expected runtime in minutes (optional)
    modules:
      - id: fastcluster
        # Inherits stage resources (20 cores, 16GB memory)
        ...
      - id: sklearn
        # Also inherits stage resources
        ...
```

### Module-Level Resources

Override stage defaults for specific modules:

```yaml
stages:
  - id: clustering
    resources:
      cores: 20        # Stage default
      mem_mb: 16000
    modules:
      - id: fastcluster
        # Uses stage default (20 cores)
        ...
      
      - id: lightweight_method
        # Override: this method only needs 4 cores
        resources:
          cores: 4
          mem_mb: 4000
        ...
```

### Metric Collector Resources

Collectors can also specify resources:

```yaml
metric_collectors:
  - id: plotting
    resources:
      cores: 4
      mem_mb: 8000
    ...
```

## Available Resource Types

| Resource | Description | Unit | Required |
|----------|-------------|------|----------|
| `cores` | CPU cores to allocate | integer | No* |
| `mem_mb` | Memory to allocate | megabytes | No |
| `disk_mb` | Disk space required | megabytes | No |
| `runtime` | Expected runtime | minutes | No |

*At least one resource must be specified if using a `resources:` block.

## Usage Examples

### Example 1: Simple Core Allocation

```yaml
stages:
  - id: data
    # No resources: uses default (2 cores)
    modules: [...]
    
  - id: analysis
    resources:
      cores: 20  # Up to 20 cores per task
    modules: [...]
```

### Example 2: Mixed Requirements

```yaml
stages:
  - id: preprocessing
    resources:
      cores: 8
      mem_mb: 8000
    modules:
      - id: normalize
        # Uses 8 cores, 8GB (inherited)
        ...
      
      - id: filter_outliers
        # Override: needs more memory
        resources:
          cores: 4
          mem_mb: 16000
        ...
```

### Example 3: Resource Pools

With 100 available cores, you can create effective "resource pools":

```yaml
stages:
  - id: lightweight_tasks
    resources:
      cores: 2
    # Can run up to 50 tasks in parallel (100 / 2 = 50)
    modules: [...]
    
  - id: medium_tasks
    resources:
      cores: 10
    # Can run up to 10 tasks in parallel (100 / 10 = 10)
    modules: [...]
    
  - id: heavyweight_tasks
    resources:
      cores: 20
    # Can run up to 5 tasks in parallel (100 / 20 = 5)
    modules: [...]
```

## Running with Resource Constraints

### Generate Snakefile

```bash
omnibenchmark run benchmark.yml --dry-run
```

This generates a `Snakefile` with resource directives.

### Execute with Snakemake

Control total available resources using Snakemake's `--resources` flag:

```bash
# Limit to 100 total cores
snakemake --cores 100 --resources cores=100

# Limit cores and memory
snakemake --cores 100 --resources cores=100 mem_mb=128000

# Unlimited cores (use all available)
snakemake --cores all
```

### How Snakemake Schedules Tasks

Snakemake uses a greedy scheduling algorithm:

1. Checks available resources in the pool
2. Finds ready tasks (dependencies satisfied)
3. Schedules tasks that fit within available resources
4. Repeats as tasks complete and free resources

**Example**: With `--resources cores=100`:
- Rule requiring 20 cores: max 5 in parallel
- Rule requiring 2 cores: max 50 in parallel
- Mix of both: Snakemake dynamically allocates

## Best Practices

### 1. Profile Your Tasks

Before setting resource requirements, profile representative tasks:

```bash
# Run a single task and measure
/usr/bin/time -v python script.py ...
```

Key metrics:
- CPU utilization (multi-threaded vs single-threaded?)
- Peak memory usage
- Disk I/O

### 2. Start Conservative

Begin with modest allocations and increase if needed:

```yaml
resources:
  cores: 4      # Start here
  mem_mb: 8000  # Typical for most tasks
```

### 3. Use Stage Defaults

Set sensible defaults at stage level, override only when necessary:

```yaml
stages:
  - id: analysis
    resources:
      cores: 8  # Default for all modules
    modules:
      - id: most_modules
        # Inherits 8 cores
        ...
      - id: special_case
        resources:
          cores: 16  # Needs more
        ...
```

### 4. Consider Memory vs Cores

Some tasks are:
- **CPU-bound**: Need many cores, moderate memory
- **Memory-bound**: Need less cores, lots of memory

Example:

```yaml
# CPU-intensive task
resources:
  cores: 20
  mem_mb: 8000

# Memory-intensive task
resources:
  cores: 1
  mem_mb: 64000
```

### 5. Runtime Estimates

The `runtime` field helps Snakemake prioritize long-running tasks:

```yaml
resources:
  cores: 4
  runtime: 480  # 8 hours
```

Snakemake may schedule long tasks earlier to reduce total wall time.

## Advanced: Cluster Execution

When running on HPC clusters, resource directives map to scheduler parameters:

```bash
# SLURM example
snakemake --cluster "sbatch --cpus-per-task={resources.cores} --mem={resources.mem_mb}"
```

This automatically requests appropriate resources from the cluster scheduler.

## Troubleshooting

### Tasks Not Running in Parallel

**Problem**: All tasks run sequentially despite high `--cores`.

**Solution**: Check resource directives. If rules request more cores than available, they'll queue:

```yaml
# Bad: requests 200 cores per task on 100-core system
resources:
  cores: 200  # Too high!

# Good: requests 20 cores per task
resources:
  cores: 20   # Can run 5 in parallel
```

### Out of Memory Errors

**Problem**: Tasks crash with OOM (out of memory).

**Solution**: Increase `mem_mb` allocation:

```yaml
resources:
  mem_mb: 16000  # Increase from 8000
```

On clusters, ensure memory is also requested from scheduler.

### Underutilized CPUs

**Problem**: System has 100 cores but only 10-20 are used.

**Solution**: 
1. Check if tasks are actually multi-threaded
2. Reduce per-task core allocation if tasks are single-threaded:

```yaml
# If tasks are single-threaded, don't allocate many cores
resources:
  cores: 1  # Single-threaded tasks
```

## Complete Example

See `examples/resource_allocation_example.yml` for a complete benchmark with:
- Stage-level defaults
- Module-level overrides
- Metric collector resources
- Comments explaining the strategy

## References

- [Snakemake Resources Documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources)
- [Snakemake Scheduling](https://snakemake.readthedocs.io/en/stable/executing/cli.html#useful-command-line-arguments)
- OmniBenchmark Issue #XXX: Resource allocation feature discussion
