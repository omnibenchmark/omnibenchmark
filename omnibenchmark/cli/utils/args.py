def parse_extra_args(args):
    """
    args is ctx.args where ctx is a click.Context object.
    Parse extra arguments after the -- separator into a dictionary.

    POSIX convention: Arguments before -- are ob-specific flags (handled by Click).
    Arguments after -- are passed to Snakemake.

    Note: Click automatically strips the -- separator, so ctx.args only contains
    arguments that were specified after --.

    Example:
        ob run bench.yaml --cores 4 -- --dryrun --verbose

        In this case, ctx.args will be ['--dryrun', '--verbose']

    If a flag is provided without a value (e.g. --flag), it's set to True.
    """
    extra_kwargs = {}
    i = 0
    while i < len(args):
        arg = args[i]
        if arg.startswith("--"):
            key = arg[2:]
            values = []
            i += 1
            while i < len(args) and not args[i].startswith("--"):
                values.append(args[i])
                i += 1
            if not values:
                extra_kwargs[key] = True
            elif len(values) == 1:
                extra_kwargs[key] = values[0]
            else:
                extra_kwargs[key] = values
        else:
            i += 1
    return extra_kwargs
