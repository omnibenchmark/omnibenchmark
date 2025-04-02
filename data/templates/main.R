library(argparse)

parser <- ArgumentParser(description="")

# Inputs
    {INPUTS_ARGPARSE}

# Parameters
    {PARAMETERS_ARGPARSE}

# Outputs
parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved", default=getwd())
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")

args <- parser$parse_args()


    {VARIABLES_INPUTS}

    {VARIABLES_PARAMETERS}

    {VARIABLES_OUTPUTS}

