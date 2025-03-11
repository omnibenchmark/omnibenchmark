import argparse
import sys
import os


def main():
    parser = argparse.ArgumentParser(description="")

    # Inputs
    INPUTS_ARGPARSE

    # Parameters
    PARAMETERS_ARGPARSE

    # Outputs
    parser.add_argument(
        "--output_dir", type=str, help="output directory to store data files."
    )
    parser.add_argument("--name", type=str, help="name of the dataset", default="")

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    VARIABLES_INPUTS

    VARIABLES_PARAMETERS

    VARIABLES_OUTPUTS


if __name__ == "__main__":
    main()
