#!/bin/bash

while [[ $# -gt 0 ]]; do
  case $1 in
{INPUTS_ARGPARSE}
{PARAMETERS_ARGPARSE}
    --output_dir)
      output_dir="$2"
      shift
      shift
      ;;
    --name)
      name="$2"
      shift
      shift
      ;;
    --*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      shift
      ;;
  esac
done

{VARIABLES_INPUTS}
{VARIABLES_PARAMETERS}
{VARIABLES_OUTPUTS}