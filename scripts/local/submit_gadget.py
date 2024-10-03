#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile
import subprocess
import shutil

import re



def main():
    args = parse_args()

    params = retrieve_parameters("param.txt")
    if not args.resubmit:
        clean_dir(params['OutputDir'])

    submit_job(args)


def parse_args():
    parser = argparse.ArgumentParser(description='Submit a job to SLURM with default or specified parameters.')
    parser.add_argument('--cores', type=int, default=1, 
                        help='number of nodes')
    parser.add_argument('--script', default=None, type=int, 
                        help='The script to run. Defaults to either run.sh for normal submition or rerun.sh for a resubmition (not implemented yet)')
    parser.add_argument('--resubmit', default=None, type=int, 
                        help='if specified, then we are resubmitting a job. Options are -1 (for last restartfile) or a nonnegative integer for the snapshot to restart from.')
    args = parser.parse_args()

    return args


def clean_dir(directory):
    """
    check if a directory exists and if it does, ask the user if they want to overwrite it
    """
    if os.path.exists(directory):
        ans = input(f"Overwrite {directory}? (y/n) ")
        if ans.lower() == 'y':
            shutil.rmtree(directory)
        else:
            sys.exit(1)
    os.makedirs(directory)



def submit_job(args):
    result = subprocess.run(
        f"export SLURM_NTASKS={args.cores};" 
         "nohup bash ./run.sh > log.out 2> err.out &", 
         shell=True
    )

    print("submitted nohup job")






def _extract_path_after_substring(full_path, substring):
    # Find the index of the substring in the path
    index = full_path.find(substring)
    if index != -1:
        # Add the length of substring to index to start after the substring
        return full_path[index + len(substring):]
    else:
        # Return None or an appropriate message if the substring is not found
        return full_path


def retrieve_parameters(filename):
    """
    Retrieve parameters from a Gadget parameter file.
    The first word is the key, and the value can be seperated by arbitrary whitespace.
    Comments are denoted by '#' or '%' anywhere in the line.

    Parameters
    ----------
    filename : str

    Returns
    -------
    parameters : dict
        A dictionary of parameters
    """

    parameters = {}

    with open(filename, 'r') as file:
        for line in file:

            key, value = _read_param_line(line)

            if key:
                parameters[key] = value

    return parameters



def _read_param_line(line):
    """
    Read a line from a Gadget parameter file.

    Parameters
    ----------
    line : str
        The line to read

    Returns
    -------
    key, value: str, str
        The name and value of the parameter. Returns None, None if the line is empty.
    """

    key = None
    value = None

    line_no_comments = line.split('#')[0].split('%')[0].strip()

    if line_no_comments:
        parts = line_no_comments.split()
        if len(parts) == 2:
            key = parts[0].strip()
            value = parts[1].strip()
        else:
            print(f"Warning: Ignoring line '{line_no_comments}'")

    return key, value







if __name__ == '__main__':
    main()
