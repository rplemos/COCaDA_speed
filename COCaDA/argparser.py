"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

from sys import exit
from argparse import ArgumentParser, ArgumentError, ArgumentTypeError


def cl_parse():
    """
    Parses command-line arguments for a PDB/mmCIF parser and contact detection tool.

    Returns:
        tuple: A tuple containing the parsed values:
            - files (list): List of input files.
            - ncores (int): Number of cores to use.
            - mode (str): Selected processing mode.
            - selcore (str): Selected core (if provided).
            - output (bool): Whether to output results to files.

    Raises:
        ArgumentError: If there's an issue with the command-line arguments.
        ValueError: If an invalid processing mode is specified.
        Exception: For any other unexpected errors during argument parsing.
    """
    
    try:
        parser = ArgumentParser(description='PDB/mmcif parser and fast contact detection using flexible CA distances.')
        parser.add_argument('-f', '--files', nargs='+', required=True, type=validate_file, help='List of files in pdb/cif format (at least one required).')
        parser.add_argument('-m', '--mode', required=False, default='Single', help='Select "SingleCore" or "MultiCore" mode.')
        parser.add_argument('-c', '--ncores', type=int, required=False, default=0, help='Number of cores to use (only needed on Multi mode). Default runs with all available cores')
        parser.add_argument('-s', '--selcore', required=False, help='Select specific core.')
        parser.add_argument('-o', '--output', required=False, action='store_true', help='Outputs the results to files in ./outputs.')

        args = parser.parse_args()

        files = args.files
        ncores = args.ncores
        mode = args.mode
        modes = ["Single", "Multi"]
        if mode not in modes:
            raise ValueError("Invalid Mode!")
        selcore = args.selcore
        output = args.output
        
    except ArgumentError as e:
        print(f"Argument Error: {str(e)}")
        exit(1)

    except ValueError as e:
        print(f"Error: {str(e)}")
        exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        exit(1)
    
    return files, ncores, mode, selcore, output
        
        
def validate_file(value):
    """
    Validates a file path to ensure it has a proper extension for PDB or mmCIF files.

    If the file has a valid extension, the function returns the file path. Otherwise, it raises an `ArgumentTypeError`.

    Args:
        value (str): The file path to validate.

    Returns:
        str: The validated file path.

    Raises:
        ArgumentTypeError: If the file does not have a valid extension.
    """
    
    if value.endswith('.pdb') or value.endswith('.cif'):
        return value
    else:
        raise ArgumentTypeError(f"{value} is not a valid file. File must end with '.pdb' or '.cif'")
