"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

import os
from timeit import default_timer as timer
from concurrent.futures import ProcessPoolExecutor, as_completed
from psutil import Process

import parser
import argparser
import contacts

def main():
    """
    Main function for the script.

    This function parses command-line arguments, sets up the environment based on the specified mode,
    and runs the appropriate processing function (single-core or multi-core) on the input files.
    It also manages core affinity, output folder creation, and timing of the entire process.
    """
    global_time_start = timer()
    
    file_list, core, output = argparser.cl_parse()
    
    if core is not None:  # Set specific core affinity
        Process(os.getpid()).cpu_affinity(core)
        print("Multicore mode selected.")
        if len(core) == 1: # One specific core
            print(f"Running on core {core[0]}.")
        elif core[-1] - core[0] == len(core) - 1:  # Check if it's a range
            print(f"Running on cores {core[0]} to {core[-1]}\nTotal number of cores: {len(core)}.")
        else: # List
            print(f"Running on cores: {', '.join(map(str, core))}\nTotal number of cores: {len(core)}.")
    else:
        print("Running on single mode with no specific core.") 
               
    if output:
        output_folder = "./outputs/"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
    else:
        output_folder = None
        
    process_func = single if core is None else multi
    process_func(file_list, output_folder, core)

    print(f"Total time elapsed: {(timer() - global_time_start):.3f}s\n")


def single(file_list, output_folder, core=None):
    """
    Processes a list of files in single-core mode.

    Args:
        file_list (list): List of file paths to process.
        output_folder (str): The directory where output files will be saved (or None if no output).

    This function processes each file in the list sequentially, detects contacts, and outputs
    the results to the console or to a file, depending on the 'output' flag.
    """
    for file in file_list:
        try:
            result = process_file(file)
            process_result(result, output_folder)
        except Exception as e:
            print(f"Error: {e}")


def multi(file_list, output_folder,core):
    """
    Processes a list of files in multi-core mode using parallel processing.

    Args:
        file_list (list): List of file paths to process.
        output_folder (str): The directory where output files will be saved (or None if no output).

    This function processes the files in the list using a process pool with the specified number of cores.
    """

    with ProcessPoolExecutor(max_workers=len(core)) as executor:

        futures = {executor.submit(process_file, file): file for file in file_list}
        for future in as_completed(futures):
            try:
                process_result(future.result(), output_folder)
            except Exception as e:
                print(f"Error: {e}")
            finally:
                del futures[future]  # Clean memory

def process_file(file_path):
    """
    Processes a single file for contact detection.

    Args:
        file_path (str): Path to the file to be processed.

    Returns:
        tuple: A tuple containing the processed Protein object, the list of detected contacts, and the processing time.
        None: If the file cannot be processed or an error occurs.

    This function parses the PDB or mmCIF file, detects contacts, and returns the results. 
    If an error occurs during processing, it logs the error and returns None.
    """
    start_time = timer()

    try:
        parsed_data = parser.parse_pdb(file_path) if file_path.endswith(".pdb") else parser.parse_cif(file_path)
            
        if parsed_data.true_count() > 25000:  # Skip very large proteins (customizable)
            print(f"Skipping ID '{parsed_data.id}'. Size: {parsed_data.true_count()} residues")
            return None

        contacts_list = contacts.contact_detection(parsed_data)
        process_time = timer() - start_time
        return parsed_data, contacts_list, process_time

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None


def process_result(result, output_folder):
    """
    Handles the result of processing a file.

    Args:
        result (tuple): A tuple containing the processed Protein object, contacts list, and processing time.
        output_folder (str): The directory where output files will be saved.
    """
    if result:
        protein, contacts_list, process_time = result
        output_data = f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.3f}s"

        print(output_data)
        if output_folder:
            with open(f"{output_folder}/{protein.id}_contacts.txt", "w") as f:
                f.write(output_data + "\n")
                f.write(contacts.show_contacts(contacts_list))


if __name__ == "__main__":
    main()