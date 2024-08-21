"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

import parser
import argparser
import contacts

import os
from timeit import default_timer as timer
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from psutil import Process


def main():
    """
    Main function for the script.

    This function parses command-line arguments, sets up the environment based on the specified mode,
    and runs the appropriate processing function (single-core or multi-core) on the input files.
    It also manages core affinity, output folder creation, and timing of the entire process.
    """

    global_time_start = timer()
    file_list, core, mode, cnum, output = argparser.cl_parse()
    
    if cnum: # specific core to run
        p = Process(os.getpid())
        p.cpu_affinity([int(cnum)])
        print(f"Running on core {cnum}")
        
    if output:
        output_folder = "./outputs/"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
    else:
        output_folder = None

        
    if mode == "Single":
        single(file_list, output, output_folder)
    elif mode == "Multi":
        multi(file_list, output, output_folder, core)

    print(f"Total time elapsed: {timer() - global_time_start}\n")


def single(file_list, output, output_folder):
    """
    Processes a list of files in single-core mode.

    Args:
        file_list (list): List of file paths to process.
        output (bool): Whether to save the results to a file.
        output_folder (str): The directory where output files will be saved.

    This function processes each file in the list sequentially, detects contacts, and outputs
    the results to the console or to a file, depending on the 'output' flag.
    """
    
    for file in file_list:
        try:
            result = process_file(file)
            if result:
                protein, contacts_list, process_time = result
                if output:
                    with open(f"{output_folder}/{protein.id}_contacts.txt", "w") as f:
                        print(f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.4f}")
                        f.write(f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.4f}\n")
                        f.write(contacts.show_contacts(contacts_list))
                else:
                    print(f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.4f}")
                    
        except Exception as e:
            print(f"Error: {e}")
            
            
def multi(file_list, output, output_folder, core):
    """
    Processes a list of files in multi-core mode using parallel processing.

    Args:
        file_list (list): List of file paths to process.
        output (bool): Whether to save the results to a file.
        output_folder (str): The directory where output files will be saved.
        core (int): The number of cores to use for parallel processing. If set to 0, all available cores are used.

    This function processes the files in the list using a process pool with the specified number of cores.
    It handles the results in parallel and outputs them to the console or a file, depending on the 'output' flag.
    """
    
    core = cpu_count() if core == 0 else core
    print(f"Starting processing with {core} cores") 
    
    with ProcessPoolExecutor(max_workers=core) as executor:
        futures = {executor.submit(process_file, file): file for file in file_list}
        
        for future in as_completed(futures):
            try:
                result = future.result()
                if result:
                    protein, contacts_list, process_time = result
                    if output:
                        with open(f"{output_folder}/{protein.id}_contacts.txt", "w") as f:
                            print(f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.4f}")
                            f.write(f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.4f}\n")
                            f.write(contacts.show_contacts(contacts_list))
                    else:
                        print(f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.4f}")
                        
            except Exception as e:
                print(f"Error: {e}")
                
            finally:
                del futures[future] # cleans memory to avoid bloating
           
                
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
            
        if parsed_data.true_count() > 25000: # skips very big proteins (customizable)
            print(f"Skipping ID '{parsed_data.id}'. Size: {parsed_data.true_count()} residues")
            return None

        contacts_list = contacts.contact_detection(parsed_data)
        
        process_time = timer() - start_time
        return parsed_data, contacts_list, process_time
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None


if __name__ == "__main__":
    main()