import click
import sys
import time

from module03.src.io.read_in import read_in
from module03.src.algorithm.signal_analysis import analyse_homo_signals
from module03.src.algorithm.create_histograms import create_homo_signal_histograms
from module03.src.io.write_out import write_output

from utils.utils import *


@click.command()
@click.option("-i", "--input-filepath", default="data/out/unique_protein_list/sample_data.sqcs")
@click.option("-o", "--output-directory", default="data/out/homo_signal/")
@click.option("-v", "--verbose-level", default=2)
def main(input_filepath, output_directory, verbose_level):
    verbose_print("Start Homo-signal analysis", 0, verbose_level)
    start_time = time.time()

    # Translate eventual windows paths and evaluate value of boolean inputs
    input_filepath, output_directory = translate_windowsos_path([input_filepath, output_directory])

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, output_directory):
        # Check output directory
        output_directory = create_out_path(output_directory, input_filepath)

        # Read dataset and add columns for results
        verbose_print("Read input", 0, verbose_level)
        data = read_in(input_filepath)

        # Analyse homo signals
        verbose_print("Analyse homo signals", 0, verbose_level)
        data = analyse_homo_signals(data)

        # Clean dataset for output
        data = clean_dataset(data)

        # Create Homo-signal statistic histograms
        verbose_print("Create homo-signal histograms", 0, verbose_level)
        create_homo_signal_histograms(data, input_filepath.split('/')[-1], output_directory)

        # Write Output
        verbose_print("Write output", 0, verbose_level)
        write_output(data, input_filepath.split('/')[-1], output_directory)

    verbose_print(f"\nEnd script (Elapsed time: {round_self(time.time() - start_time, 2)}s)", 0, verbose_level)
    verbose_print("===================================", 0, verbose_level)
    sys.exit()


def inputs_valid(input_filepath, output_directory):
    # check validity of inputted parameters
    #
    # input input_filepath: str, output_directory: str
    # return inputs_valid: bool

    # check whether an inputfile with the extension .sqcs is specified
    if input_filepath and input_filepath.endswith(".sqcs"):
        return True
    else:
        print(f"Error! The parameter \"input-filepath\" was not given or was not ending with the '.sqcs' extension "
              f"(given: {input_filepath}).")
    return False