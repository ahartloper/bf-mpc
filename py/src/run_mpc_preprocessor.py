from __future__ import absolute_import, print_function
import os
import sys
from .abaqus_inp_reader import InpReader
from .global_mpc_model import GlobalMpcModel
from .abaqus_analysis_writer import AbaAnalysisWriter
from .sections_reader import section_reader


def run_mpc_preprocessor(input_file, section_file):
    """ Runs the pre-processor for Abaqus input files that contain inelastic best-fit MPCs.
    : param str input_file: Path to the input (.inp) file to be processed.

    Notes:
        - Section definitions need to be preceded by the comment: ** Section: <flange/web>_<section name>
    """
    path_split = os.path.split(input_file)
    # output_file = os.path.join(path_split[0], path_split[1][:-4] + '_interf_props.txt')
    output_file = os.path.join(path_split[0], 'interface_props.txt')

    model = GlobalMpcModel()
    reader = InpReader(model)
    reader.read(input_file)
    model.section_props = section_reader(section_file)
    writer = AbaAnalysisWriter(model)
    writer.write(output_file)
    return


if __name__ == '__main__':
    if len(sys.argv) == 2:
        run_mpc_preprocessor(sys.argv[1])
    else:
        print('Usage:\n\tpython run_mpc_preprocessor <path_to_input_file.inp>')
