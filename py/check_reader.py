from __future__ import absolute_import
import os
from src.abaqus_inp_reader import InpReader
from src.global_mpc_model import GlobalMpcModel
from src.abaqus_analysis_writer import AbaAnalysisWriter

input_file = './testing/test_model.inp'


path_split = os.path.split(input_file)
output_file = os.path.join(path_split[0], path_split[1][:-4] + '_interf_props.txt')

model = GlobalMpcModel()
reader = InpReader(model)
reader.read(input_file)
writer = AbaAnalysisWriter(model)
writer.write(output_file)

pass
