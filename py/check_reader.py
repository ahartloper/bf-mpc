from __future__ import absolute_import
from abaqus_inp_reader import InpReader
from global_mpc_model import GlobalMpcModel

input_file = './testing/test_model.inp'

model = GlobalMpcModel()
reader = InpReader(model)
reader.read(input_file)

pass
