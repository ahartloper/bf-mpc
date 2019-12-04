MODIFICATION_FOR_1_INDEX = 1


class AbaAnalysisWriter:
    def __init__(self, global_model):
        """ Constructor.
        :param GlobalMpcModel global_model: Model that will be output to file.
        """
        self.model = global_model
        self.ELEMENTS_HEADER = '*elements'
        self.INTERFACES_HEADER = '*interface'
        pass

    def write(self, output_file):
        """ Writes the file to be read by the Abaqus analysis.
        :param str output_file: Path where the file is to be written.

        Notes:
            - The indices for the elements are output in 1-indexed format (node-to-element maps)
        """
        with open(output_file, 'w') as f:
            l = ','.join([self.ELEMENTS_HEADER, str(len(self.model.interface_elements))]) + '\n'
            f.write(l)
            for i, elem in enumerate(self.model.interface_elements):
                l = ','.join([str(elem), str(self.model.elem_directions[i])]) + '\n'
                f.write(l)
            for interf in self.model.interfaces:
                l = ','.join([self.INTERFACES_HEADER, str(len(interf.interface_nodes))]) + '\n'
                f.write(l)
                f.write(str(interf.global_id) + '\n')
                for n in interf.interface_nodes:
                    n2e_map = interf.node_to_elem_ind_map[n]
                    n2e_map = [ei + MODIFICATION_FOR_1_INDEX for ei in n2e_map]
                    if len(n2e_map) == 1:
                        l = '{0:d},0,0\n'.format(n2e_map[0])
                    elif len(n2e_map) == 2:
                        l = '{0:d},{1:d},0\n'.format(n2e_map[0], n2e_map[1])
                    else:
                        l = '{0:d},{1:d},{2:d}\n'.format(n2e_map[0], n2e_map[1], n2e_map[2])
                    f.write(l)
        return
