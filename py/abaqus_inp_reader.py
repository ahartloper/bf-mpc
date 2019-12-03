class InpReader:
    """ Reads an Abaqus input file. """

    def __init__(self, model_with_mpc):
        """ Constructor.

        :param GlobalMpcModel model_with_mpc: Container to hold the Abaqus model information.
        """
        self.active_interf = None
        self.first_mpc_line = False
        self.global_model = model_with_mpc
        self.element_sets = dict()
        self.sections = {'web': [], 'flange': []}
        self.active_section = ''
        return

    def read(self, input_file):
        """ Reads the contents of the input file into the global model and returns that model.

        :param str input_file: Path to the input file.
        :return GlobalMpcModel: Model with all the interface information established.
        """
        elem_str_id = '*Element'
        mpc_str_id = '*MPC'
        node_str_id = '*Node'
        elset_str_id = '*Elset'
        section_header = '** Section'
        read_list = {'mpc': 1, 'element': 2, 'node': 3, 'elset': 4, 'section': 5, 'none': 99}

        interfaces = []
        with open(input_file, 'r') as f:
            for line in f:
                l = line.strip()
                # Check active block
                if self._check_line(l, elem_str_id):
                    active_read = read_list['element']
                elif self._check_line(l, mpc_str_id):
                    active_read = read_list['mpc']
                    self.active_interf = InterfaceProperties()
                    interfaces.append(self.active_interf)
                elif self._check_line(l, elset_str_id):
                    active_read = read_list['elset']
                    self.active_elset = self._extract_elset_name(l)
                    self.element_sets[self.active_elset] = []
                elif self._check_line(l, section_header):
                    active_read = read_list['section']
                    self.active_section = l.split(':')[-1].split('_')[0].strip()
                elif l[0] == '*' or l[0] == '':
                    active_read = read_list['none']

                # Work on active blocks
                elif active_read == read_list['element']:
                    self._read_element(l)
                elif active_read == read_list['mpc']:
                    self._read_mpc(l)
                elif active_read == read_list['elset']:
                    self._read_elset(l)
                elif active_read == read_list['section']:
                    self._read_section(l)

        # Initialize the node to element maps
        [interf.calc_interf_elem_map() for interf in interfaces]
        return self.global_model

    def _check_line(self, line, check):
        """ Checks if the start of the line is equal to check. """
        if line[:len(check)] == check:
            return True
        else:
            return False

    def _read_element(self, line):
        l = line.split(',')
        elem = int(l[0])
        connectivity = [int(li) for li in l[1:]]
        self.all_interfaces.add_element(elem, connectivity)
        return

    def _read_mpc(self, line):
        l = line.split(',')
        if self.first_mpc_line:
            self.first_mpc_line = False
            self.active_interf.global_id = int(l[1])
            [self.active_interf.add_interface_node[int(li)] for li in l[2:]]
        else:
            [self.active_interf.add_interface_node[int(li)] for li in l[1:]]
        return

    def _extract_elset_name(self, line):
        """ Returns the name of the element set from the header line. """
        l = line.split(',')
        return l[1].split('=')[-1]

    def _read_elset(self, line):
        """ Reads the data lines of an elset.

        Assumes that the type is "generate".
        """
        l = line.split(',')
        self.element_sets[self.active_elset] = [l[0], l[1]]
        return

    def _read_section(self, line):
        " Adds the section name to the flange or web group."""
        l = line.split(',')
        elset_name = l[1].split('=')[-1]
        self.sections[self.active_section].append(elset_name)
        return
