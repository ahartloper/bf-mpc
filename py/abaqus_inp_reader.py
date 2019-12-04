from interface_properties import InterfaceProperties


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
        self.active_elset = ''
        return

    def read(self, input_file):
        """ Reads the contents of the input file into the global model and runs the model processing.

        :param str input_file: Path to the input file.

        Notes:
            - Elements are defined following: *Element, type=<elem_type>
            - Each section definition is preceded by: ** Section: <flange/web>_<arbitrary text without underscores>
            - All the element sets that contain either flange/web elements are defined using the "generate" option
                with a step size of 1 e.g.,
                    *Elset, elset=<set_name>, generate
                    <start_elem>, <end_elem>, 1
        """
        # todo: make the reader more elegant by check the keyword if there is a "*" then taking the first part of the
        #  line before the comma as the input (so don't have to use *Element, type=S4R)
        element_keyword_id = '*Element,'
        mpc_keyword_id = '*MPC'
        elset_keyword_id = '*Elset'
        generate_id = 'generate'
        section_header_id = '** Section'
        section_keyword_id = '*Shell Section'
        read_list = {'mpc': 1, 'element': 2, 'elset': 3, 'section': 4, 'none': 99}
        accepted_element_types = ['S4R']

        interfaces = []
        with open(input_file, 'r') as f:
            for line in f:
                l = line.strip()
                # Check active block
                if self._check_line(l, element_keyword_id):
                    l2 = l.split(',')
                    if l2[1].split('=')[-1] in accepted_element_types:
                        active_read = read_list['element']
                elif self._check_line(l, mpc_keyword_id):
                    active_read = read_list['mpc']
                    self.active_interf = self.global_model.add_mpc(InterfaceProperties())
                    self.first_mpc_line = True
                elif self._check_line(l, elset_keyword_id):
                    l2 = l.split(',')
                    if l2[-1].strip() == generate_id:
                        active_read = read_list['elset']
                        self.active_elset = self._extract_elset_name(l)
                elif self._check_line(l, section_header_id):
                    active_read = read_list['section']
                    self.active_section = self._extract_section_name(l)
                elif self._check_line(l, section_keyword_id):
                    self._read_section(l)
                    active_read = read_list['none']
                elif l[0] == '*' or l[0] == '':
                    active_read = read_list['none']

                # Work on active blocks
                elif active_read == read_list['element']:
                    self._read_element(l)
                elif active_read == read_list['mpc']:
                    self._read_mpc(l)
                elif active_read == read_list['elset']:
                    self._read_elset(l)
                # elif active_read == read_list['section']:

        # Finalize the global model
        self.global_model.element_sets = self.element_sets
        self.global_model.sections = self.sections
        self.global_model.process_model()
        return

    def _check_line(self, line, check):
        """ Returns the result of testing if the start of line is equal to check. """
        if line[:len(check)] == check:
            return True
        else:
            return False

    def _read_element(self, line):
        """ Adds an element to the global model. """
        l = line.split(',')
        elem = int(l[0].strip())
        connectivity = [int(li.strip()) for li in l[1:]]
        self.global_model.add_element(elem, connectivity)
        return

    def _read_mpc(self, line):
        """ Adds the properties of an interface to the global model. """
        l = line.split(',')
        if self.first_mpc_line:
            # If first line we have the global ID
            self.first_mpc_line = False
            self.active_interf.global_id = int(l[1].strip())
        else:
            # If later lines we just have a shell node
            self.active_interf.add_interface_node(int(l[1].strip()))
        for n in l[2:]:
            n2 = n.strip()
            if not n2 == '':
                self.active_interf.add_interface_node(int(n2))
        return

    def _extract_elset_name(self, line):
        """ Returns the name of the element set from the header line. """
        l = line.split(',')
        return l[1].split('=')[-1]

    def _extract_section_name(self, line):
        """ Returns the name of the section.

        Notes:
            - Assumes that the header is as follows: ** Section: <flange/web>_<arbitrary text without underscores>
        """
        return line.split(':')[-1].split('_')[0].strip()

    def _read_elset(self, line):
        """ Reads the data lines of an elset.

        Notes:
            - Assumes that the type is "generate" and the interval is 1.
        """
        l = line.split(',')
        self.element_sets[self.active_elset] = [int(l[0]), int(l[1])]
        return

    def _read_section(self, line):
        """ Adds the section name to the flange or web group."""
        l = line.split(',')
        elset_name = l[1].split('=')[-1]
        self.sections[self.active_section].append(elset_name)
        return
