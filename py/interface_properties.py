class InterfaceProperties:
    def __init__(self):
        self.global_id = None
        self.interface_nodes = list()
        self.n2e_map = dict()
        return

    def add_interface_node(self, node_id):
        self.interface_nodes.append(node_id)
        return

    def calc_interf_elem_map(self, all_interf):
        for n in self.interface_nodes:
            elems = all_interf[n]
            self.n2e_map[n] = []
            [n2e_map[n].append(e) for e in elems if e in all_interf.interface_elements]
        return
