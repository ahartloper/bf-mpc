class InterfaceProperties:
    """ Container for the global ID and node-to-element map for a single MPC interface. """

    def __init__(self):
        self.global_id = None
        self.interface_nodes = list()
        self.node_to_elem_ind_map = dict()
        return

    def add_interface_node(self, node_id):
        """ Adds a shell node to the set of shell nodes for this interface.
        :param int node_id: Tag of the node.
        """
        self.interface_nodes.append(node_id)
        return

    def set_node_to_elem_index_map(self, global_node_to_elem_map, interface_elements):
        """ Sets the shell elements that each shell node is connected to for this interface.
        :param dict global_node_to_elem_map: The node-to-element map for all the nodes/elements in the model.
        :param list interface_elements: All the interface elements in the model.
        """
        for n in self.interface_nodes:
            elements = global_node_to_elem_map[n]
            self.node_to_elem_ind_map[n] = []
            for e in elements:
                ind = interface_elements.index(e)
                self.node_to_elem_ind_map[n].append(ind)
        return
