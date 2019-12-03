from __future__ import absolute_import
from global_props import *


class GlobalMpcModel:
    """ Container for model that contains MPCs."""

    def __init__(self):
        self.node_to_element_map = dict()
        self.interface_elements = list()
        self.elem_directions = list()
        self.element_sets = dict()
        self.interfaces = list()
        self.all_interface_nodes = list()
        self.sections = dict()
        return

    def add_element(self, elem_id, connectivity):
        """ Adds and element to the node-to-element map.
        :param int elem_id: Element tag.
        :param list connectivity: (int) Nodes connected to the element.
        """
        for ci in connectivity:
            if ci not in self.node_to_element_map:
                self.node_to_element_map[ci] = [elem_id]
            else:
                self.node_to_element_map[ci].append(elem_id)
        return

    def establish_interface_elements(self):
        """ Collects all the interface elements and creates a sorted list of the elements. """
        interf_elems = set()
        for node in self.node_to_element_map:
            for elem in self.node_to_element_map[node]:
                if node in self.all_interface_nodes:
                    interf_elems.add(elem)
        self.interface_elements = list(interf_elems)
        self.interface_elements.sort()
        return

    def establish_directions(self):
        """ Determines if an interface element is on the flange or web, and sets the axial stress component.
        """
        for elem in self.interface_elements:
            for elset in self.element_sets:
                bounds = self.element_sets[elset]
                if bounds[0] <= elem <= bounds[-1]:
                    if elset in self.sections['web']:
                        self.elem_directions.append(WEB_DIRECTION)
                    elif elset in self.sections['flange']:
                        self.elem_directions.append(FLANGE_DIRECTION)
                    else:
                        raise ValueError('Section not found for element {0:d}'.format(elem))
        return

    def add_mpc(self, mpc):
        """ Appends an Interface to the global interface list.

        :param InterfaceProperties mpc: An interface in the model.
        :returns InterfaceProperties: The pointer to the interface passed in.
        """
        self.interfaces.append(mpc)
        return mpc

    def set_interface_node_to_element_index_maps(self):
        """ Sets the node-to-element maps for all the interfaces in the model. """
        for interf in self.interfaces:
            interf.set_node_to_elem_index_map(self.node_to_element_map, self.interface_elements)
        return

    def establish_all_interface_nodes(self):
        """ Collects all the interface nodes from the interfaces. """
        for interf in self.interfaces:
            for n in interf.interface_nodes:
                self.all_interface_nodes.append(n)
        return

    def process_model(self):
        """ Process the model to prepare for output. """
        self.establish_all_interface_nodes()
        self.establish_interface_elements()
        self.establish_directions()
        self.set_interface_node_to_element_index_maps()
