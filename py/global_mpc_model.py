from __future__ import absolute_import
from .global_props import *


class GlobalMpcModel:
    """ Container for model that contains MPCs."""

    def __init__(self):
        self.node_to_element_map = dict()
        self.interface_elements = list()
        self.elem_directions = list()
        self.element_sets = dict()
        self.interfaces = list()
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
            for elem in node:
                interf_elems.add(elem)
        self.interface_elements = list(interf_elems)
        self.interface_elements.sort()
        return

    def establish_directions(self):
        """ Determines if an interface element is on the flange or web, and sets the axial stress component. """
        for elem in self.interface_elements:
            for bounds in self.element_sets['web']:
                if bounds[0] <= elem <= bounds[-1]:
                    self.elem_directions.append(WEB_DIRECTION)
                else:
                    self.elem_directions.append(FLANGE_DIRECTION)
        return

    def add_mpc(self, mpc):
        """ Appends an Interface to the global interface list.

        :param InterfaceProperties mpc: An interface in the model.
        """
        self.interfaces.append(mpc)
        return
