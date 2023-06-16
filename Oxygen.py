class Oxygen:
    """
    Author: Gabe Miles
    Date: 6/14/23
    Description: This class represents an oxygen atom within an ice lattice.
    it keeps track of its position within the lattice, along with the bonds
    or 'links' between neighboring oxygens and the number of bonds to protons
    it has. Each oxygen should have no more than 4 neighbors accoring to normal
    ice rules. 
    """

    def __init__(self):
        """
        Attributes
        - coord: this is an array of length 3 that contains the 3D coordinate 
          of the oxygen in space.
        - links: this is an array of Link objects that keep track of
          the neighboring oxygens. This should not exceed 4 objects.
        - nbonds: this is an integer value of the number of protons associated 
          to this oxygen.
        """
        self.coord = []
        self.nneighbours = 0
        self.links_index = []
        self.nbonds = 0

    def set_coord(self, new_coord):
        """
        This function sets a new 3D coordinate for this oxygen atom. It checks
        that the new coordinate is an array of length 3.

        Input
        - new_coord: an array of length 3 that is the xyz coordinate of the oxygen atom

        """
        if len(new_coord) != 3:
            pass
        else:
            self.coord = new_coord
    
    def add_link(self, index):
        """
        """
        self.add_link_index(index)
        self.nneighbours += 1

    def add_bond(self):
        self.nbonds += 1

    def remove_bond(self):
        self.nbonds -= 1
    
    def add_link_index(self, new_link_index):
        self.links_index.append(new_link_index)

    def set_nbonds(self, new_nbonds):
        self.nbonds = new_nbonds

    def get_neighbour(self, index):
        return self.ne