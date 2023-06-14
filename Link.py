class Link:
    """
    """
    
    def __init__(self, neighbour1=-1, neighbour2=-1, bond=False):
        self.neighbour1 = neighbour1
        self.neighbour2 = neighbour2
        self.bond = bond


    def set_neighbour(self, new_neighbour):
        self.neighbour = new_neighbour
    