class Link:
    """
    """
    
    def __init__(self, oxy1=-1, oxy2=-1):
        self.oxy1 = oxy1
        self.oxy2 = oxy2
        self.oxy1_bond = False
        self.oxy2_bond = False


    def set_oxy1(self, new_oxy):
        self.oxy1 = new_oxy

    def set_oxy2(self, new_oxy):
        self.oxy2 = new_oxy

    def make_oxy1_bond(self):
        self.oxy1_bond = True
        self.oxy2_bond = False
    
    def make_oxy2_bond(self):
        self.oxy2_bond = True
        self.oxy1_bond = False

    def swap_bond(self):
        self.oxy1_bond = not self.oxy1_bond
        self.oxy2_bond = not self.oxy2_bond