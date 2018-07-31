class Vertex:
    def __init__(self, coord):
        self.coord = coord
        self.edges = set()  # references to edges which attatch. Added upon edge creation
        self.deleted = False
        self.atom_name = "C"

    @property
    def secondary_verteces(self):  # returns set of secondary (one edge away verteces)
        secondary_verteces = set()
        for edge in self.edges:
            for vertex in edge.verteces:
                secondary_verteces.add(vertex)
        try:
            secondary_verteces.remove(self)
        except:
            pass
        return secondary_verteces
    @property
    def trimers(self):
        trimers = set()
        for edge in self.edges:
            for trimer in edge.trimers:
                trimers.add(trimer)
        return trimers

    def delete(self):
        self.deleted = True