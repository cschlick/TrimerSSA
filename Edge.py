from lib import assertlength



class Edge:
    def __init__(self, verteces):  # initializes on a list of verteces
        self.verteces = verteces  # a set of only two
        self.vertex1 = verteces[0]
        self.vertex2 = verteces[1]
        self.coord = self.__generate_coord()
        self.full = False

        self.trimers = set()
        self.deleted = False
        self.__atom_name = None

    @property
    def atom_name(self):
        if self.__atom_name != None:
            return self.__atom_name
        elif self.full == False:
            return "O"
        else:
            return "N"
    @atom_name.setter
    def atom_name(self,name):
        self.__atom_name = name

    @property
    def verteces(self):
        return self.__verteces

    @verteces.setter
    def verteces(self, verteces):
        verteces = set(verteces)
        assertlength(verteces,2)
        for v in verteces:
            v.edges.add(self)
        self.__verteces = verteces

    @property
    def secondary_edges(self):
        secondary_edges = set()
        for vertex in self.verteces:
            for edge in vertex.edges:
                if edge is not self:
                    secondary_edges.add(edge)
        return secondary_edges

    def delete(self):
        self.deleted = True

        for vertex in self.verteces:
            if len(vertex.edges) == 1:
                vertex.delete()
            try:
                vertex.edges.remove(self)
            except:
                # print(vertex.edges,self)
                pass

    def __generate_coord(self):
        vertex_list = list(self.verteces)
        v1, v2 = vertex_list[0], vertex_list[1]
        midpoint = (v1.coord + v2.coord) / 2
        return midpoint


