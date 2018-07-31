import sys
import os
# Edit this path to the location of the PDBModule folder on your computer
sys.path.append("../PDBModule")

from Pseudoatom import Pseudoatom

from pathlib import Path



class Writer:
    def __init__(self, path):
        self.path = Path(path)
        self.com_path = Path(self.path.parent,self.path.stem+".com")

    def write_particle(self, particle):

        try:
            os.remove(self.path)
            os.remove(self.com_path)
        except:
            pass
        file = open(self.path, "w")

        def idf(object):
            return particle.intid_dict[id(object)]

        particle.intid_update()

        p = Pseudoatom(0,particle.centroid)
        p.atom_name = "Mg"
        file.write(p.generate_entry())

        for i, vertex in enumerate(particle.verteces):
            coord = vertex.coord
            intid = idf(vertex)
            p = Pseudoatom(intid, vertex.coord)
            p.atom_name = vertex.atom_name
            file.write(p.generate_entry())




        for i, vertex in enumerate(particle.debug_verteces):
            coord = vertex.coord
            intid = idf(vertex)
            p = Pseudoatom(intid, vertex.coord)
            if p.atom_name == "C":
                p.atom_name = "P"
            else:
                p.atom_name = vertex.atom_name
            file.write(p.generate_entry())

        # for i, edge in enumerate(particle.edges):
        #     intid = idf(edge)
        #     p = Pseudoatom(intid, edge.coord)
        #     p.atom_name = edge.atom_name
        #     file.write(p.generate_entry())

        for edge in particle.edges:
            v1_intid = idf(edge.vertex1)
            v2_intid = idf(edge.vertex2)
            connect_entry = "CONECT" + str(v1_intid).rjust(5) + str(v2_intid).rjust(5) + "\n"
            file.write(connect_entry)
        file.close()

        # chimera com file to get good representation
        file = open(self.com_path, "w")
        file.write("open "+self.path.name+"\n")
        file.write("color #6495ed #\n")
        file.write("represent bs #\n")
        file.write("setattr m stickScale 10\n")
        file.write("setattr m ballScale 5\n")
        file.close()