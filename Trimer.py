
import sys
sys.path.append("transformations.py")
from transformations import rotation_matrix, angle_between_vectors, unit_vector

sys.path.append("PDBModule")
from Pseudoatom import Pseudoatom




import random
import numpy as np

from Edge import Edge
from Top import Top
from Vertex import Vertex
from lib import assertlength

class Trimer:
    def __init__(self,edges): # a list of edges
        """
        Each trimer has three edges, and three verteces
        """

        self.edges = edges
        vertex_set = set()
        for edge in self.edges:
            for v in edge.verteces:
                vertex_set.add(v)
        self.verteces = vertex_set
        self.coord = self.__generate_coord()
        self.tops = None
        self.deleted = False


    @property
    def full(self):
        full_count = len([edge for edge in self.edges if edge.full ==True])
        if full_count <3:
            return False
        elif full_count == 3:
            return True
        else:
            print("Error: full count ", full_count)

    @property
    def adjacent_trimers(self):
        adjacent_trimers = set()
        for edge in self.edges:
            for trimer in edge.trimers:
                adjacent_trimers.add(trimer)
        return adjacent_trimers

    @property
    def open_edges(self):  # returns list
        return [edge for edge in self.edges if edge.full == False]


    @property
    def edges(self):
        return self.__edges

    @edges.setter
    def edges(self, edges):
        edges = set(edges)
        assertlength(edges,3)
        for edge in edges:
            edge.trimers.add(self)
        self.__edges = edges

    @property
    def verteces(self):
        return self.__verteces

    @verteces.setter
    def verteces(self, verteces):
        verteces = set(verteces)
        assert (len(verteces) == 3)
        self.__verteces = verteces

    @property
    def secondary_trimers(self):
        trimers = set()
        for edge in self.edges:
            for trimer in edge.trimers:
                trimers.add(trimer)
        trimers.remove(self)
        return trimers

    @property
    def tops(self):
        return self.__tops
    @tops.setter
    def tops(self,tops):
        if tops is None:
            v1,v2,v3 = list(self.verteces)
            n1 = v2.coord-v1.coord
            n2 = v3.coord - v1.coord
            vec = np.cross(n1, n2)
            d = 25
            top1 = Top(self.coord + d * unit_vector(vec))
            top2 = Top(self.coord - d * unit_vector(vec))
            top1.trimer = self
            top2.trimer = self
            tops = {top1,top2}
        else:
            assert(len(tops)==2)
        self.__tops = tops




    def delete(self):
        self.deleted = True
        for edge in self.edges:
            edge.trimers.remove(self)
            if edge.full == False:

                edge.delete()
            else:
                edge.full = False




    def __generate_coord(self):
        vertex_list = list(self.verteces)
        v1, v2, v3 = vertex_list[0], vertex_list[1], vertex_list[2]
        midpoint = (v1.coord + v2.coord + v3.coord) / 3
        return midpoint


    def opposite_vertex(self, edge):  # returns the vertex opposite a given edge
        return next(iter(self.verteces - edge.verteces))


    def get_new_point(self,reference_point,start_point,angle,axis):  # where angle is radians of the stem length (edge to opposite vertex)

        M = rotation_matrix(angle, axis)
        new_point = np.dot(start_point, M[:3, :3].T)
        return new_point+reference_point

    def add(self, particle, template,adding_edge,debug=False):  # returns a candidate new vertex coord
        #     V4,V5,V6,V7
        #
        #          V1-------V3
        #         /  \     /
        #       /     \  /
        #      V5------V2
        #
        # E1 is midpoint of V1,V2 (the adding edge)
        # V3 is the vertex opposite the adding edge
        # V5 is the new vertex we need to calculate
        #
        #V4 will be an intermediate point along the axis E1-V3 which is the appropriate stem length for the template
        #    (essentially making the E1V3 distance longer or shorter)
        #
        # The convention I used is upper case V is the Vertex object. Lower case v is the actual coordinate (a numpy vector)



        assert (self.full == False)

        E1 = adding_edge
        (V1, V2) = adding_edge.verteces
        V3 = self.opposite_vertex(adding_edge)

        # # maybe extend length a bit along the V1V3 vector to get the appropriate stem length for the template

        E1V3 = E1.coord-V3.coord
        V1V3_dist = np.linalg.norm(V3.coord - E1.coord)
        d = template.stem_length/V1V3_dist

        v4_1 = E1.coord - d * E1V3
        v4_2 = E1.coord + d * E1V3
        dist1 = np.linalg.norm(v4_1-V3.coord)
        dist2 = np.linalg.norm(v4_2-V3.coord)

        if dist1<dist2: # pick the value closer to V3 to become V4
            v4 = v4_1
        else:
            v4 = v4_2

        # set origin of rotation to E1
        e1 = E1.coord
        v4 = v4-e1


        # get two new points with consistent angle for first rotation
        axis = V1.coord-V2.coord
        axis = unit_vector(axis)
        a1 = template.angle_radians
        a2 = -1*(template.angle_radians)
        v5_1 =self.get_new_point(e1,v4,a1,axis)
        v5_2 =self.get_new_point(e1,v4,a2, axis)

        if debug is True: # add debug verteces for the chosen verteces
            dbgv1 = Vertex(v5_1)
            particle.debug_verteces.add(dbgv1)
            dbgv2 = Vertex(v5_2)
            particle.debug_verteces.add(dbgv2)
            dbgv3 = Vertex(E1.coord)
            particle.debug_verteces.add(dbgv3)


        # pick closest to centroid
        if np.linalg.norm(v5_1 - particle.centroid) < np.linalg.norm(v5_2 - particle.centroid):
            new_point = v5_1

        else:
            new_point = v5_2

        new_vertex = Vertex(new_point) # the final V5

        return new_vertex



    def trilaterate3D(self, p1, r1, p2, r2, p3, r3, centroid): # not currently used

        e_x = (p2 - p1) / np.linalg.norm(p2 - p1)

        i = np.dot(e_x, (p3 - p1))

        e_y = (p3 - p1 - (i * e_x)) / (np.linalg.norm(p3 - p1 - (i * e_x)))
        e_z = np.cross(e_x, e_y)
        d = np.linalg.norm(p2 - p1)
        j = np.dot(e_y, (p3 - p1))
        x = ((r1 ** 2) - (r2 ** 2) + (d ** 2)) / (2 * d)
        y = (((r1 ** 2) - (r3 ** 2) + (i ** 2) + (j ** 2)) / (2 * j)) - ((i / j) * (x))
        z1 = (r1 ** 2 - x ** 2 - y ** 2)
        z2 = (r1 ** 2 - x ** 2 - y ** 2)
        if (z1 >=0) and (z2 >=0):
            z1 = np.sqrt(z1)
            z2 = np.sqrt(z2)* (-1)
            ans1 = p1 + (x * e_x) + (y * e_y) + (z1 * e_z)
            ans2 = p1 + (x * e_x) + (y * e_y) + (z2 * e_z)
            dist1 = np.linalg.norm(centroid - ans1)
            dist2 = np.linalg.norm(centroid - ans2)
            if (dist1 < dist2):
                return ans1
            else:
                return ans2
        else:
            return None

