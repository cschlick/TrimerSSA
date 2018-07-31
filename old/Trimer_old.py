
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
        self.sites = set() # set of sets, associating each vertex with a site


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
            d = 84*0.3
            top1 = Top(self.coord + d * unit_vector(vec))
            top2 = Top(self.coord - d * unit_vector(vec))
            top1.trimer = self
            top2.trimer = self
            tops = set([top1,top2])
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

    def add(self, particle, template,adding_edge):  # returns a candidate new vertex coord
        #print("Adding trimer template type:",template.template_type)
        angle = template.angle
        stem_length = template.stem_length

        assert (self.full == False)
        # edge = adding_edge
        # opposite_vertex = self.opposite_vertex(edge)
        # v1 = edge.vertex1.coord
        # v2 = edge.vertex2.coord
        # v3 = opposite_vertex.coord
        # v4 = self.coord
        # a = edge.coord
        #

        #print("adding template type:",template.template_type)


        E1 = adding_edge
        (V1, V2) = adding_edge.verteces
        V3 = self.opposite_vertex(adding_edge)

        # # maybe extend length a bit along the V1V3 vector to get the appropriate stem length for the template

        vec = E1.coord-V3.coord

        V1V3_dist = np.linalg.norm(V3.coord - E1.coord)
        #print(V1V3_dist)

        d = stem_length/V1V3_dist # this does not extend

        #print(V1V3_dist,stem_length,d)
        extended_point1 = E1.coord - d * vec
        extended_point2 = E1.coord + d * vec
        dist1 = np.linalg.norm(extended_point1-V3.coord)
        dist2 = np.linalg.norm(extended_point2-V3.coord)
        if dist1<dist2:
            extended_point = extended_point1
        else:
            extended_point = extended_point2

        # set origin to E1
        e1 = E1.coord
        v1 = V1.coord-e1
        v2 = V1.coord-e1
        v3 = V3.coord-e1
        extended_point = extended_point-e1




        #########
        def get_new_point(reference_point,start_point,angle,axis):  # where angle is radians of the stem length (edge to opposite vertex)

            M = rotation_matrix(angle, axis)
            new_point = np.dot(start_point, M[:3, :3].T)
            return new_point+reference_point

        # get two new points with consistent angle for first rotation
        axis = V1.coord-V2.coord

        axis = unit_vector(axis)



        p1 =get_new_point(e1,extended_point,angle,axis)
        p2 = get_new_point(e1, extended_point,-1*angle, axis)
        #debug
        # dbgv1 = Vertex(p1)
        # dbgv2 = Vertex(p2)
        # particle.debug_verteces.add(dbgv1)
        # particle.debug_verteces.add(dbgv2)

        # pick closest to centroid
        if np.linalg.norm(p1 - particle.centroid) < np.linalg.norm(p2 - particle.centroid):
            new_point = p1

        else:
            new_point = p2


        new_vertex = Vertex(new_point)

        return new_vertex # returns the adding_edge, new_vertex

    def add2(self, particle, template,adding_edge):  # returns a candidate new vertex coord
        angle1 = template.angle1
        angle2 = template.angle2
        stem_length = template.stem_length

        assert (self.full == False)



    #     V4,V5,V6,V7
    #
    #          V1-------V3
    #         /  \     /
    #       /     \  /
    #      V7------V2
    #
    # E1 is midpoint of V1,V2 (the adding edge)
    # V4 is the first point on the way to the final point (V7). V4 is 1 stem_length above E1, 90 degrees from the Plan(V1,V2,V3)
    # V5 will be an intermeidate point after the first rotation
    # A 90 degree rotation of V5 each way will give V6, a test point to see which direction to rotate
    # A final rotation of V5 will yield V7

        E1 = adding_edge
        (V1,V2) = adding_edge.verteces
        V3 = self.opposite_vertex(adding_edge)

        # debug
        V1.atom_name = "Be"
        V2.atom_name = "Ba"
        V3.atom_name = "Ca"
        adding_edge.atom_name = "Na"


        vec_V2V3 = V2.coord - V3.coord
        vec_V2V1 = V2.coord - V1.coord
        vec_normal = np.cross(vec_V2V1,vec_V2V3)

        V4_1 = E1.coord + (stem_length * unit_vector(vec_normal))
        V4_2 = E1.coord - (stem_length * unit_vector(vec_normal))

        dbgv = Vertex(V4_1)
        dbgv.atom_name = "S"
        particle.debug_verteces.add(dbgv)
        dbgv = Vertex(V4_2)
        dbgv.atom_name = "O"
        particle.debug_verteces.add(dbgv)



        dist1 = np.linalg.norm(V4_1-particle.centroid)
        dist2 = np.linalg.norm(V4_2-particle.centroid)
        if (dist1>=dist2): # take the value further from centroid
            V4 = V4_1
            print("chose V4_1")
        else:
            V4 = V4_2
            print("chose V4_2")



        # set origin to E1
        e1 = E1.coord
        v1 = V1.coord-e1
        v2 = V2.coord-e1
        v3 = V3.coord-e1
        v4 = V4-e1


        #########
        def get_new_point(reference_point,start_point,angle,axis):  # where angle is radians of the stem length (edge to opposite vertex)

            M = rotation_matrix(angle, axis)
            new_point = np.dot(start_point, M[:3, :3].T)
            return new_point+reference_point

        # get two new points with consistent angle for first rotation
        axis = E1.coord-V3.coord

        axis = unit_vector(axis)


        V5_1 = get_new_point(e1,v4,angle1,axis)
        V5_2 = get_new_point(e1, v4, (2*np.pi)-angle1, axis)

        dbgv = Vertex(V5_1)
        dbgv.atom_name = "Mg"
        particle.debug_verteces.add(dbgv)
        dbgv = Vertex(V5_2)
        dbgv.atom_name = "Al"
        particle.debug_verteces.add(dbgv)


        #perform a test rotation of V5 to get V6 and determine where to rotate towards for V7
        axis = V2.coord - V1.coord
        V5 = V5_1 # chooseing V5_1
        v5 = V5-e1 # at some point will need to choose which V5 to use. (random choice?)




        a1 = (0.5) * np.pi
        a2 = (-0.5) * np.pi


        V6_1 = get_new_point(e1,v5,a1,axis)
        V6_2 = get_new_point(e1, v5, a2, axis)

        dbgv = Vertex(V6_1)
        dbgv.atom_name = "Cl"
        particle.debug_verteces.add(dbgv)
        dbgv = Vertex(V6_2)
        dbgv.atom_name = "Na"
        particle.debug_verteces.add(dbgv)

        dist1 = np.linalg.norm(V6_1-V3.coord)
        dist2 = np.linalg.norm(V6_2-V3.coord)
        print(dist1,dist2)
        if dist1>dist2:
            V6 = V6_1
            print("chose V6_1")
            sign = 1
        else:
            V6 = V6_2
            print("chose V6_2")
            sign = -1

        # make final rotation of V5 to V7
        a1 = (0.5 + (sign*angle2)) * np.pi
        a2 = (1 + (sign*angle2)) * np.pi

        V7_1 = get_new_point(e1, v5, a1, axis)
        V7_2 = get_new_point(e1, v5, a2, axis)

        dist1 = np.linalg.norm(V6-V7_1)
        dist2 = np.linalg.norm(V6-V7_2)

        if dist1 < dist2:
            V7 = V7_1
        else:
            V7 = V7_2

        dbgv = Vertex(V7_1)
        dbgv.atom_name = "Mn"
        particle.debug_verteces.add(dbgv)
        dbgv = Vertex(V7_2)
        dbgv.atom_name = "Li"
        particle.debug_verteces.add(dbgv)


        new_vertex = Vertex(V7)

        return new_vertex  # returns the adding_edge, new_vertex




    def add3(self,particle, template,adding_edge):
        (v1,v2) = adding_edge.verteces
        e1 = adding_edge
        V3 = self.opposite_vertex(adding_edge)
        p1 = v1.coord
        r1 = template.v1_edge_length
        p2 = v2.coord
        r2 = template.v2_edge_length
        p3 = V3.coord
        r3 = template.v3_edge_length
        centroid = particle.centroid
        new_vertex_coord = self.trilaterate3D(p1, r1, p2, r2, p3, r3, centroid)
        if new_vertex_coord is not None:
            new_vertex = Vertex(new_vertex_coord)
            return new_vertex
        else:
            return False





    def trilaterate3D(self, p1, r1, p2, r2, p3, r3, centroid):

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

