import numpy as np
from Trimer import Trimer
import random

from Edge import Edge
from Vertex import Vertex
from Writer import Writer
from Top import Top
from Rejection import Rejection
from ClashManager import ClashManager
from lib import assertlength
from transformations import unit_vector

class Particle:
    def __init__(self, seed_trimer,trimer_generator):
        """
        Particles mostly contain a set of trimers
        """
        self.trimers = set()
        self.trimers.add(seed_trimer)

        # integer accounting for pdb connect records
        self.available_intids = [x for x in range(1, 9999)]  # set of possible ids
        self.intid_dict = {}  # maps id(object) to intid

        self.rejections = []  # list of rejection objects (hold reasons for rejection)

        self.debug_verteces = set()
        self.centroid = np.array([0, 0, 0])  # always has intid 0

        self.curvature_history = []
        # self.stem_length = np.sqrt((self.edge_length ** 2) - ((self.edge_length / 2) ** 2))
        # self.curvaturefunc = curvaturefunc  # in radians of

        self.history = []

        self.trimer_generator = trimer_generator
        self.__timestep = 1
        self.add_events = 1
        self.remove_events_single_bond=0
        self.remove_events_double_bond=0




    @property
    def timestep(self):
        return self.__timestep


    @property
    def open_trimers(self):
        return [trimer for trimer in self.trimers if trimer.full == False]

    @property
    def single_bond_trimers(self):  # trimers with one open edges
        return [trimer for trimer in self.trimers if len(trimer.open_edges) == 2]
    @property
    def double_bond_trimers(self): # trimers with two open edges
        return [trimer for trimer in self.trimers if len(trimer.open_edges) == 1]

    @property
    def open_edges(self):
        return [edge for edge in self.edges if edge.full == False]

    @property
    def complete(self):
        if self.num_open_edges == 0:
            return True
        else:
            return False

    @property
    def num_open_edges(self):
        num = 0

        for edge in self.edges:
            if edge.full == False:
                num += 1
        return num

    @property
    def edges(self):
        myedges = set()


        for trimer in self.trimers:
            myedges = myedges | trimer.edges
        return myedges

    @property
    def verteces(self):
        verteces = set()
        for trimer in self.trimers:
            verteces = verteces | trimer.verteces
        return verteces

    @property
    def tops(self):
        tops = set()
        for trimer in self.trimers:
            tops = tops | trimer.tops
        return tops

    @property
    def objects(self):  # verteces and edges for now
        return self.edges | self.verteces | self.debug_verteces


    def increment_timestep(self):
        self.history.append(len(self.trimers))
        # check for single subunit discontinuity
        if len(self.trimers) > len(self.open_trimers):
            remove_set = {trimer for trimer in self.trimers if len(trimer.open_edges)==3}
            for trimer in remove_set:
                self.trimers.remove(trimer)

        self.__timestep+=1



    def intid_update(self):

        for object in self.objects:
            if id(object) not in self.intid_dict.keys():
                intid = random.choice(self.available_intids)
                self.intid_dict[id(object)] = intid
                self.available_intids.remove(intid)





    def merge(self,adding_edge,merging_vertex):

            # find hinge vertex...
            hinge_vertex_options_set = adding_edge.verteces.intersection(merging_vertex.secondary_verteces)



            if len(hinge_vertex_options_set) == 0:
                #print("probably a tip-tip merge, rejecting...")
                # for now, reject tip tip merges
                return Rejection("tip tip merge")

            elif len(hinge_vertex_options_set) > 1: # probably closing a hole. Need to merge two edges
                    #print("merging two edges, probably closing hole.")

                    verteces = list(hinge_vertex_options_set)
                    hinge_vertex = verteces[0]
                    hinge_opposite_vertex = verteces[1]
                    merging_edge_options_set = set()
                    for edge in hinge_vertex.edges:
                        if merging_vertex in edge.verteces:
                            merging_edge_options_set.add(edge)

                    assertlength(merging_edge_options_set, 1)

                    merging_edge = next(iter(merging_edge_options_set))
                    new_trimer = self.add_type_merge_two_edge(adding_edge,
                                                              merging_edge,
                                                              merging_vertex,
                                                              hinge_vertex,
                                                              hinge_opposite_vertex)
                    return new_trimer

            elif len(hinge_vertex_options_set) == 1: # probably adding with a single edge merge
                #print("Adding with single edge merge")
                hinge_vertex = next(iter(hinge_vertex_options_set))
                merging_edge_options_set = set()
                for edge in hinge_vertex.edges:
                    if merging_vertex in edge.verteces:
                        merging_edge_options_set.add(edge)

                assertlength(merging_edge_options_set, 1)
                merging_edge = next(iter(merging_edge_options_set))

                hinge_opposite_vertex_options_set = adding_edge.verteces - set([hinge_vertex])
                hinge_opposite_vertex = next(iter(hinge_opposite_vertex_options_set))


                new_trimer = self.add_type_merge_one_edge(adding_edge,
                                                          merging_edge,
                                                          merging_vertex,
                                                          hinge_vertex,
                                                          hinge_opposite_vertex)
                return new_trimer
            else:
                print("Error: Unhandled number of hinge verteces")

                return Rejection("Error: Unhandled number of hinge verteces")

    def add_type_merge_one_edge(self,
                                adding_edge,
                                merging_edge,
                                merging_vertex,
                                hinge_vertex,
                                hinge_opposite_vertex):



        new_edge = Edge([hinge_opposite_vertex, merging_vertex])
        #print(adding_edge, merging_edge, new_edge)
        new_trimer = Trimer([adding_edge, merging_edge, new_edge])
        adding_edge.full = True
        merging_edge.full = True
        return new_trimer


    def add_type_merge_two_edge(self,adding_edge,merging_edge,merging_vertex,hinge_vertex,hinge_opposite_vertex):

        new_edge_vertex_set = adding_edge.verteces - set([hinge_vertex])
        assert (len(new_edge_vertex_set) == 1)
        new_edge_vertex = next(iter(new_edge_vertex_set))


        possible_final_edges = [edge for edge in merging_vertex.edges if hinge_opposite_vertex in edge.verteces]
        assert (len(possible_final_edges) == 1)
        final_edge = possible_final_edges[0]

        new_trimer = Trimer([adding_edge, merging_edge, final_edge])
        adding_edge.full = True
        merging_edge.full = True
        final_edge.full = True
        return new_trimer


    def add(self,specific_edge=False,trimer_type_request=None):

        self.add_events+=1
        rejection = None





        # get a random trimer to add to
        if specific_edge is False:
            adding_trimer = random.choice(self.open_trimers)
            adding_edge = random.choice(adding_trimer.open_edges)
        elif type(specific_edge) is Edge:
            assert(len(specific_edge.trimers)==1)
            (adding_trimer,) = specific_edge.trimers
            adding_edge = specific_edge
        else:
            print("Error with add( , specify_edge=")


        # ask the trimer for a candidate trimer

        template = self.trimer_generator.choose(trimer_type_request=trimer_type_request)

        new_vertex = adding_trimer.add(self, template, adding_edge=adding_edge)
        clashmanager = ClashManager(self,new_vertex)
        clash_flag, merge_flag, merging_vertex, rejection = clashmanager.check_clash_vertex(new_vertex,adding_edge)


        # check for clashes with tops
        if clash_flag is False:
            clash_flag, rejection = clashmanager.check_clash_tops(adding_edge,new_vertex)


        if (merge_flag == False) and (clash_flag == False): # adding without merge

            #print("Adding without merge")
            adding_edge_verteces = list(adding_edge.verteces)
            new_edge1 = Edge([adding_edge_verteces[0],new_vertex])
            new_edge2 = Edge([new_vertex,adding_edge_verteces[1]])
            new_trimer = Trimer([adding_edge,new_edge1,new_edge2])

            self.trimers.add(new_trimer)
            adding_edge.full = True


        elif (merge_flag == True) and (clash_flag == False):  # adding with merge
            # check for self merge
            assert (len(adding_edge.trimers) == 1)
            adding_trimer = next(iter(adding_edge.trimers))


            if merging_vertex not in adding_trimer.verteces:
                new_trimer = self.merge(adding_edge,merging_vertex)
                if type(new_trimer) is Rejection:
                    rejection = new_trimer
                elif type(new_trimer) is Trimer:
                    self.trimers.add(new_trimer)
            else:
                rejection = Rejection("self clash")

        elif (clash_flag == True):
            # the timestep resulted in rejection
            self.rejections.append(rejection)

        else:
            print("big problem..."," Merge flag:",merge_flag," Clash flag:",clash_flag)


    def remove(self,trimer):
        open_edges = trimer.open_edges
        if len(open_edges) == 1:
            self.remove_events_double_bond +=1
        elif len(open_edges) == 2:
            self.remove_events_single_bond+=1
        else:
            print("error, number of open eges in trimer is not 1 or 2")


        self.trimers.remove(trimer)
        trimer.delete()




    def summarize(self):
        print("Particle Summary:")
        print("\tTimesteps: ",self.timestep)
        print("\tAdd events: ",self.add_events)
        print("\tRemove events single bonded:",self.remove_events_single_bond)
        print("\tRemove events double bonded:", self.remove_events_double_bond)
        print("\tTrimers: ", len(self.trimers))
        print("\tEdges: ", len(self.edges))
        print("\tVerteces: ", len(self.verteces))
        print("\tDebug Verteces: ", len(self.debug_verteces))
        print("\tOpen trimers: ", len(self.open_trimers))
        print("\tOpen edges: ", len(self.open_edges))
        print("\tRejections: ", len(self.rejections))
