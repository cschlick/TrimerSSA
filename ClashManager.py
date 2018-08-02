from Rejection import Rejection
import numpy as np
from transformations import unit_vector
from Top import Top
from Vertex import Vertex
from Edge import Edge



class ClashManager:
    def __init__(self,particle,new_vertex,merge_distribution):
        self.particle = particle
        self.new_vertex = new_vertex
        self.clash_tolerance_top_top = particle.options["clash_tolerance_top_top"]
        self.clash_tolerance_top_vertex = particle.options["clash_tolerance_top_vertex"]
        self.clash_tolerance_vertex = particle.options["clash_tolerance_vertex"]
        # initialize the merge tolerance distribution
        self.merge_distribution = merge_distribution

    @property
    def merge_tolerance(self):
        return np.random.choice(self.merge_distribution)




    def check_clash_vertex(self, candidate_vertex,
                           adding_edge):  # return Bool,Bool,object for Clash, Merge, Merge Vertex
        rejection = None
        clash_flag = False
        merge_flag = False
        merging_vertex = None

        for trimer in self.particle.trimers:

            (e1, e2, e3) = trimer.edges
            (v1, v2, v3) = trimer.verteces
            # dist1 = np.linalg.norm(trimer.coord - candidate_vertex.coord)
            # dist2 = np.linalg.norm(e1.coord - candidate_vertex.coord)
            # dist3 = np.linalg.norm(e2.coord - candidate_vertex.coord)
            # dist4 = np.linalg.norm(e3.coord - candidate_vertex.coord)
            dist5 = np.linalg.norm(v1.coord - candidate_vertex.coord)
            dist6 = np.linalg.norm(v2.coord - candidate_vertex.coord)
            dist7 = np.linalg.norm(v3.coord - candidate_vertex.coord)
            tolv = self.clash_tolerance_vertex
            tolm = self.merge_tolerance
            # check for clashes
            # if (dist1 <= tol):
            #     rejection = Rejection("vertex clash with trimer center")
            #     clash_flag = True
            # if (dist2 <= tol) and (e1 not in disallowed_edges):
            #     rejection = Rejection("vertex clash with edge")
            #     clash_flag = True
            # elif (dist3 <= tol) and (e2 not in disallowed_edges):
            #     rejection = Rejection("vertex clash with edge")
            #     clash_flag = True
            # elif (dist4 <= tol) and (e3 not in disallowed_edges):
            #     rejection = Rejection("vertex clash with edge")
            #     clash_flag = True
            if (tolm < dist5 <= tolv):
                rejection = Rejection(self.particle.timestep,1,clash_dist=dist5,clash_partner=v1)
                clash_flag = True
            elif (tolm < dist6 <= tolv):
                rejection = Rejection(self.particle.timestep,1,clash_dist=dist6,clash_partner=v2)
                clash_flag = True
            elif (tolm < dist7 <= tolv):
                rejection = Rejection(self.particle.timestep,1,clash_dist=dist6,clash_partner=v3)
                clash_flag = True

            if (dist5 <= tolm) and (clash_flag == False):  # check for merges
                merging_vertex = v1
                merge_flag = True
            elif (dist6 <= tolm) and (clash_flag == False):
                merging_vertex = v2
                merge_flag = True
            elif (dist7 <= tolm) and (clash_flag == False):
                merging_vertex = v3
                merge_flag = True

        return clash_flag, merge_flag, merging_vertex, rejection

    def check_clash_tops(self,adding_edge,candidate_vertex):
        adding_trimer = next(iter(adding_edge.trimers))

        # make tops
        (v1,v2) = adding_edge.verteces
        v3 = candidate_vertex
        n1 = v2.coord - v1.coord
        n2 = v3.coord - v1.coord
        vec = np.cross(n1, n2)
        d = 25
        center = (adding_edge.coord+candidate_vertex.coord)/2
        top1 = Top(center + d * unit_vector(vec))
        top2 = Top(center - d * unit_vector(vec))
        new_tops = {top1, top2}

        # get a set of trimers to check for clashes with (not closesly connected to current trimer)
        check_trimers = self.exclude_nearby_trimers(adding_trimer)

        # check for clashes
        top_clash_flag = False
        rejection = None
        for trimer in check_trimers:
            if top_clash_flag is False:

                (v1,v2,v3) = trimer.verteces
                (t1,t2) = trimer.tops

                results = [] # a list of True false, if any are True, there was a clash
                toltt = self.clash_tolerance_top_top
                toltv = self.clash_tolerance_top_vertex
                if results.append(self.check_clash_between_two_objects(top1,v1, toltv)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top1, v2, toltv)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top1, v3, toltv)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top1, t1, toltt)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top1, t2, toltt)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top2, v1, toltv)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top2, v2, toltv)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top2, v3, toltv)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top2, t1, toltt)):
                    top_clash_flag = True
                elif results.append(self.check_clash_between_two_objects(top2, t2, toltt)):
                    top_clash_flag = True

        if top_clash_flag is True:
            rejection = Rejection(self.particle.timestep,3)

        return top_clash_flag, rejection

    def check_clash_between_two_objects(self,obj1,obj2,tolerance):
        dist = np.linalg.norm(obj1.coord - obj2.coord)
        if (dist <= tolerance):
            return True
        else:
            return False

    def exclude_nearby_trimers(self,adding_trimer): # returns a set of trimers to actually check

        exclude_trimers = set()


        exclude_trimers.add(adding_trimer)

        # collect all secondary trimers 5 layers deep
        for secondary1 in adding_trimer.secondary_trimers:
            for secondary2 in secondary1.secondary_trimers:
                for secondary3 in secondary2.secondary_trimers:
                    for secondary4 in secondary3.secondary_trimers:
                        for secondary5 in secondary4.secondary_trimers:
                            exclude_trimers.add(secondary5)
                        exclude_trimers.add(secondary4)
                    exclude_trimers.add(secondary3)
                exclude_trimers.add(secondary2)
            exclude_trimers.add(secondary1)

        # make set of objects to actually check
        check_trimers = self.particle.trimers - exclude_trimers

        return check_trimers
