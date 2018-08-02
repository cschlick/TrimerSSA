from TrimerTemplate import TrimerTemplate
import random
import numpy as np
import scipy.stats as stats

class TrimerGenerator:


    def __init__(self,template_list,options):
        self.template_list = template_list
        self.weights_list = self.__generate_weights_list()
        self.type_dict = self.__generate_type_dict()
        # initialize another gamma distribution to decide whether an edge is too long
        alpha, loc, beta, N = 1, 3, options["edge_length_beta"], 10000
        self.padding_range = stats.gamma.rvs(alpha, loc=loc, scale=beta, size=N)

    def choose(self,adding_edge,trimer_type_request=None):
        if trimer_type_request is None:
            found_right_base_length = False
            search_counter = 0
            padding = np.random.choice(self.padding_range)
            low = adding_edge.length - padding
            high = adding_edge.length + padding
            range = high-low
            while (found_right_base_length is False) and (search_counter < 20):   # try to randomly get a trimer with correct adding_edge length
                choice = np.random.choice(self.template_list,1,p=self.weights_list)[0]
                if  (low <=choice.adding_edge_length<=high):
                    found_right_base_length = True
                search_counter+=1
            if (found_right_base_length is True) and type(choice) is TrimerTemplate:
                return choice
            else:                                                               # if we did not find one, try to brute force fine one
                for template in self.template_list:
                    if (low<=template.adding_edge_length<=high):
                        choice = template
                        found_right_base_length = True
            if (found_right_base_length is True):
                return choice

            else:                                                                  # can't do it, will return nothing which should become a rejection object
                #print("No template with correct adding edge length found, adding edge length:",adding_edge.length," low:",low," high:",high," range:",range)
                template_lengths = [t.adding_edge_length for t in self.template_list]
                #print("Template adding edge lengths:",template_lengths)

                return None

        else:
            return self.type_dict[trimer_type_request]



    def __generate_weights_list(self):
        weights_list = []
        sum_weights = 0
        for template in self.template_list:
            sum_weights+= template.weight
        for template in self.template_list:
            weight = template.weight/sum_weights
            weights_list.append(weight)
        return weights_list

    def __generate_type_dict(self):
        type_dict = {}
        for template in self.template_list:
            type_dict[template.template_type] = template
        return type_dict
