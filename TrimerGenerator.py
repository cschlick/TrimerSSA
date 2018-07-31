import random
from numpy.random import choice

class TrimerGenerator:
    def __init__(self,template_list):
        self.template_list = template_list
        self.weights_list = self.__generate_weights_list()
        self.type_dict = self.__generate_type_dict()

    def choose(self,trimer_type_request=None):
        if trimer_type_request is None:
            return choice(self.template_list,1,p=self.weights_list)[0]

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
