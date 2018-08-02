import numpy as np

class TrimerTemplate:
    def __init__(self,angle_degrees=0,stem_length=10,adding_edge_length=10,template_type=0,weight=1):
        self.stem_length = stem_length
        self.adding_edge_length=adding_edge_length
        self.angle_degrees = angle_degrees
        self.angle_radians = angle_degrees*(np.pi/180)
        self.template_type=template_type
        self.weight=weight