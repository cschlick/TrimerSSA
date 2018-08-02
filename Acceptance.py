

class Acceptance: # the alternative to creating a rejection object
    def __init__(self,timestep,add=False,remove=False,template_type=None,merging=None):
        self.timestep = timestep
        self.add=add
        self.remove = remove
        self.template_type = template_type
        self.merging = merging