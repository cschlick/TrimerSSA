

class Timestep: # organizational tool. Should generate one of these objects each time step
    def __init__(self,
                 step,
                 trimer_number,
                 add_outcome,
                 remove_single_outcome,
                 remove_double_outcome,
                 tried_to_add,
                 tried_to_remove_single,
                 tried_to_remove_double):

        self.step = step
        self.trimer_number = trimer_number
        self.add = add_outcome
        self.remove_single = remove_single_outcome
        self.remove_double = remove_double_outcome

        self.tried_to_add = tried_to_add
        self.tried_to_remove_single = tried_to_remove_single
        self.tried_to_remove_double = tried_to_remove_double