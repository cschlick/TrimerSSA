class Rejection:


    rejection_types = {
        1: "Add Error: vertex-vertex clash",
        2: "Add Error: top-vertex clash",
        3: "Add Error: top-top clash",
        4: "Add Error: tip-tip merge",
        5: "Add Error: unhandled number of hinge vertices",
        6: "Add Error: self clash",
        7: "Remove Error: number of open eges in trimer is not 1 or 2",
        8: "Add Error: too long of edge for any trimer templates",

    }


    def __init__(self,timestep,rejection_type,clash_dist=None,clash_partner=None,text=None):
        self.timestep = timestep
        self.rejection_type = rejection_type
        self.clash_distance = clash_dist
        self.clash_partner = clash_partner
        self.text = text

        # types of rejection:

