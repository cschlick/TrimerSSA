from Particle import Particle
from Trimer import Trimer
from Edge import Edge
from Vertex import Vertex
from Writer import Writer
from TrimerGenerator import TrimerGenerator
from TrimerTemplate import TrimerTemplate
from Timestep import Timestep
import random
import hashlib
import numpy as np

options = {
    "clash_tolerance_top_top": 150,
# top points closer than this value will clash with each other (excluding closely connected neighbors)
    "clash_tolerance_top_vertex": 70,
# top points closer than this value to another vertex will clash (excluding closely connected neighbors)
    "clash_tolerance_vertex": 60,
# two vertex points closer than this value (and greater than the merge tolerance) will clash
    "merge_tolerance_beta": 4,
# vertex points closer a cutoff value will trigger a merge. The value is chosen from a gamma distribution with this beta value
    "edge_length_beta": 5,
# Trying to add a template to an adding edge that is outside of a tolerance will fail. That tolerance can be adjusted here, another gamma distribution
    "centroid_update_interval": 200
# the number of timesteps between recalculating the particle centroid (used to determine the curvature direction)
}


# manually give the coordinates of a seed trimer, and initialize a particle object
def seed(trimer_generator, options):
    seed_vertex1 = Vertex(np.array([198.829, 170.530, 360.401]))
    seed_vertex2 = Vertex(np.array([161.790, 244.135, 376.374]))
    seed_vertex3 = Vertex(np.array([124.752, 199.392, 315.765]))
    seed_edge1 = Edge([seed_vertex1, seed_vertex2])
    seed_edge2 = Edge([seed_vertex2, seed_vertex3])
    seed_edge3 = Edge([seed_vertex3, seed_vertex1])
    seed_trimer = Trimer([seed_edge1, seed_edge2, seed_edge3])
    seed_center = np.array([231.149, 245.5, 262.1035])
    particle = Particle(seed_trimer, trimer_generator=trimer_generator, options=options)
    particle.centroid = seed_center
    return particle


# define trimer add types
type1 = TrimerTemplate(angle_degrees=171, stem_length=71, adding_edge_length=91, template_type=1, weight=1)
type2 = TrimerTemplate(angle_degrees=149, stem_length=77, adding_edge_length=84, template_type=2, weight=1)
type3 = TrimerTemplate(angle_degrees=149, stem_length=77, adding_edge_length=84, template_type=3, weight=1)
type4 = TrimerTemplate(angle_degrees=171, stem_length=79, adding_edge_length=91, template_type=4, weight=3)
type5 = TrimerTemplate(angle_degrees=148, stem_length=77, adding_edge_length=84.5, template_type=5, weight=1)
type6 = TrimerTemplate(angle_degrees=148, stem_length=77, adding_edge_length=84.5, template_type=6, weight=1)
type7 = TrimerTemplate(angle_degrees=180, stem_length=77, adding_edge_length=88, template_type=7,
                       weight=3)  # flat hexamer

# collect the types we want to use together in an object
trimer_generator = TrimerGenerator([type1, type2, type3, type4, type5, type6, type7], options)

# execute the seed function (above) to get a particle object
particle = seed(trimer_generator, options)

## simulation
prob_max = 1500  # the range in which to select random numbers to decide add/remove. If the prob_on==prob_max, every step will try to add
prob_on = 1000  # the probability
on_off_ratio_remove_single = 1 / 0.6
on_off_ratio_remove_double = 1 / 0.01
prob_off_single = prob_on * (1 / on_off_ratio_remove_single)
prob_off_double = prob_on * (1 / on_off_ratio_remove_double)

print("Prob on:", prob_on)
print("Prob of single bonded trimer removal:", prob_off_single)
print("Prob of double bonded trimer removal:", prob_off_double)
print("All probs out of ", prob_max)


def simulate(particle, steps):
    print("Running ", steps, " timesteps")
    rands_on = np.random.randint(0, prob_max, steps)
    rands_off = np.random.randint(0, prob_max, steps)

    complete_flag = False
    for i, rand_on in enumerate(rands_on):

        if particle.complete is False:
            tried_to_add = False
            tried_to_remove_single = False
            tried_to_remove_double = False
            add_outcome = None
            remove_single_outcome = None
            remove_double_outcome = None

            if (rand_on <= prob_on):
                add_outcome = particle.add()
                tried_to_add = True

            if (rands_off[i] <= prob_off_single):
                tried_to_remove_single = True
                if (len(particle.open_trimers) > 2) and (len(particle.single_bond_trimers) > 0):
                    removing_trimer = random.choice(particle.single_bond_trimers)
                    remove_single_outcome = particle.remove(removing_trimer)

            if (rands_off[i] <= prob_off_double):
                tried_to_remove_double = True
                if (len(particle.open_trimers) > 2) and (len(particle.double_bond_trimers) > 0):
                    removing_trimer = random.choice(particle.double_bond_trimers)
                    remove_double_outcome = particle.remove(removing_trimer)

            timestep = Timestep(particle.timestep,
                                len(particle.trimers),
                                add_outcome,
                                remove_single_outcome,
                                remove_double_outcome,
                                tried_to_add,
                                tried_to_remove_single,
                                tried_to_remove_double)

            particle.increment_timestep()
            particle.timesteps.append(timestep)
        elif complete_flag is False:
            complete_flag = True
            print("Complete particle")
        else:
            pass

    return particle


###############################################################################################################################################

# start new particle
particle = seed(trimer_generator,options)
particle = simulate(particle,1000)



particle.write("output/",summarize=True)



