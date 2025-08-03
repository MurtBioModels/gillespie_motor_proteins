import motorgillespie.simulation.variable_loops as vl
import time
import cProfile
import pstats

### Motor parameters ###
plus_params = {
 'family': 'Kinesin-1',
 'member': 'plus',
 'step_size': 8,
 'k_m': 0.2,
 'v_0': 740,
 'alfa_0': 92.5,
 'f_s': 7,
 'epsilon_0': 0.66,
 'f_d': 2.1,
 'bind_rate': 5,
 'direction': 'anterograde',
 'init_state': 'unbound',
 'calc_eps': 'exponential',
}

minus_params = {
 'family': 'Kinesin-1',
 'member': 'minus',
 'step_size': 8,
 'k_m': None,
 'v_0': 740,
 'alfa_0': 92.5,
 'f_s': 7,
 'epsilon_0': 0.66,
 'f_d': 2.1,
 'bind_rate': 5,
 'direction': 'retrograde',
 'init_state': 'unbound',
 'calc_eps': 'exponential',
}


### Simulation parameters ###
sim_params = {
 'k_t': 0,
 'f_ex': 0
}

### Simulation settings ###
gill_set = {
    'n_motors': None,
    'n_it': 1000, # set this to 1 to test if the scrip works without errors
    't_max': 100, # set this to 10 to test if the scrip works without errors
    'dimension': '1D',
    'single_run': False
}

n_motors = gill_set['n_motors']
t_end = gill_set['t_max']
singlerun = gill_set['single_run']
init_state = 'notbound'
date = time.strftime("%Y%m%d_%H%M%S")
dir = f'{date}_{t_end}_{singlerun}_{init_state}'

team_comb = [[1, 1], [2, 2], [3, 3], [4, 4]] # this isn't the order in which it is saved within the filesystem!
retro_km = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4] # this isn't the order in which it is saved within the filesystem!

# Simulate per combination of parameters
for i in team_comb:
    for j in retro_km:
        gill_set['n_motors'] = i
        minus_params['k_m'] = j

        ### Data storage ###
        subdir = f'{i}n_{j}minuskm'
        short_description = ''

        # Initiate motor motor_team an run simulation n_it times for t_end seconds each
        with cProfile.Profile() as profile:
            out = output_gillespie = vl.init_run(sim_params, gill_set, minus_params,  plus_params,  sd=short_description, dirct=dir, subdir=subdir)
            ps = pstats.Stats(profile)
            ps.sort_stats('calls', 'cumtime')
            ps.print_stats(15)
