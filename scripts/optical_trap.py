import motorgillespie.simulation.variable_loops as vl
import time
import cProfile
import pstats

### Motor parameters ###
plus_params = {
 'family': 'Kinesin-1',
 'member': 'antero',
 'step_size': 8,
 'k_m': None,
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

### Simulation parameters ###
sim_params = {
 'k_t': None,
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
init_state = 'unbound'

date = time.strftime("%Y%m%d_%H%M%S")
dirct = f'{date}_optical_trap_{t_end}_{singlerun}_{init_state}'

team_comb = [[1], [2], [3], [4]]
kt = [0.02, 0.05, 0.1] # this isn't the order in which it is saved within the filesystem!
km = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4] # this isn't the order in which it is saved within the filesystem!

# Simulate per combination of parameters
for i in team_comb:
    for j in kt:
        for k in km:
            gill_set['n_motors'] = i
            sim_params['k_t'] = j
            plus_params['k_m'] = k

            ### Data storage ###
            subdir = f'{i}n_{j}kt_{k}km'
            short_description = ''

            # Initiate motor motor_team an run simulation n_it times for t_end seconds each
            with cProfile.Profile() as profile:
                out = output_gillespie = vl.init_run(sim_params, gill_set, plus_params, sd=short_description, dirct=dirct, subdir=subdir)
                ps = pstats.Stats(profile)
                ps.sort_stats('calls', 'cumtime')
                ps.print_stats(15)
