import motorgillespie.simulation.initiate_motors as im
import motorgillespie.simulation.gillespie_simulation as gsim
import pickle
import os


def init_run(sim_params, gill_set, *motor_params, dirct, subdir, sd=None):
    """ This function handles the motor_team and fixed_motor initiating,
    execution of the Gillespie simulation and the saving of the resultant objects.

    For more information on the required dictionaries, see the 'templates' folder
    of this project.

    Parameters
    ----------
    sim_params : dictionary
        Simulation parameters are required arguments for the 'init_motor_fixed'
        function and include the constant external force and strap stiffness.
    gill_set : dictionary
        Gillespie settings required by the 'gillespie_2D_walk' function,
        including 't_max', 'n_runs', 'dimension', and 'single_run'. 'n_motors'
        is needed by the 'init_mixed_team' function to create the appropriate amount of
        each included motor species, and should be parsed as a list.
    *motor_params : dictionary
        Motor parameters to create user-defined motor species. Each dictionary
        corresponds to the parameters of a motor protein species. The order
        should match the numbers in the 'n_motors' list.

        Example: When gill_set['n_motors'] = [1,3], init_run(sim_params, gill_set, motor_params_1,
        motor_params_2, dirct, subdir, sd=None) will create, by calling 'init_mixed_team', one motor
        with parameters motor_params_1 and three with motor_params_2.
    dirct : string
        Name of the subdirectory to be created within the 'motor_objects'
        directory. It should be within the current working directory, which
        should be this project or a subdirectory of this project.
    subdir : string
        Name of the subdirectory within 'motor_objects/dirct', which represents
        a single simulation. It can be used when looping through values of
        specific simulation parameters, Gillespie settings, or motor parameters.
        See the 'end_report' directory for an example of simulation structure.
    sd : string
        Short description of the simulation for in the metadata file.

    """

    # Set local variables
    n_motors = gill_set['n_motors']
    t_max = gill_set['t_max']
    n_it = gill_set['n_it']
    dimension = gill_set['dimension']
    single_run = gill_set['single_run']

    # Create subdirectories in current working directory
    if not os.path.isdir(f'.\motor_objects\{dirct}'):
        os.makedirs(f'.\motor_objects\{dirct}')
    os.makedirs(f'.\motor_objects\{dirct}\{subdir}')

    ### Create motor_team of motor proteins ###
    print('Initiating motor motor_team..')
    motor_team = im.init_mixed_team(n_motors, *motor_params)
    ### Create fixed motor ###
    print('Initiating motor_fixed...')
    motor_fixed = im.init_motor_fixed(sim_params)

    ##############################################
    ### Simulate motor dynamics with Gillespie ###
    ##############################################
    print('Call Gillespie simulation...')
    team_out, motor_fixed_out = gsim.gillespie_2D_walk(motor_team=motor_team, motor_fixed=motor_fixed, t_max=t_max, n_runs=n_it, dimension=dimension, single_run=single_run)
    print('Done simulating')

    del motor_team
    del motor_fixed

    ### Pickle (linearize) objects ###
    print('Pickling motor motor_team...')
    while team_out:
        print(f'{team_out[-1].id}_{team_out[-1].family}_{team_out[-1].direction}')
        pickleTeam = open(f'.\motor_objects\{dirct}\{subdir}\{team_out[-1].id}_{team_out[-1].family}_{team_out[-1].direction}', 'wb')
        pickle.dump(team_out[-1], pickleTeam)
        pickleTeam.close()
        team_out.pop()

    print('Pickling motor_fixed...')
    pickleMotorFixed = open(f'.\motor_objects\{dirct}\{subdir}\motor0', 'wb')
    pickle.dump(motor_fixed_out, pickleMotorFixed)
    pickleMotorFixed.close()
    del motor_fixed_out

    ### Write metadata file in subdirectory ###
    print('Saving metadata to .txt filename...')
    with open(f".\motor_objects\{dirct}\{subdir}\parameters.txt", "w") as par_file:
        par_file.write(f"Simulation description: {sd}: \n")
        for index, dict in enumerate(motor_params):
            par_file.write(f"Motor{index+1} parameters: \n")
            for mp, value in dict.items():
                par_file.write(f"{mp}={value}\n")
        par_file.write("Simulation parameters: \n")
        for sp, value in sim_params.items():
            par_file.write(f"{sp}={value}\n")
        par_file.write("Gillespie settings: \n")
        for gs, value in gill_set.items():
            par_file.write(f"{gs}={value}\n")
    print('Done')

    return
