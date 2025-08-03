import motorgillespie.simulation.gs_functions as gsf
import motorgillespie.errors as em
import numpy as np


def gillespie_2D_walk(motor_team, motor_fixed, t_max=100, n_runs=1000, dimension='1D', single_run=False):
    """ Function that executes the Tau-leaping Gillespie algorithm for studying stochastic motor protein dynamics.

    This is a Tau-leaping Gillespie algorithm for simulating (collective) motor protein dynamics.
    This function relies on calls to functions in the 'gs_functions' module, and requires the parameters 'motor_team'
    and 'motor_fixed' to be created priorly by the 'initiate_motors' module, which create instances of
    the classes MotorProtein and FixedMotor. This can be done all in one function call by the 'init_run' function within the 'variable_loop' module
    file.

    For n_runs:
        While t < t_max:
            Calculate cargo position (motorgillespie.simulation.gs_functions.calc_force_1D/calc_force_2D)
            Calculate motor forces (motorgillespie.simulation.gs_functions.calc_force_1D/calc_force_2D)
            Calculate stepping and detachment rates
            Draw waiting time
            Draw event (motorgillespie.simulation.gs_functions.draw_event)
            Execute event
            Check how many motors are bound
            Repeat

    More information on the Gillespie algorithm: https://en.wikipedia.org/wiki/Gillespie_algorithm

    Parameters
    ----------
    motor_team : list of motorgillespie.simulation.motor_class.MotorProtein
        Team of motor protein objects.
        Collects per motor data like locations and forces over time.
    motor_fixed : motorgillespie.simulation.motor_class.MotorFixed
        Instance of class 'MotorFixed'.
        Fixed location at x = 0.
        Holds simulation parameters and collects cargo data.
    t_max : int
        Duration of one iteration of the Gillespie simulation. Default is 100.
    n_runs : int
        Number of Gillespie iterations. Default is 1000.
    dimension : str
        Dimensionality of the simulation. Possible values: '1D' or '2D'. Default is '1D'.
    single_run : bool
        If True, each iteration terminates when the cargo unbinds, creating an amount of cargo runs equal to n_iteration.
        If False, the run will continue until t_max is reached, allowing more than one cargo run per iteration and trajectories
        where the cargo falls back to the starting position x = 0. Default is False.

    Returns
    -------
    motor_team : list of motorgillespie.simulation.motor_class.MotorProtein
        Updated MotorProtein objects, hold motor data.
    motor_fixed : motorgillespie.simulation.motor_class.MotorFixed
        Updated MotorFixed object, holds cargo data.
    """


    # Checking simulation settings
    print('Checking Gillespie settings...')
    dimension_options = ['1D', '2D']
    if dimension not in dimension_options:
        raise ValueError("Invalid function. Expected one of: %s" % dimension_options)
    # Checking motor parameters
    print('Checking motor parameters...')
    for motor in motor_team:
        motor.valid_motor(dimension=dimension)

    # Set local variables and necessary class attributes once (for speed purposes)
    print('Setting local variables...')
    radius = motor_fixed.radius
    rest_length = motor_fixed.rest_length
    dp_v1 = motor_fixed.dp_v1
    dp_v2 = motor_fixed.dp_v2
    temp = motor_fixed.temp
    if radius and rest_length is not None:
        motor_fixed.calc_angle()
    k_t = motor_fixed.k_t
    x_motor0 = motor_fixed.x_m
    f_ex = motor_fixed.f_ex

    # Init cargo settings once:
    motor_fixed.init_valid_once(motor_team=motor_team)

    ####################################
    ## Simulate n_runs Gillespie runs ##
    ####################################
    print('Begin simulation...')
    for i in range(0, n_runs):
        print(f'{i}th iteration')
        # Time
        t = 0

        # Set motors and cargo (motor_fixed) to initial state and update list attributes
        motor_fixed.fixed_init()
        for motor in motor_team:
            motor.motor_init(dimension=dimension, f_ex=f_ex)
        # This turns True when the last motor detaches
        end_run = False

        #####################################################################
        ### Gillespie run until end time (or optional: when bead unbinds) ###
        #####################################################################
        while t < t_max:

            #################################################
            ### Calculate cargo position and motor forces ###
            #################################################
            if dimension == '1D':
                try:
                    gsf.calc_force_1D(motor_team=motor_team, motor_fixed=motor_fixed, k_t=k_t, f_ex=f_ex, t=t, i=i, end_run=end_run, x_motor0=x_motor0)
                except em.NetForceError as err:
                    print("Force calculation error:", str(err))
                if end_run is True:
                    end_run = False
            else:
                gsf.calc_force_2D(team=motor_team, motor_0=motor_fixed, k_t=k_t, rest_length=rest_length, radius=radius, i=i)

            ##########################################################
            ### Update stepping- and unbinding rates of all motors ###
            ##########################################################
            for motor in motor_team:
                motor.stepping_rate()
                if dimension == '1D':
                    motor.unbind_rate_1D()
                else:
                    motor.unbind_rate_2D(dp_v1=dp_v1, dp_v2=dp_v2, T=temp, i=i)

            # Create lists of rates and corresponding ID's of motors for bookkeeping
            list_rates = []
            list_ids = []
            # Append rates of only the possible events (+ corresponding ID's)
            for motor in motor_team:
                id = motor.id
                if motor.unbound:
                    list_rates.append(motor.binding_rate)
                    list_ids.append(id)
                else:
                    list_rates.append(motor.epsilon)
                    list_ids.append(id)
                    list_rates.append(motor.alfa)
                    list_ids.append(id)
            # Sum of rates
            sum_rates = sum(list_rates)

            #############################
            ### Draw waiting time Tau ###
            #############################
            if sum_rates > 0:
                # draw tau, i.e. draw time before the next event will happen based on the sum of rates
                tau = np.random.exponential(scale=(1/sum_rates))
                # Update time
                motor_fixed.time_points[i].append(motor_fixed.time_points[i][-1] + tau)
                t += tau
                if t != motor_fixed.time_points[i][-1]:
                    raise AssertionError('Something wrong with time keeping')
            else:
                raise AssertionError(f'sum of rates is 0, time is {motor_fixed.time_points[i][-1]}, ln(time): {len(motor_fixed.time_points[i])},t={t}')

            ##############################
            ### Draw and execute event ###
            ##############################
            # Draw event to happen
            index = gsf.draw_event(list_rates=list_rates, sum_rates=sum_rates)
            # For control
            id_match = 0
            event_match = 0
            # Look up which event corresponds to 'index'
            for motor in motor_team:
                # Which motor
                if list_ids[index] == motor.id:
                    id_match += 1
                    # Which event
                    if list_rates[index] == motor.epsilon:
                        event_match += 1
                        # Save unbinding data
                        if dimension == '1D':
                            motor.forces_unbind.append(motor.f_current) # 1D
                            motor.run_length.append(motor.xm_rel)
                            if motor.direction == 'anterograde':
                                motor_fixed.antero_unbind_events[i] += 1
                            elif motor.direction == 'retrograde':
                                motor_fixed.retro_unbind_events[i] += 1
                            else:
                                raise AssertionError('Something wrong with motor direction')
                        elif dimension == '2D':
                            motor.fx_unbind.append(motor.f_x[i][-1]) # 2D
                            motor.fz_unbind.append(motor.f_z[i][-1]) # 2D
                            motor.run_length.append(motor.xm_rel)
                        else:
                            raise AssertionError('Something wrong with dimension settings')

                        ### INITIATE UNBINDING EVENT ###
                        motor.unbinding_event()

                    elif list_rates[index] == motor.alfa:
                        event_match += 1

                        ### INITIATE STEPPING EVENT ###
                        motor.stepping_event()

                    elif list_rates[index] == motor.binding_rate:
                        event_match += 1

                        ### INITIATE BINDING EVENT ###
                        if f_ex != 0:
                            # Make sure the cargo starts at X0
                            if motor_fixed.antero_bound[i][-1] + motor_fixed.retro_bound[i][-1] == 0:
                                cargo_distance = f_ex/motor.k_m
                                motor_bind = 0 - cargo_distance
                                motor.binding_event(x_cargo=motor_bind)
                            else:
                                motor.binding_event(x_cargo=motor_fixed.x_cargo[i][-1])
                        else:
                            motor.binding_event(x_cargo=motor_fixed.x_cargo[i][-1])
                    else:
                        raise AssertionError('Event index not found')
                else:
                    pass

            # There should be only one id match and one event match
            if id_match != 1:
                raise AssertionError(f'Motor protein ID should have 1 index match not {id_match}')
            if event_match != 1:
                raise AssertionError(f'Event should have 1 index match not {event_match}')

            # Save all current motor locations (also for the motors that didn't move: important for data analysis_plotting)
            for motor in motor_team:
                motor.x_m_abs[i].append(motor.xm_abs)

            #######################################
            ### Check how many motors are bound ###
            #######################################
            antero_bound = 0
            retro_bound = 0
            for motor in motor_team:
                if not motor.unbound:
                    if motor.direction == 'anterograde':
                        antero_bound += 1
                    elif motor.direction == 'retrograde':
                        retro_bound += 1
                    else:
                        raise AssertionError('Something wrong with motor direction')

            motor_fixed.antero_bound[i].append(antero_bound)
            motor_fixed.retro_bound[i].append(retro_bound)

            # If there are currently zero anterograde- and retrograde motors bound, the cargo has detached from the microtubule
            if antero_bound + retro_bound == 0:

                ## Save cargo run data ##
                motor_fixed.runlength_cargo[i].append(motor_fixed.x_cargo[i][-1]) # Save current cargo location = run length
                motor_fixed.time_unbind[i].append(t)  # Save current time: time of cargo unbind event
                motor_fixed.stall_time.append(tau) # Stall time: time before last motor unbinds
                end_run = True

                ## Stop the current Gillespie iteration i if single_run parameter is True, otherwise cargo will (re)attach next time step ##
                if single_run is True:
                    break

    return motor_team, motor_fixed

