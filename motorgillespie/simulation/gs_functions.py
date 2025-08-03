import motorgillespie.errors as em
import numpy as np


def calc_force_1D(motor_team, motor_fixed, k_t, f_ex, t, i, end_run, x_motor0):
    """ Calculates the cargo position and motor forces in a one-dimensional system.

    This function first calculates the cargo position 'cargo_loc' based on the location
    ('MotorProtein.xm_abs') and stiffness ('MotorProtein.k_m') of all currently
    bound motors, and when applicable the external force ('f_ex') or trap stiffness ('k_t').
    After this, current motor forces are calculated for each individual motor based on their
    location, their motor stiffness and the cargo location.

    Parameters
    ----------
    motor_team : list [motorgillespie.simulation.motor_class.MotorProtein]
          List of MotorProtein instances representing the team of motor proteins.
    motor_fixed : motorgillespie.simulation.motor_class.MotorFixed
              Cargo object
    k_t : int or float
        Trap stiffness.
    f_ex: int or float
        Constant external force.
    i : int
        Current iteration for list indexing
    t : int
        Current time
    end_run : bool
             If last time step the cargo detached, end_run == True
    x_motor0 : mc.MotorFixed.__x_m
              Location of fixed motor, private MotorProtein attribute, == 0.

    """

    # Calculate the numerator: sum of km_i * xm_i of bound motors in the motor_team
    xm_km_sum = sum([(motor.xm_abs * motor.k_m) for motor in motor_team if not motor.unbound])
    # Subtract external force in the numerator
    if f_ex != 0:
        xm_km_sum += f_ex

    # Calculate the denominator: sum of km_i of bound motors in the motor_team
    km_sum = sum([motor.k_m for motor in motor_team if not motor.unbound])
    # Add trap stiffness in the denominator
    if k_t != 0:
        km_sum += k_t

    ##############################
    ## Calculate position cargo ##
    ##############################

    ## If no external force is present ##
    if f_ex == 0:
        net_force = 0
        # Start of Gillespie run
        if t == 0:
            # No motors are bound at t=0, denominator is 0, artificially start cargo at x=0
            if motor_fixed.antero_bound[i][-1] == 0 and motor_fixed.retro_bound[i][-1] == 0:
                cargo_loc = 0
            # At least one motor starts bound at t=0
            else:
                cargo_loc = xm_km_sum/km_sum
        # Cargo detached last time step, denominator is 0, artificially but cargo back at x=0
        elif end_run is True:
            cargo_loc = 0
        # t>0, no cargo detachment last time step, i.e. at least one motor is bound
        else:
            cargo_loc = xm_km_sum/km_sum

    ## If external force is present ##
    else:
        # Start of Gillespie run, there is always at least one motor bound at t=0 and f_ex != 0
        if t == 0:
            net_force = f_ex
            # Artificially start cargo at x=0
            cargo_loc = 0
            # Calculate initial distance from the cargo to bound motor(s)
            cargo_distance = xm_km_sum/km_sum
            for motor in motor_team:
                if motor.unbound is False:
                    xm = 0 - cargo_distance
                    # Change from 0 to correct position
                    motor.xm_abs = xm
                    # First entry of motor.x_m_abs[i] is the calculated motor position
                    motor.x_m_abs[i].append(xm)
        # Cargo detached last time step, denominator is 0, artificially but cargo back at x=0
        elif end_run is True:
            net_force = 0
            cargo_loc = 0
        # t>0, no cargo detachment last time step, i.e. at least one motor is bound
        else:
            net_force = f_ex
            cargo_loc = xm_km_sum/km_sum

    # Save and append current cargo location
    motor_fixed.x_cargo[i].append(cargo_loc)

    ############################
    ## Calculate motor forces ##
    ############################
    for motor in motor_team:
        if motor.unbound is True:
            f = float('nan')
        else:
            f = motor.k_m * (motor.xm_abs - cargo_loc)
            net_force += f

        motor.f_current = f
        motor.forces[i].append(f)

    # If trap stiffness is present, calculate trap force using MotorFixed location (this is 0),
    # and thus k_t * - cargo locations is added to net_force.
    if k_t != 0:
        f0 = k_t*(x_motor0 - cargo_loc)
        net_force += f0

    # net_force is used to check if forces balance out: sum(motor_force_i) == f_ex,
    # sum(motor_force_i) == k_t * cargo_loc or sum(motor_force_i) == 0.
    # Due to rounding errors, a margin is used.
    if (net_force**2)**0.5 > 10**-10:
        print(f'Net force = {net_force}')
        raise em.NetForceError('Net force on cargo should be zero. Check for problems in the code.')

    return


def calc_force_2D(team, motor_0, k_t, rest_length, radius, i):
    """ Calculates the cargo position and motor forces in a two-dimensional system.

    This function updates the current force acting on the motor protein (ONLY ONE MOTOR!)
    in the Gillespie Stochastic simulation gillespie_motor_team

    Parameters
    ----------
    team : list of MotorProtein
           Team of motor protein objects for simulation simulation
    motor_0 : MotorFixed
              Instance of 'MotorFixed' class replacing the optical trap.
              Fixed location at x = 0.
    k_t : float
    rest_length : float or integer
    radius : float or integer
    i : integer
        Current iteration for list indexing

    """

    # List of Km_i*Xm_i of motor proteins in motor motor_team
    if len(team) != 1:
        raise ValueError("The 2D simulation is currently only available for 1 motor simulations")
    xm_km_list = [(motor.xm_abs * motor.k_m) for motor in team]
    xm_km_list.append(motor_0.k_t * motor_0.xm_abs)
    # List of Km's of bound motor proteins in motor motor_team
    km_list = [motor.k_m for motor in team if not motor.init_state]
    km_list.append(k_t)
    # Calculate position cargo
    cargo_loc = sum(xm_km_list) / sum(km_list)
    motor_0.x_cargo[i].append(cargo_loc)
    if km_list == 0:
        raise ValueError("Invalid sum of Km's. Deliminator can not be zero")

    # Update forces acting on each individual motor protein
    net_force = 0
    for motor in team:
        if motor.unbound:
            motor.f_current = 0
            motor.f_x[i].append(0)
            motor.f_z[i].append(0)
            net_force += motor.f_current
        else:
            fx = motor.k_m * (motor.xm_abs - cargo_loc)
            motor.f_current = fx
            motor.f_x[i].append(fx)
            motor.f_z[i].append( fx / np.sqrt( ( (1 + (rest_length / radius) )**2 ) - 1)  ) # Dit nog checken met de simpele?
            net_force += motor.f_current

    f_fixed = motor_0.k_t*(motor_0.xm_abs - cargo_loc)
    net_force += f_fixed

    if abs(net_force) > 10**-10:
        print(net_force)
        raise AssertionError('Net force on cargo should be zero, look for problems in code')
    return


def draw_event(list_rates, sum_rates):
    """ Randomly selects an event to occur based on the given rates.

    This function determines which event will happen by drawing a random number from a uniform distribution and
    iterating through the list of rates 'list_rates'. Every iteration the conditional statement checks if the sum_prob
    is equal or smaller then the randomly chosen number. If this condition is true, the rate at index 'i' in 'list_rates'
    is divided by the sum of rates and added to 'sum_prob', and the index is raised by one. When the beforementioned condition
    is false, the event associated with the previous rate, for which the condition was true, is chosen for execution.

    This function is part of a Gillespie algorithm-based simulation, where 'list_rates' contains the rates
    of all possible events for all proteins. These rates are updated each iteration based on the
    current force acting on each motor and the bound or unbound state of each motor.

    Code based on: http://be150.caltech.edu/2019/handouts/12_stochastic_simulation_all_code.html

    Parameters
    ----------
    list_rates : list of floats
        List of rates of all possible events for all proteins.
        Example:
        - 1 bound + 1 unbound: [stepping rate motor1, unbinding rate motor1, attachment rate motor2]
        - Corresponding ID list: [motor1, motor1, motor2]

    sum_rates: float
        The total sum of rates in 'list_rates'.

    Returns
    -------
    int
        Index of the event to happen in 'list_rates'.
    """

    # Generate random number
    rand = np.random.rand() # >>> save random seed for reproducibility

    # Find event
    i = 0
    probs_sum = 0.0
    while probs_sum <= rand: # < vs <=
        probs_sum += list_rates[i]/sum_rates
        i += 1

    return i - 1


