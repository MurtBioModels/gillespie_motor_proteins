import motorgillespie.simulation.motor_class as mc

def init_motor_fixed(sim_params):
    """
    Initialize a MotorFixed instance. It is mainly used for collecting cargo data and to hold
    simulation parameters. The MotorFixed instance merely has a position, which is 0 always (hence the name MotorFixed).
    When a trap stiffness is added, the MotorFixed object serves as the optical trap, with x = 0
    and k_t being used to calculate the cargo location, and pulling the cargo back to x = 0 when the
    last motor unbinds. Here, the opposing trap force = k_t * x_cargo (see the 'calc_force_1D' function in thw 'gs_functions' module).

    Parameters
    ----------
    sim_params : dict
        A dictionary containing simulation parameters.

    Returns
    -------
    motor_fixed :  motorgillespie.simulation.motor_class.MotorFixed
        Instance of class 'MotorFixed'.
        Fixed location at x = 0.
        Holds simulation parameters and collects cargo data.

    Notes
    -----
      The function creates a MotorFixed object based on the provided simulation parameters.
      The simulation parameters should include at least the following keys:
        - 'k_t' : float
            Trap stiffness [pN/nm], k_t = 0 when excluded.
        - 'f_ex' : float
            Constant external force [pN], f_ex = 0 when excluded.
      and when initiating a two-dimensional simulation:
        - 'dp_v1' : float
            Displacement vector strongly bound state [d_1_x, d_1_z]
        - 'dp_v2' : float
            Displacement vector weakly bound state [d_2_x, d_2_z]
        - 'radius' : float
            Radius of the cargo [nm].
        - 'rest_length' : float
            Rest length of motor linkers [nm].
        - 'temp' : float
            Temperature in Kelvin [K]
      The MotorFixed class is defined in the 'motor_class' module, see it's documentation
      for more information.

    """
    # Create fixed motor
    if all(k in sim_params for k in ("dp_v1", "dp_v2", 'radius', 'rest_length', 'temp')):
        motor_fixed = mc.MotorFixed(sim_params['k_t'],  sim_params['f_ex'], sim_params['dp_v1'], sim_params['dp_v2'], sim_params['radius'], sim_params['rest_length'], sim_params['temp'])
    else:
        motor_fixed = mc.MotorFixed(sim_params['k_t'],  sim_params['f_ex'])
    return motor_fixed


def init_mixed_team(mnr, *motor_params):
    """
    Initialize a list of MotorProtein instance(s).

    Parameters
    ----------
    mnr : list
     A list containing the number of motors per species/kind (biological species or
     own set of parameters) to be simulated.
    motor_params : dict
     A dictionary per species/kind. The order in which the dictionaries are parsed
     has to be the same as the order of quantities in nmr.


    Returns
    -------
    motor_team : list [motorgillespie.simulation.motor_class.MotorProtein]
          List of MotorProtein instances representing the team of motor proteins.

    Notes
    -----
      The function creates a list of MotorProtein object(s) based on the provided simulation parameters.
      The MotorProtein class is defined in the 'motor_class' module, see it's documentation
      for more information on the keys that need to be included.

    Example
    -----
    Example: You want to simulate two Kinesin-1 and three Kinesin-2.
    init_mixed_team(mnr = [2,3], kinesin_1_params, kinesin_2_params)
    or
    init_mixed_team(mnr = [3,2], kinesin_2_params, kinesin_1_params)

    """

    motor_team = []
    for index, params in enumerate(motor_params):
        for i in range(mnr[index]):
            motor_team.append(mc.MotorProtein(params['family'], params['member'], params['k_m'], params['alfa_0'], params['f_s'], params['epsilon_0'], params['f_d'], params['bind_rate'], params['step_size'], params['direction'], params['init_state'], params['calc_eps'], len(motor_team)))
            print(f'{motor_team[-1].id} - {motor_team[-1].direction}')

    return motor_team


