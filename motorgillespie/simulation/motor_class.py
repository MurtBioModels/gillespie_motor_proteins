import numpy as np


class MotorProtein(object):
    """ Creates motor protein instance(s) to parse through 'simulation.gillespie_simulation.gillespie_2D_walk()'.

    This class is specifically designed to use in combination with 'simulation.gillespie_simulation.gillespie_2D_walk()'.
    As the Gillespie simulation calls the correct methods at the correct time, using this class outside this project is not
    recommended as this can lead to serious errors. Conveniently, all this is handled by 'simulation.variable_loops.init_run()'.

    Attributes
    ----------
    family : string
            - For the Kinesin superfamily there are 14 families
            (Lawrence CJ, Dawe RK, Christie KR, et al. A standardized kinesin nomenclature.
            J Cell Biol. 2004;167(1):19-22. doi:10.1083/jcb.200408113)
    member : string
            For example KIF5A, a member of the family Kinesin-1.
    step_size : float or integer
                Step size [nm]
    k_m : float or integer
         Motor stiffness [pN/nm]
    alfa_0 : float or integer
           Zero-force stepping rate [1/s]
    f_s = 7 : float or integer
            Stall force [pN]
    epsilon_0 : tuple of floats or integers or integer or float
               Zero-force unbinding rate [1/s]
    f_d : float or integer
          Detachment force [pN]
    bind_rate : float or integer
               Binding rate [1/s]
    direction : string, 'retrograde' or 'anterograde'
                Stepping direction
    init_state : string, 'bound' or 'unbound'
                Determines if the motor starts bound or unbound,
                parsed through motor_init().
    calc_eps : string, 'gaussian', 'exponential' or 'constant'
               Function to calculate the unbinding rate in unbinding_rate_1D().
    id : int
            Numerical id for bookkeeping during simulation

    unbound : bool
             Keeps track every time step if the motor currently bound (False) or unbound (True)
    xm_abs : float
            Current location of motor, updated by stepping_event(). After calling unbinding_event() this attribute
            is set to NaN, and after calling binding_event() the attribute is set to the current cargo location.
            When initiating a 'MotorProtein' instance, this value is None until motor_init() is called.
    xm_rel : int
            Current location of motor, updated by stepping_event(). After calling unbinding_event() this value is
            added to 'self.run_length' and the attribute is set to NaN. After calling binding_event() the attribute is set to 0.
            When initiating a 'MotorProtein' instance, this value is None until motor_init() is called.
    f_current : float
                Current motor force [pN], updated every time step.
                When initiating a 'MotorProtein' instance, this value is None until calc_force_1D() or calc_force_2D()
                is called for the first time. For the 2D simulation this is the f_x.
    alfa : float
                Current stepping rate [1/s], updated every time step based on f_current.
                When initiating a 'MotorProtein' instance, this value is None until stepping_rate()
                is called for the first time.
    epsilon : float
                Current unbinding rate [1/s], updated every time step based on f_current.
                When initiating a 'MotorProtein' instance, this value is None until unbinding_rate_1D()
                or unbinding_rate_2D() is called for the first time.

    x_m_abs : list of list of floats
            Nested list (one list per Gillespie run) of motor locations per time step.
    forces : list of list of floats
            Nested list (one list per 1D Gillespie run) of motor forces per time step.
    f_x = : list of list of floats
            Nested list (one list per 2D Gillespie run) of x-component of force per time step.
    f_z : list of list of floats
            Nested list (one list per 2D Gillespie run) of z-component of force per time step.
    run_length : list of floats
                Collects the motor runlength values for the whole simulation (all Gillespie runs).
                Before unbinding, the current value of 'xm_rel' is appended to this list.
    forces_unbind : list of floats
                Collects the motor unbinding forces for the whole 1D simulation (all Gillespie runs).
                Before unbinding, the current value of 'f_current' is appended to this list.
    fx_unbind : list of floats
                Collects the x-component of the motor unbinding forces for the whole 2D simulation (all Gillespie runs).
                Before unbinding, the current value of 'f_x' is appended to this list.
    fz_unbind : list of floats
                Collects the z-component of the motor unbinding forces for the whole 2D simulation (all Gillespie runs).
                Before unbinding, the current value of 'f_z' is appended to this list.

    """

    ### Motor protein attributes ###
    def __init__(self, family, member, k_m, alfa_0, f_s, epsilon_0, f_d, bind_rate, step_size, direction, init_state, calc_eps, id):
        """ Constructor for the MotorProtein class.

        Parameters
        ----------
        family : string
                - For the Kinesin superfamily there are 14 families
                (Lawrence CJ, Dawe RK, Christie KR, et al. A standardized kinesin nomenclature.
                J Cell Biol. 2004;167(1):19-22. doi:10.1083/jcb.200408113)
        member : string
                For example KIF5A, a member of the family Kinesin-1.
        step_size : float or integer
                    Step size [nm]
        k_m : float or integer
             Motor stiffness [pN/nm]
        alfa_0 : float or integer
               Zero-force stepping rate [1/s]
        f_s : float or integer
                Stall force [pN]
        epsilon_0 : tuple of floats or integers or integer or float
                   Zero-force unbinding rate [1/s]
        f_d : float or integer
              Detachment force [pN]
        bind_rate : float or integer
                   Binding rate [1/s]
        direction : string
                    Stepping direction, 'retrograde' or 'anterograde'
        init_state : string
                    This attribute determines if the motor starts bound or unbound,
                    parsed through motor_init(), 'bound' or 'unbound'
        calc_eps : string
                   Function to calculate the unbinding rate in unbinding_rate_1D().
                   'exponential', 'gaussian' or 'constant'.
        id : int
                Numerical id for bookkeeping during simulation
        """
        ## Motor protein parameters ##
        self.family = family
        self.member = member
        self.k_m = k_m
        self.alfa_0 = alfa_0
        self.f_s = f_s
        self.epsilon_0 = epsilon_0
        self.f_d = f_d
        self.binding_rate = bind_rate
        self.step_size = step_size
        self.direction = direction
        self.init_state = init_state
        self.calc_eps = calc_eps
        self.id = int(id) + 1
        ## Updating variables ##
        self.unbound = None
        self.xm_abs = None
        self.xm_rel = None
        self.f_current = None
        self.alfa = None
        self.epsilon = None
        ## Motor data ##
        self.x_m_abs = []
        self.forces = []
        self.f_x = []
        self.f_z = []
        self.run_length = []
        self.forces_unbind = []
        self.fx_unbind = []
        self.fz_unbind = []

    ### Print info ###
    def __str__(self):
        """
        Returns
        -------
        str
            Short description of motor.
        """
        return f'Motor number {self.id}: {self.member} of {self.family}, steps {self.direction}'

    ### Check parameter values are correct ###
    def valid_motor(self, dimension):
        """
        Called once within by 'simulation.gillespie_simulation.gillespie_2D_walk()'
        before starting the simulation to check if all motors have valid attributes.

        Parameters
        ----------
        dimension: str, 1D or 2D
                 One of the Gillespie setting, defined by the user and parsed through
                 gillespie_2D_walk().
        """
        # Check if the force dependent unbinding relation
        calc_epsilon_options = ['gaussian', 'exponential', 'constant']
        if self.calc_eps not in calc_epsilon_options:
            raise ValueError("Not a valid motor: invalid force-unbinding relation. Expected one of: %s. Equations can be added in MotorProtein code" % calc_epsilon_options)
        # Check initial state
        initial_options = ['bound', 'unbound']
        if self.init_state not in initial_options:
            raise ValueError("Not a valid motor: invalid initial state. Expected one of: %s." % initial_options)
        # Check motor direction
        direction_options = ['anterograde', 'retrograde']
        if self.direction not in direction_options:
            raise ValueError("Not a valid motor: invalid polarity. Expected one of: %s." % direction_options)
        # Check motor parameters
        if self.f_d == 0:
            raise ValueError("Invalid detachment force (Fd) force. Deliminator can not be zero in equations")
        if self.f_s == 0:
            raise ValueError("Invalid stall force (Fs). Deliminator can not be zero in equations")
        # Check if dimensional restrictions are met
        if dimension == '1D':
            if type(self.epsilon_0) != int and type(self.epsilon_0) != float:
                raise ValueError(f"Not a valid motor: zero force unbinding rate should be an integer or float, not {type(self.epsilon_0)}")
        if dimension == '2D':
            if type(self.epsilon_0) != list:
                raise ValueError(f"Not a valid motor: zero force unbinding rate should be a list, not {type(self.epsilon_0)}")
            if len(self.epsilon_0) != 2:
                raise ValueError(f"Not a valid motor: zero force unbinding rate list should contain two values, not {len(self.epsilon_0)}")
        return

    ### Initiate motor to default state every Gillespie run (t=0) ###
    def motor_init(self, dimension, f_ex):
        """
        Called before every Gillespie run to set
        initial 'unbound', 'xm_abs', 'xm_rel' and 'x_m_abs' based
        on 'dimension' and 'f_ex'.

        Parameters
        ----------
        dimension: str
                 Determines if the simulation is one dimensional or two dimensional (angle included).
                 This is one of the Gillespie setting, defined by the user and parsed through
                 gillespie_2D_walk(). '1D' or '2D'.
        f_ex: float
                 Constant external force [pN]. This is one of the simulation parameters,
                 defined by the user and parsed through gillespie_2D_walk().
        """
        if self.init_state == 'unbound':
            self.unbound = True
            self.xm_abs = float('nan')
            self.xm_rel = float('nan')
            self.x_m_abs.append([float('nan')])
        elif self.init_state == 'bound':
            self.unbound = False
            self.xm_abs = 0
            self.xm_rel = 0
            if f_ex != 0:
                self.x_m_abs.append([])
            elif f_ex == 0:
                self.x_m_abs.append([0])
            else:
                raise ValueError("Something wrong with f_ex")
        else:
            raise ValueError("Motor proteins can either be in an bound or unbound initial state")

        if dimension == '1D':
            self.forces.append([])
        elif dimension == '2D':
            self.f_x.append([])
            self.f_z.append([])
        else:
            raise ValueError(f"The simulation settings except only 1D or 2D, not {dimension}")

        return

    ### Force dependent rate equations ###
    def stepping_rate(self):
        """
        Calculates the stepping rate based on the current motor force
        'f_current'. Since the force is calculated by k_m*(x_m - x_c), the motors (and thus cargo)
        start at x=0, and calling stepping_event() on retrograde motors subtract its 'step_size' from
        the current location 'xm_abs', this would result in retrograde motors keeping 'alfa' = 'alfa_0'.
        To correct this, a condition is added to correct for the negative values of retrograde motors.

        See Also
        --------
        simulation.gs_functions.calc_force_1D
        """
        f_current = self.f_current
        direction = self.direction
        f_s = self.f_s
        if direction == 'anterograde':
            if f_current < 0:
                self.alfa = self.alfa_0
            elif f_current > f_s:
                self.alfa = 0
            else:
                self.alfa = self.alfa_0 * (1 - (f_current / f_s))
        elif direction == 'retrograde':
            if f_current > 0:
                self.alfa = self.alfa_0
            elif (-1*f_current) > f_s:
                self.alfa = 0
            else:
                self.alfa = self.alfa_0 * (1 - ((-1 * f_current) / f_s))
        else:
            raise ValueError("Transport must be either retrograde or anterograde")

        return

    def unbind_rate_1D(self):
        """
        Calculates the unbinding rate based on the current motor force
        'f_current'. Based on 'calc_eps', either a gaussian, exponential or
        constant force-unbinding relation is used. This methods is called when
        dimension = '1D' within gillespie_2D_walk().

        """

        f_current = (self.f_current**2)**0.5
        calc_eps = self.calc_eps
        #
        if calc_eps == 'exponential':
            self.epsilon = self.epsilon_0 * (np.e**(f_current / self.f_d))
        elif calc_eps == 'gaussian':
            self.epsilon = 1 + (6 * np.exp(-(f_current - 2.5) ** 2))
        elif calc_eps == 'constant':
            self.epsilon = self.epsilon_0
        else:
            raise ValueError("Unbinding equation not recognized. Retry or add equation to unbinding_rate_1")

        return

    def unbind_rate_2D(self, dp_v1, dp_v2, T, i):
        """
        Calculates the unbinding rate based on the current motor force
        components 'f_x' and 'f_z'. This methods is called when
        dimension = '2D' within gillespie_2D_walk().

        Parameters
        ----------
        dp_v1 : np.array[float]
            Displacement vector strongly bound state [d_1_x, d_1_z], np.array.
        dp_v2 : np.array[float]
            Displacement vector weakly bound state [d_2_x, d_2_z], np.array.
        T : float or integer
            Temperature in Kelvin [K]
        i : int
            current iteration

        Notes
        -------
        - Khataee H, Mahamdeh M, Neufeld Z. Processivity of molecular motors under vectorial loads.
        Phys Rev E. 2020 Aug 18;102(2):022406.
        - Khataee H, Howard J. Force Generated by Two Kinesin Motors Depends on the Load Direction and
        Intermolecular Coupling. Phys Rev Lett. 2019 May 8;122(18):188101.

        """
        if T == 0:
            raise ValueError("Invalid temperature (T). Deliminator can not be zero")
        # Check if 2d conditions are given
        if len(self.epsilon_0) != 2:
            raise ValueError("Zero force unbinding rate list should contain two values")
        f_v = np.array([-1*self.f_x[i][-1], self.f_z[i][-1]])
        Boltzmann = 1.38064852e-23
        k1 = self.epsilon_0[0] * np.exp((np.dot(f_v, dp_v1)) / (Boltzmann * T))
        k2 = self.epsilon_0[1] * np.exp((np.dot(f_v, dp_v2)) / (Boltzmann * T))

        self.epsilon = k1 * k2 / (k1 + k2)

        return

    ### Events ###
    def stepping_event(self):
        """
        When called, the motor location 'xm_abs' and traveled distance 'xm_rel'
        are increased by 'step_size' for anterograde motors, and decreased by 'step_size' for
        retrograde motors. Hence, this method mimics motor stepping within the simulation.
        """
        step_size = self.step_size
        direction = self.direction
        if direction == 'anterograde':
            self.xm_abs = self.xm_abs + step_size
            self.xm_rel = self.xm_rel + step_size
        elif direction == 'retrograde':
            self.xm_abs = self.xm_abs - step_size
            self.xm_rel = self.xm_rel - step_size
        else:
            raise ValueError("Transport must be either retrograde or anterograde")

        return

    def binding_event(self, x_cargo):
        """
        When called, the motor location 'xm_abs' is set to the cargo location x_cargo
        and traveled distance 'xm_rel' is set to 0. Also, the 'unbound' attributed is set to False.
        Within 'gillespie_2D_walk()', this methods can only be called when the motor is currently unbound
        (i.e. 'unbound' = True). Hence, this method mimics motor (re)binding within the simulation.

        Parameters
        ----------
        x_cargo : float
                Current cargo position. Within  gillespie_2D_walk(), this is retrieved from the (single)
                'MotorFixed' instance.
        """

        self.xm_abs = x_cargo
        self.xm_rel = 0
        self.unbound = False

        return

    def unbinding_event(self):
        """
        When called, the motor location 'xm_abs' and traveled distance 'xm_rel' are both
        set to NaN. Also, the 'unbound' attributed is set to True. Hence, this method mimics motor
        (re)binding within the simulation. It is important that this value is set to NaN instead of 0
        for later data analysis_plotting, as this distinguishes when the motor was bound at x = 0 or not bound at all.

        """

        self.xm_abs = float('nan')
        self.xm_rel = float('nan')
        self.unbound = True

        return


class MotorFixed(object):
    """ Creates a fixed motor instance to parse through 'simulation.gillespie_simulation.gillespie_2D_walk()'.

    A single instance of this class is used within 'simulation.gillespie_simulation.gillespie_2D_walk()' to hold
    simulation parameters and to collect cargo data. When a trap stiffness is included, this instance mimics the
    optical trap by being fixed at x = 0 (hence the name) and the 'k_t' attribute is the trap stiffness  when > 0
    (see 'gs_functions' module for more information on calculation of cargo location and motor forces).

    This class is specifically designed to use in combination with 'simulation.gillespie_simulation.gillespie_2D_walk()'.
    As the Gillespie simulation calls the correct methods at the correct time, using this class outside this project is not
    recommended as this can lead to serious errors. Conveniently, all this is handled by 'simulation.variable_loops.init_run()'.

    Attributes
    ----------
    k_t : int or float
        Trap stiffness.
    f_ex: int or float
        Constant external force.

    Only for 2D simulation:

    dp_v1 : np.array[float]
        Displacement vector strongly bound state [d_1_x, d_1_z], np.array.
    dp_v2 : np.array[float]
        Displacement vector weakly bound state [d_2_x, d_2_z], np.array.
    radius : int or float
        Radius of the cargo.
    rest_length : int or float
        Rest length of the motor (2D simulations are only possible for 1 motor in the current version).
    temp : int or float
            Temperature in Kelvin [K].

    Calculate later:

    angle : None
        Angle between the horizontal stepping plane and cargo linker.
    init_antero : None
        How many anterograde motors start bound every Gillespie run,
        calculated once based on the motor instances.
    init_retro :  None
        How many retrograde motors start bound every Gillespie run,
        calculated once based on the motor instances.

    Cargo/simulation data:

    time_points : list of lists
        Time points of the events which are drawn by the Tau-leaping Gillespie algorithm.
    x_cargo : list of lists
        Cargo locations per time point.
    antero_bound : list of lists
        Number of bound anterograde motors per time point.
    retro_bound : list of lists
        Number of bound retrograde motors per time point.
    antero_unbind_events : list
        Unbind events of anterograde motors per Gillespie run.
    retro_unbind_events : list
        Unbind events of retrograde motors per Gillespie run.
    runlength_cargo : list of lists
        Run lengths per Gillespie runs (amount of run lengths per Gillespie run differs).
        Divide by k_t to get force (of optical trap)
    time_unbind : list of lists
        Time of unbinding. To calculate the run time, the difference between the entries have
        to be taken (np.diff() for example)
    stall_time : list
        Time between the the last motor being still bound and the cargo unbinding.
    """

    ## Class attributes ##
    __id = 0
    __unbound = False
    __x_m = 0

    def __init__(self, k_t, f_ex, dp_v1=None, dp_v2=None, radius=None, rest_length=None, temp=None):
        """ Constructor for the MotorFixed class.

        Parameters
        ----------
        k_t : int or float
            Trap stiffness.
        f_ex: int or float
            Constant external force.
        dp_v1 : np.array[float]
            Displacement vector strongly bound state [d_1_x, d_1_z], np.array.
            Optional, only for 2D simulations, default = None.
        dp_v2 : np.array[float]
            Displacement vector weakly bound state [d_2_x, d_2_z], np.array.
            Optional, only for 2D simulations, default = None.
        radius : int or float
            Radius of the cargo.
            Optional, only for 2D simulations, default = None.
        rest_length : int or float
            Rest length of the motor (2D simulations are only possible for 1 motor in the current version).
            Optional, only for 2D simulations, default = None.
        temp : int or float
                Temperature in Kelvin [K].
            Optional, only for 2D simulations, default = None.


        """
        ## Simulation parameters ##
        self.k_t = k_t
        self.f_ex = f_ex
        self.dp_v1 = dp_v1
        self.dp_v2 = dp_v2
        self.radius = radius
        self.rest_length = rest_length
        self.temp = temp
        # Calculate later
        self.angle = None
        self.init_antero = None
        self.init_retro = None
        ## Cargo/simulation data  ##
        self.time_points = []
        self.x_cargo = []
        self.antero_bound = []
        self.retro_bound = []
        self.antero_unbind_events = []
        self.retro_unbind_events = []
        self.runlength_cargo = []
        self.time_unbind = []
        self.stall_time = []

    ### Initiate ###

    def init_valid_once(self, motor_team):
        """
        Called once within by 'simulation.gillespie_simulation.gillespie_2D_walk()'
        before starting the simulation to check if the attributes are valid,
        and initiates the simulation based on the validated attributes/simulation parameters.

        Parameters
        ----------
        motor_team : list of motorgillespie.simulation.motor_class.MotorProtein
            Team of motor protein objects.
            Collects per motor data like locations and forces over time.

        """
        antero_init_bound = 0
        retro_init_bound = 0
        for motor in motor_team:
            if motor.init_state == 'bound':
                if motor.direction == 'anterograde':
                    antero_init_bound += 1
                if motor.direction == 'retrograde':
                    retro_init_bound += 1
        if self.f_ex != 0 and (antero_init_bound + retro_init_bound) == 0:
            raise ValueError(f'If an external force is added to the simulation, at least one motor has to start bound')
        self.init_antero = antero_init_bound
        self.init_retro = retro_init_bound

        if self.k_t != 0 and self.f_ex != 0:
            raise ValueError(f'Either add an external force or a trap stiffness to the simulation, not both.')

    def fixed_init(self):
        """
        Called before every Gillespie run to correctly append
        to the data collection attributes.

        """
        self.time_points.append([0])
        self.antero_bound.append([self.init_antero])
        self.retro_bound.append([self.init_retro])
        self.retro_unbind_events.append(0)
        self.antero_unbind_events.append(0)
        self.x_cargo.append([])
        self.runlength_cargo.append([])
        self.time_unbind.append([])

        return

    ### Set angle in 2D simulation ###
    def calc_angle(self):
        """
        Called once within by 'simulation.gillespie_simulation.gillespie_2D_walk()'
        before starting the simulation to calculate the angle between the horizontal
        stepping plane and cargo linker.

        """
        self.angle = np.arcsin(self.radius / (self.radius + self.rest_length))

    ### Getters for private attributes ###
    @property
    def x_m(self):
        """Location of fixed motor = 0 """
        return self.__x_m

    @property
    def id(self):
        """id of fixed motor = 0 """
        return self.__id

    @property
    def unbound(self):
        """Fixed motor is always bound """
        return self.__unbound

    # Static methods useful for analysing and validating data
    @staticmethod
    def stepping_rate(alfa0, f_s, f_current, direction):
        """ Calculate the motor stepping rate without any instance.

        Calculates the stepping rate based on the current motor force
        'f_current'. Since the force is calculated by k_m*(x_m - x_c), the motors (and thus cargo)
        start at x=0, and calling stepping_event() on retrograde motors subtract its 'step_size' from
        the current location 'xm_abs', this would result in retrograde motors keeping 'alfa' = 'alfa_0'.
        To correct this, a condition is added to correct for the negative values of retrograde motors.

        See Also
        --------
        simulation.gs_functions.calc_force_1D
        """

        if direction == 'anterograde':
            if f_current < 0:
                alfa = alfa0
            elif f_current > f_s:
                alfa = 0
            else:
                alfa = alfa0 * (1 - (f_current / f_s))
        elif direction == 'retrograde':
            if f_current > 0:
                alfa = alfa0
            elif (-1*f_current) > f_s:
                alfa = 0
            else:
                alfa = alfa0 * (1 - ((-1 * f_current) / f_s))
        else:
            raise ValueError("Transport must be either retrograde or anterograde")

        return alfa

    @staticmethod
    def unbind_rate_1D(epsilon0, f_d, f_current, calc_eps):
        """ Calculate the 1D motor unbinding rate without any instance.

        Calculates the unbinding rate based on the current motor force
        'f_current'. Based on 'calc_eps', either a gaussian, exponential or
        constant force-unbinding relation is used. This methods is called when
        dimension = '1D' within gillespie_2D_walk().

        """

        f_current = (f_current**2)**0.5
        calc_eps = calc_eps
        #
        if calc_eps == 'gaussian':
            epsilon = 1 + (6 * np.exp(-(f_current - 2.5) ** 2))
        elif calc_eps == 'exponential':
            epsilon = epsilon0 * (np.e**(f_current / f_d))
        elif calc_eps == 'constant':
            epsilon = epsilon0
        else:
            raise ValueError("Unbinding equation not recognized. Retry or add equation to unbinding_rate_1")

        return epsilon
