
This Python project was made in 2022/2023 during my internship (major research project) at the department
of *Cellbiology, neurobiology and biophysics* at Utrecht University. I was investigating the biophysical 
properties of motor proteins and their effect on bidirectional cargo transport.

### Abstract
Intracellular cargo is actively transported by motor proteins stochastically stepping, binding and unbinding
along the microtubule network, often by multiple motors at the same time. These motor proteins have to
cooperate to allow for efficient long-distance transport, which is vital for maintaining proper cell organisation
and -functioning. Single-motor behaviour is heavily influenced by the force acting on the motor, and this has
been extensively studied and well quantified. However, understanding how single-motor mechanical properties
contribute to collective transport dynamics is still a topic with many open questions. This is especially
emphasised by the discovery that both kinesins and dyneins are simultaneously bound to some cargo, pulling
against each other in opposite directions, a process known in literature as ’tug-of-war’, or bidirectional
transport. As motors are elastically coupled to the cargo, they exert strain induced forces on each other,
affecting each others behaviour directly. Hence, the strength of the elastic coupling, represented by the motor
stiffness, affects how motors cooperate within a team. To investigate this, we developed a simulation tool that
allows us to simulate any number and combination of motors and parameters under defined external forces.
We use this to examine the effect of elastic coupling on collective cargo transport. Our results show that
increasing motor stiffness has a negative effect on cooperation and overall cargo transport. We then used this
knowledge to bias a tug-of-war between two completely symmetrical teams in one direction by only changing
the motor stiffness of one team. We found that changing the motor stiffness is indeed enough to mechanically
regulate the direction of transport. Taken together, our results emphasize the importance of elastic coupling
strength as a mechanical regulatory factor in collective- and bidirectional transport.
Research

### Python project
This Python project contains the package **motorgillespie** which contains the following modules:
- **gillespie_simulation.py** is the actual Gillespie simulation
- **gs_functions.py** contains supporting functions 
- **motor_class.py** contains the classes used by initiate_motors.py
- **initiate_motors.py** contains the functions to create motor objects to be used in gillespie_simulation.py.
- **variable_loops.py** contains a function that calls initiate_motors.py and gillepsie_simulation.py.

See **simulation_tool.png**

The scripts in the **end_report** folder are examples I used during my project that use 
the package.





