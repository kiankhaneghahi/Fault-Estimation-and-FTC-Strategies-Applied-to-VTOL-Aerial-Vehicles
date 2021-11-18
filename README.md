# Fault Estimation and Fault Tolerant Control Strategies For VTOL Vehicles In The Case of Soft Actuator Faults

This is the final project of my fault detection and control university course. The simulation was mainly based on the works done in the following journal.

G. Ortiz-Torres et al., "Fault Estimation and Fault Tolerant Control Strategies Applied to VTOL Aerial Vehicles With Soft and Aggressive Actuator Faults," in IEEE Access, vol. 8, pp. 10649-10661, 2020,

This project is used to estimate, isolate and diagnose faults for a quadcopter and a PVTOL and also use a methods to control the system by tolerating the fault. Both quadcopter and PVTOL systems have nonlinear dynamics. The ways for fault estimation in this project consist of nonlinear AO and linear PIO for the PVTOL and qLPV PIO for the quadcopter. The nominal controller in both systems uses unit quaternions. For soft fault in both vehicles, a fault accommodation method is implemented where the estimated fault is added to the nominal control signal to cancel the addetive fault. The aggressive fault was not considered because of problems in the simulation. Because the utilized dynamics in this project differ from those of the original journal which was not completely given and incomplete expression of the fault estimation and fault control methods, the simulation results were undesirable. Overall, the best result was for the linear PIO method.

This project was simulated in MATLAB and Simulink. The main files are the simulink files but the Quadcopter_Fault_Parameters_and_Flags.m file should be run first as a dependancy for the simulink that sets the parameters and etc. Each simulink file belongs to one of the two vehicles. More info is given in the presentation pdf file (in English) and in the report file (in persian).
