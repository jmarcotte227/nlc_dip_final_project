# nlc_dip_final_project

David Chen and John Marcotte's final project for Non-Linear Control: the control of a double inverted pendulum. Some important files relevant to the report are: 

1. passivity.m - Contains the symbolic equations for the derivation of the passivity-based control law.
2. unc_model.m - Contains the simulation and symbolic derivation of the error gradient used to update the model parameters in-process. Can be used to simulate the passivity control with known parameters by disabling the parameter update, and setting the estimates to the true values.
3. pole_placement.m - Contains the simulation for feedback linearization with decoupled dynamics and pole placement.
4. state_estimator.m - Contains the simulation for tracking the system to a reference with sensor bias and noisy measurement data.
5. show.m - Can be used to visualize the systems motion after running one of the simulations.
