# Extended SLIP model

mode1.m - leg dynamics + set angle actuation in stance phase
mode2.m - leg dynamics + leg retraction in flight phase
run.m   - sequence of mode2->mode1->mode2->mode1... simulation for testing
guard_stance.m - detect liftoff
guard_flight.m - detect touchdown

Flight-Dynamics_EulerLagrangeII.mlx - calculation of eq. of motion with E-Lagrange algorithm
Stance-Dynamics_EulerLagrangeII.mlx - calculation of eq. of motion with E-Lagrange algorithm