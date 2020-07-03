// extrapolate fields to next-time step to obtain a better initial condition
p_rgh = (1.5 * p_rgh.oldTime() - 0.5* p_rgh.oldTime().oldTime());
// p_rgh.correctBoundaryConditions(); // will not work for fixedFluxPressure BC

// Update U and phi with their respective extrapolated value
phi = (1.5 * phi.oldTime() - 0.5* phi.oldTime().oldTime());

U = (1.5 * U.oldTime() - 0.5* U.oldTime().oldTime());
// U.correctBoundaryConditions(); updated at UEqn matrix formation

// MRF.correctBoundaryVelocity(U);

// -- only velocity and pressure values are extrapolated, since they are 
//    non-linearly coupled and benefits most from a better estimate.