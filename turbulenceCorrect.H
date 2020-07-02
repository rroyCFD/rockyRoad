// update the laminar and turbulent transport models.
laminarTransport.correct();
turbulence->correct();

// update modeled turbulent thermal conductivy now in case the the turbulence
// model does not update it (the OpenFOAM standard ones do not, but the SOWFA
// custom ones mostly do).

// ** kappat has calculated BC for correcponding fixedValue BC of nut, 
//    and can be estimated in postProcess: kappat = turbulence->nut()/Prt;

// update the boundary conditions on turbulent thermal conductivity, surface
// heat flux, and surface shear stress.
// kappat.correctBoundaryConditions();
qwall.correctBoundaryConditions();
Rwall.correctBoundaryConditions();
