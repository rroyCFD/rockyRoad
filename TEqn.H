{
   
    // kappat = turbulence->nut()/Prt;
    // kappat.correctBoundaryConditions();

    // kappat have calculated BC when nut has fixedValue type BC
    volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + turbulence->nut()/Prt);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
      - fvc::div(qwall)
      - temperatureSourceTerm.force()
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    Boussinesq.updateDensityField();
}
