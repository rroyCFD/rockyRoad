    // Solve the momentum equation
    
    // Update Coriolis forcing
    Coriolis.update();

    // Update buoyancy term
    // Boussinesq.updateBuoyancyTerm();
    Boussinesq.updateBuoyancyTermVol();

    fvVectorMatrix UEqn
    (
        fvm::ddt(U) 
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
      + fvc::div(Rwall)
      - momentumSourceTerm.force()
      - Coriolis.force()
      - Boussinesq.buoyancyTermVol()
     ==
        fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            - fvc::grad(p_rgh)
            // fvc::reconstruct
            // (
            //     (
            //       - fvc::snGrad(p_rgh)
            //       // + Boussinesq.buoyancyTerm()
            //     )*mesh.magSf()
            // )
        );

        fvOptions.correct(U);
    }
