    // Solve the momentum equation

    // Update Coriolis forcing
    Coriolis.update();

    // Update buoyancy term
    Boussinesq.updateBuoyancyTermVol();

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
      + fvc::div(Rwall)
      - momentumSourceTerm.force()
      - Coriolis.force()
     ==
        fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (spaece.momentumPredictor())
    {
        solve (UEqn == - fvc::grad(p_rgh) + Boussinesq.buoyancyTermVol() );
        fvOptions.correct(U);
    }
