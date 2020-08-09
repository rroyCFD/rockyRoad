    // Solve the momentum equation
    
    // Update Coriolis forcing
    Coriolis.update();

    // Update buoyancy term
    Boussinesq.updateBuoyancyTerm();

    tmp<fvVectorMatrix> tUEqn
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
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve (
        UEqn 
      == 
        - fvc::grad(p_rgh)
        + fvc::reconstruct (Boussinesq.buoyancyTerm() * mesh.magSf())
    );
    
    fvOptions.correct(U);
