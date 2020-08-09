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

    if(regOn)
    {
        // - Only applicable for A6
        //   Extrapolated or updated phi and U values are used for regularization
        solve(
            UEqn 
          == 
            - fvc::grad(p)
            + fvc::reconstruct (Boussinesq.buoyancyTerm() * mesh.magSf()) 
            + RegModel->convectionResidual(phi,U));
    }
    else
    {
        solve (
            UEqn 
          == 
            - fvc::grad(p_rgh)
            + fvc::reconstruct (Boussinesq.buoyancyTerm() * mesh.magSf())
        );
    }

    fvOptions.correct(U);
