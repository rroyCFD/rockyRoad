{
    // The inverse of the momentum equation "A" operator (the matrix diagonal) and
    // its projection from cell centers to cell faces.  Because it is a diagonal matrix,
    // The inverse is just the diagonal matrix of the reciprocals, hence the "r" in the
    // name "rAU".
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    // update mass-flux
    phi = fvc::flux(U); 

    // RhieChow flux correction for pressure gradient
    if(RhieChow)
    {
        tmp<surfaceVectorField> gradpDiff_f
        (
            (  (fvc::interpolate(fvc::grad(p_rgh)) & ed)
              - faceGradient.snGrad(p_rgh)
            ) * ed
        );

        phi += rAUf* (gradpDiff_f() & mesh.Sf());

        gradpDiff_f.clear();
    }

    // RhieChow flux correction for time-steps
    if(ddtCorr)
    {
        phi += rAUf*fvc::ddtCorr(U, phi);
    }

    // This part balances global mass flux.  It does it in a temporary field, and then
    // applies the correction indirectly by setting the pressure gradient to be used in
    // the fixedFluxPressure boundary condition on p_rgh, or directly if the zeroGradient
    // boundary condition on p_rgh is used.
    surfaceScalarField phiFixedFlux = phi;

    adjustPhi(phiFixedFlux, U, p_rgh);

    // update flux at zero-gradient BCs
    forAll(p_rgh.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p_rgh.boundaryField()[patchi]))
        {
            phi.boundaryFieldRef()[patchi] = phiFixedFlux.boundaryField()[patchi];
        }
    }

    // Update the fixedFluxPressure BCs to ensure flux consistency
    if(fixedFluxBC)
    {
        setSnGrad<fixedFluxPressureFvPatchScalarField>
        (
            p_rgh.boundaryFieldRef(),
            (
                phi.boundaryField()
              - phiFixedFlux.boundaryField()
            )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
        );
        p_rgh.correctBoundaryConditions();

        // Apply the same gradient to pp following the coresponding p_rgh BCs
        forAll(p_rgh.boundaryField(), patchi)
        {
            if(isType<fixedFluxPressureFvPatchScalarField> (p_rgh.boundaryField()[patchi]))
            {
                fixedFluxPressureFvPatchScalarField& p_rghGradPatch = refCast<fixedFluxPressureFvPatchScalarField>(p_rgh.boundaryFieldRef()[patchi]);
                fixedGradientFvPatchScalarField& ppGradPatch        = refCast<fixedGradientFvPatchScalarField>(pp.boundaryFieldRef()[patchi]);

                ppGradPatch.gradient() = p_rghGradPatch.gradient();
                Info << "maxGrad: " << gMax(p_rghGradPatch.gradient()) << endl;
            }
        }   
    }

    // Non-orthogonal corrector loop.
    while (spaece.correctNonOrthogonal())
    {
        fvScalarMatrix ppEqn
        (
            fvm::laplacian(rAUf, pp) == fvc::div(phi)
        );

        // For periodic cases at adjustPressureLevel header
        // p[pRefCell] = p_rgh[pRefCell] = pRefValue is ensured
        // thus, no correction required at the pRefCell
        ppEqn.setReference(pRefCell, 0.0);

        ppEqn.solve(mesh.solver(pp.select(spaece.finalInnerIter())));

        if (spaece.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi -= ppEqn.flux();

            // Explicitly relax pressure for momentum corrector
            pp.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rAU*fvc::grad(pp);
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            p_rgh += pp;
            p_rgh.correctBoundaryConditions();    
        }
    }

    #include "adjustPressureLevel.H"
}
