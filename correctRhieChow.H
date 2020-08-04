// RhieChow flux correction for pressure gradient
if(RhieChow)
{
    // // this Rhie chow correction uses orthogonal gradient for snGrad and
    // // aligned gradient for face-interpolated-gradient
    // tmp<surfaceVectorField> gradpDiff_f
    // (
    //     (  (fvc::interpolate(fvc::grad(p_rgh)) & ed)
    //       - faceGradient.snGrad(p_rgh)
    //     ) * ed
    // );
    // phi += rAUf* (gradpDiff_f() & mesh.Sf());
    // gradpDiff_f.clear();

	// Regardless of magnitude of the face non-orthogonality,
	// gradients are estimated in the surface normal direction
    tmp<surfaceScalarField> gradpDiffFlux
    (
          (fvc::snGrad(p_rgh) * mesh.magSf())
        - (fvc::interpolate(fvc::grad(p_rgh)) & mesh.Sf())
    );

    // following the sign convection (snGrand - interpolatedGrad)
    phi -= rAUf * gradpDiffFlux();

    gradpDiffFlux.clear();
}

// RhieChow flux correction for time-steps
if(ddtCorr)
{
    phi += rAUf*fvc::ddtCorr(U, phi);
}


surfaceScalarField buoyancyTerm_RI_Flux
(
    fvc::dotInterpolate(
        mesh.Sf(),
        fvc::reconstruct(Boussinesq.buoyancyTerm() * mesh.magSf() )
    )
);

surfaceScalarField buoyancyTermDiff (Boussinesq.buoyancyTerm()* mesh.magSf() - buoyancyTerm_RI_Flux);
Info << "\t*** buoyancyTermDiff max: " << gMax(buoyancyTermDiff) << tab << "min: " << gMin(buoyancyTermDiff) << endl;