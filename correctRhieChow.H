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
