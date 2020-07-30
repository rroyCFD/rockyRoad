// Compute the full pressure and adjust the pressure level.

// Note background pressure is updated after TEqn solve (may be lagged in a non-iterative algorithm)
// Boussinesq.updateBackgroundPressure();

p = p_rgh + Boussinesq.backgroundPressure();

if (pRefOn)
{
	// change in Buoyant background pressure may alter pRefValue at pRefCell
	// since the linear system enforces pRefCell value within the solver tolerance
	// pp[pRefCell] not-equal-to 0.0
    Info << "p_rgh[pRefCell] (before adjustment): " << getRefCellValue(p_rgh, pRefCell) << endl;
    
    dimensionedScalar pAdjust(
    	"pAdjust",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell));

    p     += pAdjust;
    p_rgh += pAdjust;

    Info << "p_rgh[pRefCell] (after  adjustment): " << getRefCellValue(p_rgh, pRefCell) << endl;
}
