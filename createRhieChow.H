// Rhie-Chow for pressure gradient correction (uses orthogonal correction only)

// Check for the RhieChow correction
bool RhieChow (spaece.dict().lookupOrDefault<bool>("RhieChow", true));
Info << "RhieChow-pGrad: " << RhieChow << endl;

// cell-centers connecting vector
const surfaceVectorField ed = mesh.delta()()/mag(mesh.delta()());

// surface normal orthogonal gradient object
Foam::fv::orthogonalSnGrad<scalar> faceGradient(mesh);

//----------------------------------------------------------------------------//
// check for RhieChow correction for time step
bool ddtCorr (spaece.dict().lookupOrDefault<bool>("ddtCorr", true));
Info << "RhieChow-ddtCorr: " << ddtCorr << endl;