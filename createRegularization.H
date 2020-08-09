// RegularizationModel
RegularizationModel* RegModel;

// Check the convection-term-regularization on or off
bool regOn (spaece.dict().lookupOrDefault<bool>("regularizationOn", false));
Info << "regularizationOn: " << regOn << endl;

if(regOn)
{
    RegModel = new RegularizationModel(p_rgh, pRefCell, pRefValue, spaece.nNonOrthCorr());
    RegModel->setRegOn(regOn);
}
else
{
    delete RegModel;
}
//------------------------------------------------------------------------//