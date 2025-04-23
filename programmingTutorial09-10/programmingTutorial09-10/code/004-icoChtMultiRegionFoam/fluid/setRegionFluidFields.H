    fvMesh& mesh = fluidRegions[i];

  //  singlePhaseTransportModel& thermo = thermoFluid[i];
    
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& p = pFluid[i];
    volScalarField& T = TFluid[i];
//    volScalarField& K = KFluid[i];
    
    incompressible::turbulenceModel& turb = turbulence[i];
    
//    dimensionedScalar& nu = nuFluid[i];
    dimensionedScalar& Pr = PrFluid[i];
    dimensionedScalar& Prt = PrtFluid[i];
    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];
    
    fv::IOoptionList& fvOptions = fluidFvOptions[i];
    
    /*
PtrList<label> pRefCellFluid(fluidRegions.size());
PtrList<scalar> pRefValueFluid(fluidRegions.size());

PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<dimensionedScalar> nuFluid(fluidRegions.size());
PtrList<dimensionedScalar> PrFluid(fluidRegions.size());
PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());

/////////////////////////////////////////////////
    thermo.validate(args.executable(), "h", "e");

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];
    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& g = gFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    radiation::radiationModel& rad = radiation[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::IOoptionList& fvOptions = fluidFvOptions[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    const bool frozenFlow = frozenFlowFluid[i];
*/
