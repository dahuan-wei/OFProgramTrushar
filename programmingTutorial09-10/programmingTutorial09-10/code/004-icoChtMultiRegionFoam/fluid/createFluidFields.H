// Initialise fluid field pointer lists

PtrList<singlePhaseTransportModel> thermoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<dimensionedScalar> nuFluid(fluidRegions.size());
PtrList<dimensionedScalar> PrFluid(fluidRegions.size());
PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());
PtrList<incompressible::turbulenceModel> turbulence(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());

PtrList<fv::IOoptionList> fluidFvOptions(fluidRegions.size());

/*
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());

PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::IOoptionList> fluidFvOptions(fluidRegions.size());
*/


// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;
/*
    Info<< "    Adding to thermoFluid\n" << endl;

    thermoFluid.set
    (
        i,
        rhoThermo::New(fluidRegions[i]).ptr()
    );
*/

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );
    
    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(UFluid[i])
                & fluidRegions[i].Sf()
        )
    );
    
    Info<< "    Adding to pFluid\n" << endl;
    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to TFluid\n" << endl;
    TFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );


    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "K",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );    

    Info<< "    Adding to thermoFluid\n" << endl;

    thermoFluid.set
    (
        i,
        new singlePhaseTransportModel(UFluid[i], phiFluid[i])
    );
    
        

    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        incompressible::turbulenceModel::New
        (
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );
    
    /*
        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<incompressible::turbulenceModel>
            (
            	incompressible::turbulenceModel::New
                (
                    UFluid[i],
                    phiFluid[i],
                    thermoFluid[i]
                )
            ).ptr()
        );
        */
    
        setRefCell
        (
            pFluid[i],
            fluidRegions[i].solutionDict().subDict("PIMPLE"),
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
        
   // const dictionary& pimpleDict =
   //     fluidRegions[i].solutionDict().subDict("PIMPLE");
    
        nuFluid.set
    (
		i,
        new dimensionedScalar(thermoFluid[i].lookup("nu"))	
    );
    
        PrFluid.set
    (
		i,
        new dimensionedScalar(thermoFluid[i].lookup("Pr"))	
    );
        PrtFluid.set
    (
		i,
        new dimensionedScalar(thermoFluid[i].lookup("Prt"))	
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::IOoptionList(fluidRegions[i])
    );
    ////////////////////////////////////////////////////////////////
    /*Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );


    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar("hRef", dimLength, 0)
        )
    );

    dimensionedScalar ghRef
    (
        mag(gFluid[i].value()) > SMALL
      ? gFluid[i]
          & (cmptMag(gFluid[i].value())/mag(gFluid[i].value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", gFluid[i].dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    radiation.set
    (
        i,
        radiation::radiationModel::New(thermoFluid[i].T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                "dpdt",
                thermoFluid[i].p().dimensions()/dimTime,
                0
            )
        )
    );

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");
    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::IOoptionList(fluidRegions[i])
    );
    */ 
}
