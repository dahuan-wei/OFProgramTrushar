    Info<< "Reading incremental displacement field dU\n" << endl;
    volVectorField dU
    (
        IOobject
        (
            "dU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    Info<< "Reading displacement field dU\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dU
    );
    
    Info<< "Reading stress field sigma\n" << endl;
    
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading incremental plasic strain dpstr\n" << endl;
    volTensorField dpstr
    (
        IOobject
        (
            "dpstr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

/*
    volScalarField* Tptr = NULL;

    if (thermalStress)
    {
        Info<< "Reading field T\n" << endl;
        Tptr = new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }

    volScalarField& T = *Tptr;
*/
