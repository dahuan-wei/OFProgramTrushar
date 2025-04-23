    if (runTime.outputTime())
    {
        volTensorField graddU = fvc::grad(dU);

        volSymmTensorField sigma_old = sigma;
        
        sigma +=  2.0*mu*symm(graddU) + lambda*I*tr(graddU);

        /*
        if (thermalStress)
        {
            sigma = sigma - I*(rho*threeKalpha*T);
        }
        */ 

     //Check the yield condition
        volScalarField sqrtJ2 = sqrt((1.0/2.0)*magSqr(dev(sigma)));

        volScalarField fac = sqrtJ2/Ys;

        forAll(fac, celli) 
        {
            if (fac[celli] > 1.0)         //Plasticity occurs
             {
                sigma[celli] = 1.0/3.0*I*tr(sigma[celli]) + dev(sigma[celli])/fac[celli];
                symmTensor dsigma = sigma[celli] - sigma_old[celli];
                tensor deps_p = 1.0/3.0*I*tr(dsigma)/(3.0*lambda+2.0*mu).value()
                                + dev(dsigma)/(2.0*mu.value());
                tensor deps = 1.0/2.0*(graddU[celli] + graddU[celli].T());
                dpstr[celli] = deps - deps_p;
             }
             else                         // only elasticity
             {
                dpstr[celli] = tensor::zero;
             }
         }  
        volScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
        );

        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;

        volScalarField sigmaxx
        (
            IOobject
            (
                "sigmaxx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XX)
        );

        volScalarField sigmayy
        (
            IOobject
            (
                "sigmayy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::YY)
        );

        volScalarField sigmazz
        (
            IOobject
            (
                "sigmazz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::ZZ)
        );

        Info<< "Max sigmazz = " << max(sigmazz).value()
            << endl;

        volScalarField sigmaxy
        (
            IOobject
            (
                "sigmaxy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XY)
        );

        volScalarField sigmaxz
        (
            IOobject
            (
                "sigmaxz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XZ)
        );

        volScalarField sigmayz
        (
            IOobject
            (
                "sigmayz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::YZ)
        );

        runTime.write();
    }
