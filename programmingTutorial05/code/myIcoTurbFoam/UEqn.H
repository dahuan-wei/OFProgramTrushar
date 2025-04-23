// Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
 //         - fvm::laplacian(nu, U)
          + turbulence->divDevReff(U)
        );

        UEqn.relax();
        
        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }
