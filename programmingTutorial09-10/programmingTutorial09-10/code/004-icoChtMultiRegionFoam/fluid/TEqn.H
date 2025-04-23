{
    //alphat = turbulence->nut()/Prt;
   // alphat.correctBoundaryConditions();

    volScalarField alphaEff
    (
    "alphaEff", 
    turb.nu()/Pr + turb.nut()/Prt
    );

  //  while(nNonOrthCorr)
    {
		fvScalarMatrix TEqn
		(
			fvm::ddt(T)
		  + fvm::div(phi, T)
		  - fvm::laplacian(alphaEff, T)
		);

		TEqn.relax();

		//fvOptions.constrain(TEqn);

		TEqn.solve();
	}

	   // radiation->correct();

	   // fvOptions.correct(T);

	  // rhok = 1.0 - beta*(T - TRef);
}
