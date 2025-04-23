
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
		fvScalarMatrix TEqn
		(
			fvm::ddt(rho*cp,T)
		  - fvm::laplacian(K, T)
		);

		TEqn.relax();
		TEqn.solve();
	}
