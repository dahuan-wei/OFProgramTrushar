    const dictionary& stressControl =
        mesh.solutionDict().subDict("nonLinearStressedFoam");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("dU")));
