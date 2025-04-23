/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    stressedFoam

Description
    steady-state segregated finite-volume solver of linear-elastic 
    + perfect plastic, small-strain deformation of a solid body

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "readMechanicalProperties.H"
//#   include "readThermalProperties.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Iteration: " << runTime.timeName() << nl << endl;

#       include "readStressedFoamControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;

        do
        {
            volTensorField graddU = fvc::grad(dU);
			/*
            if (thermalStress)
            {
                solve
                (
                    fvm::ddt(T) == fvm::laplacian(DT, T)
                );
            }
            */ 

            fvVectorMatrix dUEqn
            (
                fvm::laplacian(2*mu + lambda, dU, "laplacian(DdU,dU)")                
             ==
              - fvc::div
                (
                    mu*graddU.T() + lambda*(I*tr(graddU)) - (mu + lambda)*graddU,
                    "div(sigmaExp)"
                )
                
              + fvc::div
                (
                    2*mu*dpstr + lambda*I*tr(dpstr),
                    "div(sigmaP)"
                )
            );

            /*
            if (thermalStress)
            {
                UEqn += threeKalpha*fvc::grad(T);
            }
            */ 

            //UEqn.setComponentReference(1, 0, vector::X, 0);
            //UEqn.setComponentReference(1, 0, vector::Z, 0);

            initialResidual = dUEqn.solve().initialResidual();

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
        
        U += dU;

#       include "calculateStress.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    //deleteDemandDrivenData(Tptr);

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
