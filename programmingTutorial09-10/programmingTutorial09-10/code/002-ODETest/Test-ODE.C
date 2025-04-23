/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class myODE
:
    public ODESystem
{

public:

    myODE()
    {}

    label nEqns() const
    {
        return 2;
    }

    void derivatives
    (
        const scalar x,
        const scalarField& y,
        scalarField& dydx
    ) const
    {
        dydx[0] = y[1];
        dydx[1] = y[0];
    }

    void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
        dfdx[0] = 0.0;
        dfdx[1] = 0.0;

        dfdy[0][0] = 0.0;
        dfdy[0][1] = 1.0;

        dfdy[1][0] = 1.0;
        dfdy[1][1] = 0.0;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    
    argList::validArgs.append("ODESolver");
    argList args(argc, argv);

    // Create the ODE system
    myODE ode;

    dictionary dict;
    dict.add("solver", args[1]);

    // Create the selected ODE system solver
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);

    // Initialise the ODE system fields
    scalar xStart = 0.00001;
    scalarField yStart(ode.nEqns());
    
    yStart[0] = 0.0; // y(0) = 0
    yStart[1] = 1.0; //y'(0) = 1
    
    // Print the evolution of the solution and the time-step
    scalarField dyStart(ode.nEqns());
    ode.derivatives(xStart, yStart, dyStart);

    Info<< setw(10) << "relTol" << setw(12) << "dxEst";
    Info<< setw(13) << "dxDid" << setw(14) << "dxNext" << endl;
    Info<< setprecision(6);

    for (label i=0; i<15; i++)
    {
        scalar relTol = ::Foam::exp(-scalar(i + 1));

        scalar x = xStart;
        scalarField y(yStart);

        scalar dxEst = 0.6;
        scalar dxNext = dxEst;

        odeSolver->relTol() = relTol;
        odeSolver->solve(x, y, dxNext);

        Info<< scientific << setw(13) << relTol;
        Info<< fixed << setw(11) << dxEst;
        Info<< setw(13) << x - xStart << setw(13) << dxNext
            << setw(13) << y[0] << setw(13) << y[1]
            << setw(13) << y[2] << setw(13) << y[3]
            << endl;
    }

    scalar x = xStart;
    scalar xEnd = x + 1.0;
    scalarField y(yStart);

    scalarField yEnd(ode.nEqns());
    yEnd[0] = ::Foam::j0(xEnd);
    yEnd[1] = ::Foam::j1(xEnd);

    scalar dxEst = 0.5;

    odeSolver->relTol() = 1e-4;
    odeSolver->solve(x, xEnd, y, dxEst);

    Info      << "Numerical:  " << y << ", dxEst = " << dxEst << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
