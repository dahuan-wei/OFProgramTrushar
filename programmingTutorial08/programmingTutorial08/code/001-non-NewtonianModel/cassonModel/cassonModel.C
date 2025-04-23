/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "cassonModel.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(cassonModel, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        cassonModel,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::cassonModel::calcNu() const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            pow
              (
				pow(
						tau0_
						/
						max
						(
							strainRate(),
							dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
						)
					,
					0.5
				   )
				   +
				   pow(m_,0.5)
              ,
              scalar(2)
              )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::cassonModel::cassonModel
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    cassonModelCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nuMin_(cassonModelCoeffs_.lookup("nuMin")),
    nuMax_(cassonModelCoeffs_.lookup("nuMax")),
    m_(cassonModelCoeffs_.lookup("m")),
    tau0_(cassonModelCoeffs_.lookup("tau0")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::cassonModel::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    cassonModelCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    cassonModelCoeffs_.lookup("nuMin") >> nuMin_;
    cassonModelCoeffs_.lookup("nuMax") >> nuMax_;
    cassonModelCoeffs_.lookup("m") >> m_;
    cassonModelCoeffs_.lookup("tau0") >> tau0_;

    return true;
}


// ************************************************************************* //
