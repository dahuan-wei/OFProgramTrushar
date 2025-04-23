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

#include "temperatureDependentPowerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(temperatureDependentPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        temperatureDependentPowerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::temperatureDependentPowerLaw::calcNu() const
{
	const volScalarField& T = U_.mesh().lookupObject<volScalarField>("T");
	
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            (kzero_-mk_*(T-Tzero_))*pow
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1.0)*strainRate(),
                    dimensionedScalar("VSMALL", dimless, VSMALL)
                ),
                n_.value() - scalar(1.0)
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::temperatureDependentPowerLaw::temperatureDependentPowerLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    temperatureDependentPowerLawCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    n_(temperatureDependentPowerLawCoeffs_.lookup("n")),
    nuMin_(temperatureDependentPowerLawCoeffs_.lookup("nuMin")),
    nuMax_(temperatureDependentPowerLawCoeffs_.lookup("nuMax")),
    kzero_(temperatureDependentPowerLawCoeffs_.lookup("kzero")),
    mk_(temperatureDependentPowerLawCoeffs_.lookup("mk")),
    Tzero_(temperatureDependentPowerLawCoeffs_.lookup("Tzero")),
    
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

bool Foam::viscosityModels::temperatureDependentPowerLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    temperatureDependentPowerLawCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    temperatureDependentPowerLawCoeffs_.lookup("n") >> n_;
    temperatureDependentPowerLawCoeffs_.lookup("nuMin") >> nuMin_;
    temperatureDependentPowerLawCoeffs_.lookup("nuMax") >> nuMax_;
    temperatureDependentPowerLawCoeffs_.lookup("kzero") >> kzero_;
    temperatureDependentPowerLawCoeffs_.lookup("mk") >> mk_;
    temperatureDependentPowerLawCoeffs_.lookup("Tzero") >> Tzero_;

    return true;
}


// ************************************************************************* //
