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

\*---------------------------------------------------------------------------*/

#include "myInletPulsatileBCFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myInletPulsatileBCFvPatchVectorField::myInletPulsatileBCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umean_(0.0),
    period_(0.0)
{}


myInletPulsatileBCFvPatchVectorField::myInletPulsatileBCFvPatchVectorField
(
    const myInletPulsatileBCFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umean_(ptf.Umean_),
    period_(ptf.period_)
{}


myInletPulsatileBCFvPatchVectorField::myInletPulsatileBCFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umean_(readScalar(dict.lookup("Umean"))),
    period_(readScalar(dict.lookup("period")))
{
}


myInletPulsatileBCFvPatchVectorField::myInletPulsatileBCFvPatchVectorField
(
    const myInletPulsatileBCFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Umean_(fcvpvf.Umean_),
    period_(fcvpvf.period_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myInletPulsatileBCFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar curTime_ = this->db().time().value();

    scalar Umean_ = 0.511;
    scalar period_ = 0.928;
    scalar pi = 3.141592654;

    scalar UTime1 = ((-((Umean_) + (-0.0080467)*cos(1 *2.*pi*curTime_/period_) +
                  (0.086438)*sin(1 *2.*pi*curTime_/period_) +
                  (-0.0057042)*cos(2 *2.*pi*curTime_/period_) +
                  (0.036684)*sin(2 *2.*pi*curTime_/period_) + 
                  (-0.037873)*cos(3 *2.*pi*curTime_/period_) + 
                  (0.023462)*sin(3 *2.*pi*curTime_/period_) +
                  (-0.019703)*cos(4 *2.*pi*curTime_/period_) +
                  (0.0023512)*sin(4 *2.*pi*curTime_/period_) +
                  (-0.019301)*cos(5 *2.*pi*curTime_/period_) +
                  (0.006191)*sin(5 *2.*pi*curTime_/period_) +
                  (-0.015524)*cos(6 *2.*pi*curTime_/period_) +
                  (-0.0082778)*sin(6 *2.*pi*curTime_/period_) +
                  (-0.013376)*cos(7 *2.*pi*curTime_/period_) +
                  (-0.007184)*sin(7 *2.*pi*curTime_/period_) +
                  (8.3728e-005)*cos(8 *2.*pi*curTime_/period_) +
                  (-0.0034172)*sin(8 *2.*pi*curTime_/period_) +
                  (-0.0061556)*cos(9 *2.*pi*curTime_/period_) +
                  (-0.0056501)*sin(9 *2.*pi*curTime_/period_) +
                  (0.0014531)*cos(10*2.*pi*curTime_/period_) +
                  (-0.0035759)*sin(10*2.*pi*curTime_/period_))
           )/ Umean_ * 0.205);

    tmp<vectorField> n = patch().nf();

    vectorField::operator=(n*UTime1);
}


// Write
void myInletPulsatileBCFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Umean")
        << Umean_ << token::END_STATEMENT << nl;
    os.writeKeyword("period")
        << period_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, myInletPulsatileBCFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
