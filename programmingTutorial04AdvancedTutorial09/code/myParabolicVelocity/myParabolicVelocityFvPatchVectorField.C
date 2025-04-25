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

#include "myParabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
1)	parabolic Velocity

		u(r) = 2 Vavg (1 - (r/R)^2)
		Umax = 2 Vavg

		u(r) = Umax (1 - (r/R)^2)
*/
namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myParabolicVelocityFvPatchVectorField::myParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(0),
    Rm_(0),
    n_(1, 0, 0),
    y_(0, 1, 0)
{}


myParabolicVelocityFvPatchVectorField::myParabolicVelocityFvPatchVectorField
(
    const myParabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umax_(ptf.Umax_),
    Rm_(ptf.Rm_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


myParabolicVelocityFvPatchVectorField::myParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(readScalar(dict.lookup("Umax"))),
    Rm_(readScalar(dict.lookup("Rm"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("myParabolicVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    evaluate();
}


myParabolicVelocityFvPatchVectorField::myParabolicVelocityFvPatchVectorField
(
    const myParabolicVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Umax_(fcvpvf.Umax_),
    Rm_(fcvpvf.Rm_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myParabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get range and orientation
  //  boundBox bb(patch().patch().localPoints(), true);

  //  vector ctr = 0.5*(bb.max() + bb.min());//center 

    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
   // scalarField coord = 2*((c - ctr) & y_)/((bb.max() - bb.min()) & y_);
   
    scalarField coord_r = c & y_;//the y coord or the each  cell in center

    vectorField::operator=(n_*Umax_*(1.0 - sqr(coord_r/Rm_))); // u(r) = Umax (1 - (r/R)^2)
}


// Write
void myParabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Umax")
        << Umax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Rm")
        << Rm_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, myParabolicVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
