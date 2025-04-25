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

#include "paraboloidVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

paraboloidVelocityFvPatchVectorField::paraboloidVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(0),
    Rm_(0),
    n_(1, 0, 0),
    y_(0, 1, 0),
    z_(0, 0, 1),
    ctr_(0, 0, 0)
{}


paraboloidVelocityFvPatchVectorField::paraboloidVelocityFvPatchVectorField
(
    const paraboloidVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umax_(ptf.Umax_),
    Rm_(ptf.Rm_),
    n_(ptf.n_),
    y_(ptf.y_),
    z_(ptf.z_),
    ctr_(ptf.ctr_)
{}


paraboloidVelocityFvPatchVectorField::paraboloidVelocityFvPatchVectorField
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
    y_(dict.lookup("y")),
    z_(dict.lookup("z")),
    ctr_(dict.lookup("ctr"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL || mag(z_) < SMALL)
    {
        FatalErrorIn("paraboloidVelocityFvPatchVectorField(dict)")
            << "n or y or z given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);
    z_ /= mag(z_);

    evaluate();
}


paraboloidVelocityFvPatchVectorField::paraboloidVelocityFvPatchVectorField
(
    const paraboloidVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Umax_(fcvpvf.Umax_),
    Rm_(fcvpvf.Rm_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
    z_(fcvpvf.z_),
    ctr_(fcvpvf.ctr_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void paraboloidVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get range and orientation
   //  boundBox bb(patch().patch().localPoints(), true);

    // vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
   // scalarField coord = 2*((c - ctr) & y_)/((bb.max() - bb.min()) & y_);
   
    scalarField coord_y = (c - ctr_) & y_;

    scalarField coord_z = (c - ctr_) & z_;
    
    scalarField coord_r = sqrt (coord_y*coord_y + coord_z*coord_z);
    
    vectorField::operator=(n_*Umax_*(1.0 - sqr(coord_r/Rm_))); // u(r) = Umax (1 - (r/R)^2)
}


// Write
void paraboloidVelocityFvPatchVectorField::write(Ostream& os) const
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
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    os.writeKeyword("ctr")
        << ctr_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, paraboloidVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
