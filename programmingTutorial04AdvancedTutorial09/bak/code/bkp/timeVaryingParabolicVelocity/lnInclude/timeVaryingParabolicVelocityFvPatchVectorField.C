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

#include "timeVaryingParabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingParabolicVelocityFvPatchVectorField::timeVaryingParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(0),
    Rm_(0),
    omega_(0),
    n_(1, 0, 0),
    y_(0, 1, 0)
{}


timeVaryingParabolicVelocityFvPatchVectorField::timeVaryingParabolicVelocityFvPatchVectorField
(
    const timeVaryingParabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Umax_(ptf.Umax_),
    Rm_(ptf.Rm_),
    omega_(ptf.omega_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


timeVaryingParabolicVelocityFvPatchVectorField::timeVaryingParabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Umax_(readScalar(dict.lookup("Umax"))),
    Rm_(readScalar(dict.lookup("Rm"))),
    omega_(readScalar(dict.lookup("omega"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("timeVaryingParabolicVelocityFvPatchVectorField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    evaluate();
}


timeVaryingParabolicVelocityFvPatchVectorField::timeVaryingParabolicVelocityFvPatchVectorField
(
    const timeVaryingParabolicVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Umax_(fcvpvf.Umax_),
    Rm_(fcvpvf.Rm_),
    omega_(fcvpvf.omega_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingParabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get range and orientation
  //  boundBox bb(patch().patch().localPoints(), true);

  //  vector ctr = 0.5*(bb.max() + bb.min());

    const vectorField& c = patch().Cf();

    // Calculate local 1-D coordinate for the parabolic profile
   // scalarField coord = 2*((c - ctr) & y_)/((bb.max() - bb.min()) & y_);
   
    scalar time_curr = this->db().time().value();
    
    scalarField coord_r = c & y_;

    vectorField::operator=(n_*Umax_*(1.0 - sqr(coord_r/Rm_))*sin(omega_*time_curr)); // u(r) = Umax (1 - (r/R)^2)*sin(omega*time)
}


// Write
void timeVaryingParabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Umax")
        << Umax_ << token::END_STATEMENT << nl;
    os.writeKeyword("Rm")
        << Rm_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega")
        << omega_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, timeVaryingParabolicVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
