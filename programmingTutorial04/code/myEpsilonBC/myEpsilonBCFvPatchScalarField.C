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

#include "myEpsilonBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myEpsilonBCFvPatchScalarField::myEpsilonBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    a_(0),
    alpha_(0),
    Iu_(0),
    H_(0),
    n_(1, 0, 0),
    y_(0, 1, 0)
{}


myEpsilonBCFvPatchScalarField::myEpsilonBCFvPatchScalarField
(
    const myEpsilonBCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    a_(ptf.a_),
    alpha_(ptf.alpha_),
    Iu_(ptf.Iu_),
    H_(ptf.H_),
    n_(ptf.n_),
    y_(ptf.y_)
{}


myEpsilonBCFvPatchScalarField::myEpsilonBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    a_(readScalar(dict.lookup("a"))),
    alpha_(readScalar(dict.lookup("alpha"))),
    Iu_(readScalar(dict.lookup("Iu"))),
    H_(readScalar(dict.lookup("H"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y"))
{
    if (mag(n_) < SMALL || mag(y_) < SMALL)
    {
        FatalErrorIn("myEpsilonBCFvPatchScalarField(dict)")
            << "n or y given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);

    evaluate();
}


myEpsilonBCFvPatchScalarField::myEpsilonBCFvPatchScalarField
(
    const myEpsilonBCFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    a_(fcvpvf.a_),
    alpha_(fcvpvf.alpha_),
    Iu_(fcvpvf.Iu_),
    H_(fcvpvf.H_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myEpsilonBCFvPatchScalarField::updateCoeffs()
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
   
    scalarField coord = c & y_;
    
    scalarField Uinf = a_*pow(coord,alpha_);  
    
    scalarField TKE = 1.5*sqr(Uinf*Iu_);  
    
    scalarField epsilon = pow(0.09,0.75)*pow(TKE,1.5)/H_;  

    scalarField::operator=(epsilon); 
}


// Write
void myEpsilonBCFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("a")
        << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha")
        << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("Iu")
        << Iu_ << token::END_STATEMENT << nl;
    os.writeKeyword("H")
        << H_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, myEpsilonBCFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
