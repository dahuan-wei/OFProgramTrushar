/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "myPolynomVelocityPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myPolynomVelocityPointPatchVectorField::
myPolynomVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    origin_(vector::zero),
    p0_(p.localPoints()),
    X2_(0.0),
    X1_(0.0),
    Y2_(0.0),
    Y1_(0.0),
    Cconst_(0.0),
    xAxis_(vector::zero),
    yAxis_(vector::zero),
    periodic_(0.0),
    defTime_(0.0)
{}


myPolynomVelocityPointPatchVectorField::
myPolynomVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    origin_(dict.lookup("origin")),
    X2_(readScalar(dict.lookup("X2"))),
    X1_(readScalar(dict.lookup("X1"))),
    Y2_(readScalar(dict.lookup("Y2"))),
    Y1_(readScalar(dict.lookup("Y1"))),
    Cconst_(readScalar(dict.lookup("Cconst"))),
    xAxis_(dict.lookup("xAxis")),
    yAxis_(dict.lookup("yAxis")),
    periodic_(readScalar(dict.lookup("periodic"))),
    defTime_(readScalar(dict.lookup("defTime")))
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = p.localPoints();
    }
}


myPolynomVelocityPointPatchVectorField::
myPolynomVelocityPointPatchVectorField
(
    const myPolynomVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    p0_(ptf.p0_),
    X2_(ptf.X2_),
    X1_(ptf.X1_),
    Y2_(ptf.Y2_),
    Y1_(ptf.Y1_),
    Cconst_(ptf.Cconst_),
    xAxis_(ptf.xAxis_),
    yAxis_(ptf.yAxis_),
    periodic_(ptf.periodic_),
    defTime_(ptf.defTime_)
{}


myPolynomVelocityPointPatchVectorField::
myPolynomVelocityPointPatchVectorField
(
    const myPolynomVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    origin_(ptf.origin_),
    p0_(ptf.p0_),
    X2_(ptf.X2_),
    X1_(ptf.X1_),
    Y2_(ptf.Y2_),
    Y1_(ptf.Y1_),
    Cconst_(ptf.Cconst_),
    xAxis_(ptf.xAxis_),
    yAxis_(ptf.yAxis_),
    periodic_(ptf.periodic_),
    defTime_(ptf.defTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myPolynomVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->dimensionedInternalField().mesh()();
    const Time& t = mesh.time();
    const pointPatch& p = this->patch();

    vectorField p0Rel = p0_ - origin_;
    vector zAxis = xAxis_ ^ yAxis_;
    vector xAxisOrg = vector(1, 0, 0);//Original axis used for reference
    vector yAxisOrg = vector(0, 1, 0);
    vector zAxisOrg = vector(0, 0, 1);
   
    // Euler angles start
    vector Nline = (xAxisOrg ^ yAxisOrg) ^ (xAxis_ ^ yAxis_);
    scalar alpha = acos(xAxisOrg & Nline);///(mag(xAxisOrg)*mag(Nline)));
    scalar beta = acos(zAxisOrg & zAxis);///(mag(zAxisOrg)*mag(zAxis)));
    scalar gamma = acos(Nline & xAxis_);///(mag(Nline)*mag(xAxis_)));
    scalar Rrot1(cos(alpha)*cos(gamma)-sin(alpha)*cos(beta)*sin(gamma));
    scalar Rrot2(-cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma));
    scalar Rrot3(sin(beta)*sin(alpha));
    scalar Rrot4(sin(alpha)*cos(gamma)+cos(alpha)*cos(beta)*sin(gamma));
    scalar Rrot5(-sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma));
    scalar Rrot6(-sin(beta)*cos(alpha));
    scalar Rrot7(sin(beta)*sin(gamma));
    scalar Rrot8(sin(beta)*cos(gamma));
    scalar Rrot9(cos(beta));
    // Rotation matrix created
    tensor Rrot(Rrot1, Rrot2, Rrot3, Rrot4, Rrot5, Rrot6, Rrot7, Rrot8, Rrot9);
    tensor RrotInv = inv(Rrot);
    vector p0rot;
    vectorField sd=p0Rel;
    forAll(p0_,iter)
    {
        p0rot = p0Rel[iter] & Rrot; // p relative to new origin rotated
        // Plane from x and y values calculated and inserted into z values
        p0rot = vector(0, 0, X2_*p0rot[0]*p0rot[0]+X1_*p0rot[0]+Y2_*p0rot[1]*p0rot[1]+Y1_*p0rot[1]+Cconst_);
        sd[iter] = p0rot & RrotInv; // Plane rotated back to original position
    };
    scalar multipl = 1;
    if ( periodic_ == 1 ) // For periodic b.c.
    {
	if ((int)floor(t.value()/defTime_)% 2  != 0) multipl = -1; // Revese motion for periodic b.c.
    }
    else if ((periodic_ == 0) && (t.value()> defTime_)) multipl = 0; // No motion
    vectorField::operator=
    (
        sd *multipl / defTime_ 
    );

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void myPolynomVelocityPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeKeyword("origin")
        << origin_ << token::END_STATEMENT << nl;
    p0_.writeEntry("p0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    myPolynomVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
