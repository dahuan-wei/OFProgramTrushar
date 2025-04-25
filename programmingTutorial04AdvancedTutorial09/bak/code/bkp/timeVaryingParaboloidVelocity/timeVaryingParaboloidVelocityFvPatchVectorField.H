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

Class
    timeVaryingParaboloidVelocityFvPatchVectorField

Description
    Boundary condition specifies a parabolic velocity inlet profile
    (fixed value), given maximum velocity value (peak of the parabola),
    flow direction n and direction of the parabolic coordinate y

SourceFiles
    timeVaryingParaboloidVelocityFvPatchVectorField.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingParaboloidVelocityFvPatchVectorField_H
#define timeVaryingParaboloidVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class timeVaryingParaboloidVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingParaboloidVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Peak velocity magnitude
        
        scalar Umax_;
                
        scalar Rm_;
        
        scalar omega_;
        
        //- Flow direction
        vector n_;

        //- Direction of the y-coordinate
        vector y_;
        vector z_;
        vector ctr_;


public:

    //- Runtime type information
    TypeName("timeVaryingParaboloidVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingParaboloidVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingParaboloidVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingParaboloidVelocityFvPatchVectorField
        //  onto a new patch
        timeVaryingParaboloidVelocityFvPatchVectorField
        (
            const timeVaryingParaboloidVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingParaboloidVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingParaboloidVelocityFvPatchVectorField
        (
            const timeVaryingParaboloidVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingParaboloidVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions
               
        scalar& Umax()
        {
            return Umax_;
        }
        
        scalar& Rm()
        {
            return Rm_;
        }
        scalar& omega()
        {
            return omega_;
        }

        //- Return flow direction
        vector& n()
        {
            return n_;
        }

        //- Return y direction
        vector& y()
        {
            return y_;
        }
        
        vector& z()
        {
            return z_;
        }        
     
        vector& ctr()
        {
            return ctr_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
