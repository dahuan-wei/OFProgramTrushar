/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::dynamicWaveFvMesh

Description
    Mesh motion specifically for the "pumping" system of an ink-jet
    injector.

    The set of points in the "pumping" region are compressed and expanded
    sinusoidally to impose a sinusoidal variation of the flow at the
    nozzle exit.

SourceFiles
    dynamicWaveFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicWaveFvMesh_H
#define dynamicWaveFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "pointIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicWaveFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicWaveFvMesh
:
    public dynamicFvMesh
{
    // Private data

        dictionary dynamicMeshCoeffs_;

  //      scalar amplitude_;
  //      scalar frequency_;
  //      scalar refPlaneX_;

	// New data for wave motion
        scalar waveAmplitude_;
        scalar waveLength_;
        scalar wavePeriod_;
        scalar waveDir_; //either +1,-1 or 0
        scalar waveType_; //1:linear, 2:2nd order, 3: 3rd order
        scalar alpha_; //diffusivity parameter
        
        pointIOField stationaryPoints_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicWaveFvMesh(const dynamicWaveFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicWaveFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicWaveFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicWaveFvMesh(const IOobject& io);


    //- Destructor
    ~dynamicWaveFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
