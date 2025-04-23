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

\*---------------------------------------------------------------------------*/

#include "dynamicWaveFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicWaveFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicWaveFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicWaveFvMesh::dynamicWaveFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
 //   amplitude_(readScalar(dynamicMeshCoeffs_.lookup("amplitude"))),
//    frequency_(readScalar(dynamicMeshCoeffs_.lookup("frequency"))),
 //   refPlaneX_(readScalar(dynamicMeshCoeffs_.lookup("refPlaneX"))),  
	// New data for wave motion
    waveAmplitude_(readScalar(dynamicMeshCoeffs_.lookup("waveAmplitude"))),
    waveLength_(readScalar(dynamicMeshCoeffs_.lookup("waveLength"))),
    wavePeriod_(readScalar(dynamicMeshCoeffs_.lookup("wavePeriod"))),
    waveDir_(readScalar(dynamicMeshCoeffs_.lookup("waveDir"))),
    waveType_(readScalar(dynamicMeshCoeffs_.lookup("waveType"))),
    alpha_(readScalar(dynamicMeshCoeffs_.lookup("alpha"))),                  
    stationaryPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
{
    Info<< "Performing a dynamic mesh calculation for wave motion: " << endl ;
    //    << "amplitude: " << amplitude_
    //    << " frequency: " << frequency_
   //     << " refPlaneX: " << refPlaneX_ << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicWaveFvMesh::~dynamicWaveFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicWaveFvMesh::update()
{
    
    pointField newPoints = stationaryPoints_;

	if(waveType_==2)
	{
		 newPoints.replace
		(
			vector::Z,
			stationaryPoints_.component(vector::Z)
		 +
	     	Foam::exp(
						-alpha_*stationaryPoints_.component(vector::Z)
					  )
		 *
		    (
		    waveAmplitude_*cos(
								(constant::mathematical::twoPi)*stationaryPoints_.component(vector::X)/waveLength_
								-
								(constant::mathematical::twoPi)*time().value()/wavePeriod_
							  )
			+
			(
			0.5*(waveAmplitude_*(constant::mathematical::twoPi)/waveLength_)*
			cos(
				2*(
					(constant::mathematical::twoPi)*stationaryPoints_.component(vector::X)/waveLength_
								-
					(constant::mathematical::twoPi)*time().value()/wavePeriod_
				  )
				)
			)
		    )  					  
		);
	
	}
	else if(waveType_==3)
	{
		 newPoints.replace
		(
			vector::Z,
			stationaryPoints_.component(vector::Z)
		 +
	     	Foam::exp(
						-alpha_*stationaryPoints_.component(vector::Z)
					  )
		 *
		    (
		    waveAmplitude_*cos(
								(constant::mathematical::twoPi)*stationaryPoints_.component(vector::X)/waveLength_
								-
								(constant::mathematical::twoPi)*time().value()/wavePeriod_
							  )
			+
			(
			0.5*(waveAmplitude_*(constant::mathematical::twoPi)/waveLength_)*
			cos(
				2*(
					(constant::mathematical::twoPi)*stationaryPoints_.component(vector::X)/waveLength_
								-
					(constant::mathematical::twoPi)*time().value()/wavePeriod_
				  )
				)
			)
			+
			(
			 (3/8)*pow(
					   (0.5*(waveAmplitude_*(constant::mathematical::twoPi)/waveLength_)),2
					  )*
					  cos(
						3*(
							(constant::mathematical::twoPi)*stationaryPoints_.component(vector::X)/waveLength_
										-
							(constant::mathematical::twoPi)*time().value()/wavePeriod_
						  )
						)
			)
		    )  					  
		);
    }
    else
    {
		 newPoints.replace
		(
			vector::Z,
			stationaryPoints_.component(vector::Z)
		 +
	     	Foam::exp(
						-alpha_*stationaryPoints_.component(vector::Z)
					  )
		 *
		    (
		    waveAmplitude_*cos(
								(constant::mathematical::twoPi)*stationaryPoints_.component(vector::X)/waveLength_
								-
								(constant::mathematical::twoPi)*time().value()/wavePeriod_
							  )

		    )  					  
		);		
	}		
   
    fvMesh::movePoints(newPoints);

    volVectorField& U =
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"));
    U.correctBoundaryConditions();

    return true;
}


// ************************************************************************* //
