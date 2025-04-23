/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::arTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    arTransportI.H
    arTransport.C

\*---------------------------------------------------------------------------*/

#ifndef arTransport_H
#define arTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class arTransport;

template<class Thermo>
inline arTransport<Thermo> operator+
(
    const arTransport<Thermo>&,
    const arTransport<Thermo>&
);

template<class Thermo>
inline arTransport<Thermo> operator-
(
    const arTransport<Thermo>&,
    const arTransport<Thermo>&
);

template<class Thermo>
inline arTransport<Thermo> operator*
(
    const scalar,
    const arTransport<Thermo>&
);

template<class Thermo>
inline arTransport<Thermo> operator==
(
    const arTransport<Thermo>&,
    const arTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const arTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class arTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class arTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline arTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline arTransport(const word&, const arTransport&);

        //- Construct from Istream
        arTransport(Istream&);

        //- Construct from dictionary
        arTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<arTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<arTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<arTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "ar<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline arTransport& operator=(const arTransport&);

        inline void operator+=(const arTransport&);

        inline void operator-=(const arTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend arTransport operator+ <Thermo>
        (
            const arTransport&,
            const arTransport&
        );

        friend arTransport operator- <Thermo>
        (
            const arTransport&,
            const arTransport&
        );

        friend arTransport operator* <Thermo>
        (
            const scalar,
            const arTransport&
        );

        friend arTransport operator== <Thermo>
        (
            const arTransport&,
            const arTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const arTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "arTransportI.H"

#ifdef NoRepository
#   include "arTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
