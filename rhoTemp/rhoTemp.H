/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::rhoTemp

Description
    RhoConst (rho = const) of state.

SourceFiles
    rhoTempI.H
    rhoTemp.C

\*---------------------------------------------------------------------------*/

#ifndef rhoTemp_H
#define rhoTemp_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class rhoTemp;

template<class Specie>
inline rhoTemp<Specie> operator+
(
    const rhoTemp<Specie>&,
    const rhoTemp<Specie>&
);

template<class Specie>
inline rhoTemp<Specie> operator-
(
    const rhoTemp<Specie>&,
    const rhoTemp<Specie>&
);

template<class Specie>
inline rhoTemp<Specie> operator*
(
    const scalar,
    const rhoTemp<Specie>&
);

template<class Specie>
inline rhoTemp<Specie> operator==
(
    const rhoTemp<Specie>&,
    const rhoTemp<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const rhoTemp<Specie>&
);


/*---------------------------------------------------------------------------*\
                        Class rhoTemp Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class rhoTemp
:
    public Specie
{
    // Private data

        //- Density
        scalar rho0_;

        //- Reference Temperature
        scalar T0_;


public:

    // Constructors

        //- Construct from components
        inline rhoTemp
        
        (
            const Specie& sp, 
            const scalar rho0,
            const scalar T0
            
        );

	//- Construct from Boussinesq
        inline rhoTemp(const rhoTemp& sp);

        //- Construct from Istream
        rhoTemp(Istream&);

        //- Construct from dictionary
        rhoTemp(const dictionary& dict);

        //- Construct as named copy
        inline rhoTemp(const word& name, const rhoTemp&);

        //- Construct and return a clone
        inline autoPtr<rhoTemp> clone() const;

        // Selector from Istream
        inline static autoPtr<rhoTemp> New(Istream& is);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "rhoTemp<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Added section -
	    inline scalar h(const scalar p, const scalar T) const;

	    inline scalar cp(scalar p, scalar T) const;
	    //-----------------

	    //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kmol K]
            inline scalar cpMcv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const rhoTemp&);
        inline void operator-=(const rhoTemp&);

        inline void operator*=(const scalar);


    // Friend operators

        friend rhoTemp operator+ <Specie>
        (
            const rhoTemp&,
            const rhoTemp&
        );

        friend rhoTemp operator- <Specie>
        (
            const rhoTemp&,
            const rhoTemp&
        );

        friend rhoTemp operator* <Specie>
        (
            const scalar s,
            const rhoTemp&
        );

        friend rhoTemp operator== <Specie>
        (
            const rhoTemp&,
            const rhoTemp&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const rhoTemp&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rhoTempI.H"

#ifdef NoRepository
#   include "rhoTemp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
