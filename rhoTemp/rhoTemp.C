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

\*---------------------------------------------------------------------------*/

#include "rhoTemp.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::rhoTemp<Specie>::rhoTemp(Istream& is)
:
    Specie(is),
    rho0_(readScalar(is)),
    T0_(readScalar(is))
{
    is.check("rhoTemp<Specie>::rhoTemp(Istream& is)");
}


template<class Specie>
Foam::rhoTemp<Specie>::rhoTemp(const dictionary& dict)
:
    Specie(dict),
    rho0_(readScalar(dict.subDict("equationOfState").lookup("rho0"))),
    T0_(readScalar(dict.subDict("equationOfState").lookup("T0")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::rhoTemp<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("rho0", rho0_);
    dict.add("T0", T0_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const rhoTemp<Specie>& ico)
// os: A reference to an Ostream object, which represents the output stream.
// ico: A constant reference to a rhoTemp object that needs to be serialized.
{
    os  << static_cast<const Specie&>(ico)
        << token::SPACE << ico.rho0_;

    os.check("Ostream& operator<<(Ostream& os, const rhoTemp<Specie>& ico)");
    return os;
}


// ************************************************************************* //
