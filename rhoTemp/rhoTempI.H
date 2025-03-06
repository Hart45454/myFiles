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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoTemp<Specie>::rhoTemp
(
    const Specie& sp,
    const scalar rho0,
    const scalar T0
)
:
    Specie(sp),
    rho0_(rho0),
    T0_(T0)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoTemp<Specie>::rhoTemp
(
    const word& name,
    const rhoTemp<Specie>& ico
)
:
    Specie(name, ico),
    rho0_(ico.rho0_),
    T0_(ico.T0_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoTemp<Specie> >
Foam::rhoTemp<Specie>::clone() const
{
    return autoPtr<rhoTemp<Specie> >(new rhoTemp<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::rhoTemp<Specie> >
Foam::rhoTemp<Specie>::New(Istream& is)
{
    return autoPtr<rhoTemp<Specie> >(new rhoTemp<Specie>(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoTemp<Specie>::rho(scalar p, scalar T) const
{
    return rho0_ + ((T - T0_) - 0.833 * pow(10,-3));
}


template<class Specie>
inline Foam::scalar Foam::rhoTemp<Specie>::psi(scalar, scalar T) const
{
    return 0.0;
}

template<class Specie>
inline Foam::scalar Foam::rhoTemp<Specie>::h(scalar p, scalar T) const
{
    return 0.0;
}

template<class Specie>
inline Foam::scalar Foam::rhoTemp<Specie>::cp(scalar p, scalar T) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTemp<Specie>::Z(scalar, scalar) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::rhoTemp<Specie>::cpMcv(scalar, scalar) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoTemp<Specie>::operator+=(const rhoTemp<Specie>& ico)
{
    scalar molr1 = this->nMoles();

    Specie::operator+=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho0_ = molr1*rho0_ + molr2*ico.rho0_;
}


template<class Specie>
inline void Foam::rhoTemp<Specie>::operator-=(const rhoTemp<Specie>& ico)
{
    scalar molr1 = this->nMoles();

    Specie::operator-=(ico);

    molr1 /= this->nMoles();
    scalar molr2 = ico.nMoles()/this->nMoles();

    rho0_ = molr1*rho0_ - molr2*ico.rho0_;
}


template<class Specie>
inline void Foam::rhoTemp<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoTemp<Specie> Foam::operator+
(
    const rhoTemp<Specie>& ico1,
    const rhoTemp<Specie>& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return rhoTemp<Specie>
    (
        static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2),
        molr1*ico1.rho0_ + molr2*ico2.rho0_
    );
}


template<class Specie>
inline Foam::rhoTemp<Specie> Foam::operator-
(
    const rhoTemp<Specie>& ico1,
    const rhoTemp<Specie>& ico2
)
{
    scalar nMoles = ico1.nMoles() + ico2.nMoles();
    scalar molr1 = ico1.nMoles()/nMoles;
    scalar molr2 = ico2.nMoles()/nMoles;

    return rhoTemp<Specie>
    (
        static_cast<const Specie&>(ico1)
      - static_cast<const Specie&>(ico2),
        molr1*ico1.rho0_ - molr2*ico2.rho0_
    );
}


template<class Specie>
inline Foam::rhoTemp<Specie> Foam::operator*
(
    const scalar s,
    const rhoTemp<Specie>& ico
)
{
    return rhoTemp<Specie>(s*static_cast<const Specie&>(ico), ico.rho0_);
}


template<class Specie>
inline Foam::rhoTemp<Specie> Foam::operator==
(
    const rhoTemp<Specie>& ico1,
    const rhoTemp<Specie>& ico2
)
{
    return ico2 - ico1;
}


// ************************************************************************* //
