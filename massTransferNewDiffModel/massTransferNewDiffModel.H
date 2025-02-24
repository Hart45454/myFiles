/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria

    Description
    Correlation for Sherwood number according to
    Gunn, D. J. International Journal of Heat and Mass Transfer 21.4 (1978)

\*---------------------------------------------------------------------------*/

#ifndef massTransferNewDiffModel_H
#define massTransferNewDiffModel_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "massTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class massTransferNewDiffModel Declaration
\*---------------------------------------------------------------------------*/

class massTransferNewDiffModel
:
    public massTransferModel
{
protected:

    // Configuration and Flags:

    dictionary propsDict_;

    bool multiTypes_;

    //bool expSherwood_;

    bool interpolation_;

    bool verbose_;

    bool implicit_;

    // Field Names and Properties:

    word SPartFluidName_;

    volScalarField SPartFluid_;

    word SPartFluidCoeffName_;

    volScalarField SPartFluidCoeff_;

    // Transport Properties:

    volScalarField ReField_;

    //volScalarField ShField_;

    // Reference Fields:

    const word tempFieldName_;

    const volScalarField& tempField_;           // ref to temperature field

    word concFieldName_;

    const volScalarField& concField_;           // ref to concentration field

    word satConcFieldName_;

    const volScalarField& satConcField_;        // ref to saturation concentration field

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;        // ref to voidfraction field

    // Other Parameters:

    scalar maxSource_;                          // max (limited) value of src field

    // Fluid Properties:
    
    word velFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;

    // Other Parameters:

    bool coupleDEM_;

    const word partMassFluxName_;

    const word partMassFluxCoeffRegName_;

    const word partReRegName_;

    const word partShRegName_;

    mutable scalar scaleDia_;

    scalarList typeCG_;

    // Utility Functions:

    //scalar Sherwood(scalar, scalar, scalar) const;

    virtual void giveData();

    virtual void massFlux(label, scalar, scalar, scalar, scalar, scalar cg3 = 1.0);

public:

    //- Runtime type information
    TypeName("massTransferNewDiffModel");

    // Constructors

        //- Construct from components
        massTransferNewDiffModel
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~massTransferNewDiffModel();


    // Member Functions

        void addMassContribution(volScalarField&) const;

        void addMassCoefficient(volScalarField&) const;

        void calcMassContribution();

        void postFlow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

