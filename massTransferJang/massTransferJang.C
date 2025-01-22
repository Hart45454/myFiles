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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "massTransferJang.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(massTransferJang, 0);

addToRunTimeSelectionTable(massTransferModel, massTransferJang, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
massTransferJang::massTransferJang
(
    const dictionary& dict,
    // Contains the configuration and properties for the mass transfer model, typically 
    // specified in a simulation setup file.
    cfdemCloudEnergy& sm
    // A reference to the simulation manager object, which contains information about the 
    // CFD-DEM system (e.g., mesh, time steps, fields, etc.).
)
:
// The constructor uses an initialization list to set up the massTransferGunn object:
    massTransferModel(dict,sm),
    // Calls the constructor of the base class massTransferModel, passing the dictionary (dict) and simulation manager (sm).

    propsDict_(dict.subDict(typeName + "Props")),
    // Extracts a sub-dictionary named <typeName>Props from the input dictionary, which contains specific properties 
    // for the massTransferGunn model.

    multiTypes_(false),
    // Indicates whether multiple particle types are supported, default is false

    expSherwood_(propsDict_.lookupOrDefault<bool>("expSherwood",false)),
    // Determines whether to use an experimental Sherwood number. Value is read from the dictionary or defaults to false.

    interpolation_(propsDict_.lookupOrDefault<bool>("interpolation",false)),
    // Indicates whether fluid concentration is interpolated to particle positions.

    verbose_(propsDict_.lookupOrDefault<bool>("verbose",false)),
    // Enables verbose logging for debugging or detailed output.

    implicit_(propsDict_.lookupOrDefault<bool>("implicit",true)),
    // Indicates whether the mass transfer is treated implicitly. Value is read from the dictionary or defaults to true.

    // Scalar Fields:
    SPartFluidName_(propsDict_.lookupOrDefault<word>("SPartFluidName","SPartFluid")),
    // Name of the scalar field representing the particle-fluid coupling term.
    // Defaults to "SPartFluid" if not specified in the dictionary.
    SPartFluid_
    (   IOobject
        (
            SPartFluidName_,
            sm.mesh().time().timeName(),
            // The current simulation time from sm.mesh().time().timeName()
            sm.mesh(),
            // The simulation mesh (sm.mesh()).
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
            // IOobject::READ_IF_PRESENT (read if available, otherwise skip) 
            // IOobject::AUTO_WRITE (save to disk automatically).
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(1,-3,-1,0,0,0,0), 0.0) // [kg/(m3 s)]
    ),

    SPartFluidCoeffName_(propsDict_.lookupOrDefault<word>("SPartFluidCoeffName","SPartFluidCoeff")),
    // Name of the scalar field representing coefficients for the particle-fluid coupling term
    SPartFluidCoeff_
    (   IOobject
        (
            SPartFluidCoeffName_,
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,-1,0,0,0,0), 0.0) // [1/(m s)]
    ),
    ReField_
    (   IOobject
        (
            "ReField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    ShField_
    (   IOobject
        (
            "ShField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    concFieldName_(propsDict_.lookupOrDefault<word>("concFieldName","C")),
    concField_(sm.mesh().lookupObject<volScalarField> (concFieldName_)),

    // - add the temperature field to calculate the saturation conc.
    tempFieldName_(propsDict_.lookupOrDefault<word>("tempFieldName","T")),
    tempField_(sm.mesh().lookupObject<volScalarField> (tempFieldName_)),
    // -----------

    satConcFieldName_(propsDict_.lookupOrDefault<word>("satConcFieldName","Cs")),
    satConcField_(sm.mesh().lookupObject<volScalarField> (satConcFieldName_)), // set in the file transportProperties
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    maxSource_(1e30),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    coupleDEM_(propsDict_.lookupOrDefault<bool>("coupleDEM",false)),
    partMassFluxName_(propsDict_.lookupOrDefault<word>("partMassFluxName","convectiveMassFlux")),
    partMassFluxCoeffRegName_(typeName + "partMassFluxCoeff"),
    partReRegName_(typeName + "partRe"),
    partShRegName_(typeName + "partSh"),
    scaleDia_(1.),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0)))
{
    // The codes following are to:
    // register particle properties; (e.g., mass flux, Sherwood/Reynolds numbers)
    // initialize fields and limits; (e.g., source field maximum, scale factors)
    // control writing options for fields;
    // handle verbose and implicit modes appropriately.

    particleCloud_.registerParticleProperty<double**>(partMassFluxName_,1);
    // Registers properties to store specific data for each particle in the particle cloud.
    // double** indicates that the properties are two-dimensional arrays (e.g., per-particle values over time or subproperties).
    // The second argument (1) likely specifies the number of components or entries for this property.
    // Here registers the particle mass flux
    if(implicit_)
    {
        particleCloud_.registerParticleProperty<double**>(partMassFluxCoeffRegName_,1);
        // Here registers the mass flux coefficient, used in implicit schemes.
    }
    if(verbose_)
    {
        particleCloud_.registerParticleProperty<double**>(partReRegName_,1);
        particleCloud_.registerParticleProperty<double**>(partShRegName_,1);
        // Here registers the Reynolds and Sherwood numbers for each particle.
    }


    if (propsDict_.found("maxSource"))
    {
        maxSource_=readScalar(propsDict_.lookup ("maxSource"));
        Info << "limiting eulerian source field to: " << maxSource_ << endl;
    }

    if (!implicit_)
    {
        SPartFluidCoeff_.writeOpt() = IOobject::NO_WRITE;
        // If the scheme is not implicit (!implicit_), this field is not written to disk (IOobject::NO_WRITE).
    }

    if (verbose_)
    {
        ReField_.writeOpt() = IOobject::AUTO_WRITE;
        ShField_.writeOpt() = IOobject::AUTO_WRITE;
        ReField_.write();
        ShField_.write();
        if (expSherwood_)
        {
          FatalError <<"Cannot read and create ShField at the same time!\n" << abort(FatalError);
        }
    }

    if (propsDict_.found("scale") && typeCG_.size()==1)
    {
        // if "scale" is specified and there's only one single type, use "scale"
        scaleDia_=scalar(readScalar(propsDict_.lookup("scale")));
        typeCG_[0] = scaleDia_;
    }
    else if (typeCG_.size()>1)
    {
        multiTypes_ = true;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

massTransferJang::~massTransferJang()
{
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void massTransferJang::calcMassContribution()
{
    double**& partMassFlux_ = particleCloud_.getParticlePropertyRef<double**>(partMassFluxName_);

    // reset Scalar field
    SPartFluid_.primitiveFieldRef() = 0.0;
    
    // Defines the fluid's viscosity, accounting for turbulence.
    #ifdef compre
       const volScalarField mufField = particleCloud_.turbulence().mu();
    #else
       const volScalarField mufField = particleCloud_.turbulence().nu()*rho_;
    #endif
    
    // Retrieves a reference to the diffusion coefficient field
    const volScalarField& D0Field_ = D0Field();

    if (typeCG_.size()>1 || typeCG_[0] > 1)
    {
        Info << "massTransferJang using scale = " << typeCG_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "massTransferJang using scale from liggghts cg = " << scaleDia_ << endl;
        // it is the coarse-graining factor in liggghts, i.e. cg > 1 means particle grouping
    }

    // calc La based heat flux
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    scalar Tfluid(0); //added
    scalar Cfluid(0);
    scalar Csfluid(0);
    label cellI=0;
    vector Us(0,0,0);
    scalar ds(0);
    // Solid diameter
    scalar ds_scaled(0);
    // Solid diameter scaled by cg
    scalar scaleDia3 = scaleDia_*scaleDia_*scaleDia_;
    scalar muf(0);
    scalar rhof(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Sc(0);
    // Sc: Schmidt number
    scalar Shp(0);

    scalar cg = scaleDia_;
    label partType = 1;

    interpolationCellPoint<scalar> voidfractionInterpolator_(voidfraction_);
    interpolationCellPoint<vector> UInterpolator_(U_);
    interpolationCellPoint<scalar> CInterpolator_(concField_);
    interpolationCellPoint<scalar> TInterpolator_(tempField_);

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                if(interpolation_)
                {
                    vector position = particleCloud_.position(index);
                    voidfraction = voidfractionInterpolator_.interpolate(position,cellI);
                    Ufluid = UInterpolator_.interpolate(position,cellI);
                    Cfluid = CInterpolator_.interpolate(position,cellI);
                    Tfluid = TInterpolator_.interpolate(position,cellI);  //added
                }
                else
                {
                    voidfraction = voidfraction_[cellI];
                    Ufluid = U_[cellI];
                    Cfluid = concField_[cellI];
                    Tfluid = tempField_[cellI];  //added
                }

                if (voidfraction < 0.01)
                    voidfraction = 0.01;

                if (multiTypes_)
                {
                    partType = particleCloud_.particleType(index);
                    cg = typeCG_[partType - 1];
                    scaleDia3 = cg*cg*cg;
                }

                // calc relative velocity
                Us = particleCloud_.velocity(index);
                magUr = mag(Ufluid - Us);
                ds = 2.*particleCloud_.radius(index);
                ds_scaled = ds/cg;
                muf = mufField[cellI];
                rhof = rho_[cellI];
                Csfluid = 1.34 + 0.00257 * Tfluid * 67; 
                // The unit here is wt%, so it needs to be times by c to convert into kg/m3
                // c is rho_fluid/100, so here is 67 generally

                Rep = ds_scaled * magUr * voidfraction * rhof/ muf;

                scalar D0 = D0Field_[cellI]; //mass diffusivity, m2/s

                Sc = max(SMALL, muf / (rhof*D0));

                Shp = Sherwood(voidfraction, Rep, Sc);

                scalar h = D0 * Shp / ds_scaled;  //particle-to-fluid mass transfer coefficient, Kpf in the paper, m/s
                scalar As = ds_scaled * ds_scaled * M_PI; // surface area of sphere

                // calc convective heat flux [W]
                massFlux(index, h, As, Cfluid, Csfluid, scaleDia3);

                if(verbose_)
                {
                    double**& partRe_ = particleCloud_.getParticlePropertyRef<double**>(partReRegName_);
                    double**& partSh_ = particleCloud_.getParticlePropertyRef<double**>(partShRegName_);
                    partRe_[index][0] = Rep;
                    partSh_[index][0] = Shp;
                }

                if(verbose_ && index >=0 && index <2)
                {
                    Pout << "partMassFlux = " << partMassFlux_[index][0] << endl;
                    Pout << "magUr = " << magUr << endl;
                    Pout << "D0 = " << D0 << endl;
                    Pout << "rho = " << rhof << endl;
                    Pout << "h = " << h << endl;
                    Pout << "ds = " << ds << endl;
                    Pout << "ds_scaled = " << ds_scaled << endl;
                    Pout << "As = " << As << endl;
                    Pout << "muf = " << muf << endl;
                    Pout << "Rep = " << Rep << endl;
                    Pout << "Sc = " << Sc << endl;
                    Pout << "Shp = " << Shp << endl;
                    Pout << "voidfraction = " << voidfraction << endl;
                    Pout << "Cfluid = " << Cfluid << endl;
                    Pout << "Csfluid = " << Csfluid << endl;
                }
            }
    }

    particleCloud_.averagingM().setScalarSum
    (
        SPartFluid_,
        partMassFlux_,
        particleCloud_.particleWeights(),
        NULL
    );

    SPartFluid_.primitiveFieldRef() /= -SPartFluid_.mesh().V();

    if(implicit_)
    // This code seems to be part of the implicit solution process where mass flux coefficients between particles
    // and the fluid are computed and averaged.
    // The SPartFluidCoeff_ field is initialized and then updated with the average of particle-to-fluid mass 
    // flux coefficients, weighted by particle weights.
    // Finally, the result is normalized by the mesh volume, likely to adjust the coefficients for the simulation's 
    // spatial resolution.
    {
        double**& partMassFluxCoeff_ = particleCloud_.getParticlePropertyRef<double**>(partMassFluxCoeffRegName_);
        SPartFluidCoeff_.primitiveFieldRef() = 0.0; //initialize the field to 0

        particleCloud_.averagingM().setScalarSum 
        // Computes an averaged scalar sum, which likely sums up the mass flux coefficients across the particles
        // and weighted by their respective weights
        (
            SPartFluidCoeff_,
            partMassFluxCoeff_,
            particleCloud_.particleWeights(),
            NULL
        );

        SPartFluidCoeff_.primitiveFieldRef() /= -SPartFluidCoeff_.mesh().V();
        // Normalized the field by dividing it by the volume of the mesh
    }

    if(verbose_)
    // The code checks if verbose output is enabled. If so, it calculates the average Reynolds number (ReField_) 
    // and Sherwood number (ShField_) for the particles in the simulation.
    {
        double**& partRe_ = particleCloud_.getParticlePropertyRef<double**>(partReRegName_);
        // This line retrieves a reference to the particle property partRe_ (Reynolds number) from the particleCloud_ 
        // object. The property is expected to be a 2D array (double**), which likely represents the Reynolds number 
        // for each particle in the simulation.
        double**& partSh_ = particleCloud_.getParticlePropertyRef<double**>(partShRegName_);
        // Similarly, this is for the sherwood number
        ReField_.primitiveFieldRef() = 0.0;
        ShField_.primitiveFieldRef() = 0.0;
        // Initiallized the fields
        particleCloud_.averagingM().resetWeightFields();
        // Reset the weight fields for the particle averaging process
        particleCloud_.averagingM().setScalarAverage
        // This line computes the scalar average of the Reynolds number across the particles and stores it in ReField_.
        (
            ReField_,
            partRe_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            ShField_,
            partSh_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );
    }

    // limit source term in explicit treatment
    if(!implicit_)
    {
        forAll(SPartFluid_,cellI)
        {
            scalar EuFieldInCell = SPartFluid_[cellI];

            if(mag(EuFieldInCell) > maxSource_ )
            {
                 Pout << "limiting source term\n"  << endl  ;
                 SPartFluid_[cellI] = sign(EuFieldInCell) * maxSource_;
            }
        }
    }

    SPartFluid_.correctBoundaryConditions();

    volScalarField minParticleWeights = particleCloud_.averagingM().UsWeightField();
    Info << "Minimum Particle Weight " << gMin(minParticleWeights) << endl;
    Info << "Minimum Fluid Concentration: " << gMin(concField_) << endl;
    Info << "Maximum Fluid Concentration: " << gMax(concField_) << endl;
}

void massTransferJang::addMassContribution(volScalarField& Sm) const
{
    Sm += SPartFluid_;
}

void massTransferJang::addMassCoefficient(volScalarField& Smi) const
{
    if(implicit_)
    {
        Smi += SPartFluidCoeff_;
    }
}

scalar massTransferJang::Sherwood(scalar voidfraction, scalar Rep, scalar Sc) const
{
    return (7 - 10 * voidfraction + 5 * voidfraction * voidfraction) *
                        (1 + 0.7 * Foam::pow(Rep,0.2) * Foam::pow(Sc,0.33)) +
                        (1.33 - 2.4 * voidfraction + 1.2 * voidfraction * voidfraction) *
                        Foam::pow(Rep,0.7) * Foam::pow(Sc,0.33);
}

void massTransferJang::massFlux(label index, scalar h, scalar As, scalar Cfluid, scalar Csfluid, scalar cg3)
// Calculates and updates the mass flux between a particle and the surrounding fluid, considering the 
// transfer coefficient, surface area, concentrations, and correction factors.
{
    scalar hAs = h * As * cg3;
    // h: particle to fluid mass transfer coefficient, m/s
    // As: particle surface area, m2
    double**& partMassFlux_ = particleCloud_.getParticlePropertyRef<double**>(partMassFluxName_);
    // Retrieves a reference to the partMassFlux_ array from the particleCloud_ object.

    if (particleCloud_.getParticleEffVolFactors())
    // Checks if effective volume factors are being used for the particles
    {
        scalar effVolFac = particleCloud_.particleEffVolFactor(index);
        hAs *= effVolFac;
    }

    partMassFlux_[index][0] = - hAs * Csfluid; 
    // Csfluid: satuarted concentration in fluid
    // Here the initial mass flux for the particle is calculated.
    // The negative sign might indicate that this is the mass flux out of the particle. 

    if(!implicit_)
    // If implicit: Adds hAs * Cfluid to the mass flux. This introduces a term that depends on the 
    // current concentration in the fluid.
    {
        partMassFlux_[index][0] += hAs * Cfluid;
    }
    else
    {
        double**& partMassFluxCoeff_ = particleCloud_.getParticlePropertyRef<double**>(partMassFluxCoeffRegName_);
        partMassFluxCoeff_[index][0] = hAs;
    }
}

void massTransferJang::giveData()
{
    Info << "total convective particle-fluid mass flux [kg/s] (Eulerian) = " << gSum(SPartFluid_*1.0*SPartFluid_.mesh().V()) << endl;

    double**& partMassFlux_ = particleCloud_.getParticlePropertyRef<double**>(partMassFluxName_);
    particleCloud_.dataExchangeM().giveData(partMassFluxName_,"scalar-atom", partMassFlux_);
}

void massTransferJang::postFlow()
// Post-process the mass flux (partMassFlux_) for each particle after fluid-related calculations are complete.
// Transfers the computed flux data to the DEM solver using the giveData() function.
{
    // send mass flux to DEM
    if (coupleDEM_)
    // The mass flux is processed only if DEM coupling is enabled.
    {
        if(implicit_)
        {
            label cellI;
            // The index of the fluid cell corresponding to the current particle.

            scalar Cfluid(0.0);
            // The fluid concentration near the particle, initialized to 0.

            interpolationCellPoint<scalar> CInterpolator_(concField_);
            // An interpolator object is created for the fluid concentration field (concField_), allowing 
            // precise retrieval of concentration values at particle positions.

            double**& partMassFlux_ = particleCloud_.getParticlePropertyRef<double**>(partMassFluxName_);
            double**& partMassFluxCoeff_ = particleCloud_.getParticlePropertyRef<double**>(partMassFluxCoeffRegName_);
            // Retrieve references to the particle mass flux array (partMassFlux_) and mass flux coefficient array (partMassFluxCoeff_).

            for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
            // Particle loop
            {
                cellI = particleCloud_.cellIDs()[index][0];
                // Retrieves the fluid cell ID (cellI) corresponding to the particle at the current index.
                if(cellI >= 0)
                // Ensures the particle is associated with a valid fluid cell (negative IDs likely represent invalid or boundary conditions).
                {
                    if(interpolation_)
                    // If interpolation is enabled: 
                    // use the particle’s position to interpolate the concentration from the fluid field (CInterpolator_),
                    // otherwise use the concentration value directly from the fluid cell (concField_[cellI]).
                    {
                        vector position = particleCloud_.position(index);
                        Cfluid = CInterpolator_.interpolate(position,cellI);
                    }
                    else
                    {
                        Cfluid = concField_[cellI];
                    }

                    partMassFlux_[index][0] += Cfluid * partMassFluxCoeff_[index][0];
                    // Updates the mass flux for the current particle:
                    // The mass flux is increased by the fluid concentration and the 
                    // The coefficient computed during earlier steps.

                    if(verbose_ && index >=0 && index <2)
                    // Debugging Output (Optional):
                    {
                        Pout << "partMassFlux = " << partMassFlux_[index][0] << endl;
                        Pout << "Cfluid = " << Cfluid << endl;
                    }
                }
            }
        }

        giveData();
        // A call to giveData() transfers the updated mass flux data to the DEM solver.
        // This step ensures that the particle mass flux data is available for DEM calculations in the next simulation stage.
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

