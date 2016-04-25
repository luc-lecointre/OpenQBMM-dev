/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.
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

#include "Fuchs.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(Fuchs, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        Fuchs,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::Fuchs
::Fuchs
(
    const dictionary& dict
)
:
    aggregationKernel(dict),
    MCarbon_
    (
        "carbonMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.012
    ),
    df_(readScalar(dict.lookup("df"))),
    lambdaf_(dict.lookup("lambdaf")),
    nCarbonSoot_(readScalar(dict.lookup("nCarbonSoot"))),
    rhoSoot_(
        "rhoSoot",
        Foam::dimensionSet(1,-3,0,0,0,0,0),
        1800.0
    ),
    abscissa0_(6.0/Foam::constant::mathematical::pi*pow(MCarbon_*nCarbonSoot_/(rhoSoot_*Foam::constant::physicoChemical::NA),1.0/3.0))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::Fuchs
::~Fuchs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::Kn
(
    const volScalarField& abscissa
) const
{
    /*const fluidThermo& flThermo = abscissa.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
    
    tmp<volScalarField> lambdaf = flThermo.mu()/flThermo.p()*sqrt(Foam::constant::mathematical::pi*Foam::constant::physicoChemical::k*flThermo.T()/(2*rhoSoot_*Foam::constant::mathematical::pi/6.0*pow3(abscissa0_)));*/
    
    
    return 2.0*lambdaf_/dc(abscissa);
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::D
(
    const volScalarField& abscissa
) const
{  
    const fluidThermo& flThermo = abscissa.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
    
    return Foam::constant::physicoChemical::k*flThermo.T()/(3.0*Foam::constant::mathematical::pi*flThermo.mu()*dc(abscissa))*(5.0+4.0*Kn(abscissa)+6.0*Kn(abscissa)*Kn(abscissa)+8.0*pow3(Kn(abscissa)))/(5.0-Kn(abscissa)+(8.0+Foam::constant::mathematical::pi)*Kn(abscissa)*Kn(abscissa));
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::dc
(
    const volScalarField& abscissa
) const
{       
    return abscissa*pow(abscissa/abscissa0_,3/df_);
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::velocity
(
    const volScalarField& abscissa
) const
{      
    const fluidThermo& flThermo = abscissa.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
    
    return sqrt(8.0*Foam::constant::physicoChemical::k*flThermo.T()/(sqr(Foam::constant::mathematical::pi)*rhoSoot_/6.0*pow3(abscissa)));
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::l
(
    const volScalarField& abscissa
) const
{  
    return 8.0*D(abscissa)/(Foam::constant::mathematical::pi*velocity(abscissa));
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::g
(
    const volScalarField& abscissa
) const
{   
    return 1.0/(3.0*dc(abscissa)*l(abscissa))*(pow3(dc(abscissa)+l(abscissa))
    -pow(sqr(dc(abscissa))+sqr(l(abscissa)),3.0/2.0))-dc(abscissa);
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{   
    
    /*tmp<volScalarField> aggregationKernel
    (
        new volScalarField
        (
            IOobject
            (
                "aggregationKernel",
                abscissa1.mesh().time().timeName(),
                abscissa1.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            abscissa1.mesh(),
            dimensionedScalar("zero", dimless , 0.0)
        )
    );
    
    volScalarField& kernel = aggregationKernel.ref();
    
    forAll(abscissa1,cellI)
    {
        if ( abscissa1[cellI]==0 )
        {
            kernel[cellI] = 0;
        }
        else
        {
            kernel[cellI] = 2.0*Foam::constant::mathematical::pi*(D(abscissa1)+D(abscissa2))*(dc(abscissa1)+dc(abscissa2))*pow((dc(abscissa1)+dc(abscissa2))/(dc(abscissa1)+dc(abscissa2)+2.0*sqrt(g(abscissa1)*g(abscissa1)+g(abscissa2)*g(abscissa2)))+8.0*(D(abscissa1)+D(abscissa2))/((dc(abscissa1)+dc(abscissa2))*sqrt(velocity(abscissa1)*velocity(abscissa1)+velocity(abscissa2)*velocity(abscissa2))),-1.0);
        }
    }
    kernel.dimensions().reset(dimensionSet(0,-3,1,0,0,0,0));*/
    
    volScalarField a1=max(abscissa1,abscissa0_);
    volScalarField a2=max(abscissa2,abscissa0_);
    
    tmp<volScalarField> aggregationKernel = 2.0*Foam::constant::mathematical::pi*(D(a1)+D(a2))*(dc(a1)+dc(a2))*pow((dc(a1)+dc(a2))/(dc(a1)+dc(a2)+2.0*sqrt(g(a1)*g(a1)+g(a2)*g(a2)))+8.0*(D(a1)+D(a2))/((dc(a1)+dc(a2))*sqrt(velocity(a1)*velocity(a1)+velocity(a2)*velocity(a2))),-1.0);
    
    return aggregationKernel;
}

// ************************************************************************* //
