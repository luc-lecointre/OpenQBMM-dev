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
    nCarbonSoot_(readScalar(dict.lookup("nCarbonSoot"))),
    rhoSoot_(
        "rhoSoot",
        Foam::dimensionSet(1,-3,0,0,0,0,0),
        1800.0
    ),
    abscissa0_(MCarbon_*nCarbonSoot_/(rhoSoot_*Foam::constant::physicoChemical::NA))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::Fuchs
::~Fuchs()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField Foam::populationBalanceSubModels::aggregationKernels::Fuchs::Kn
(
    const volScalarField& abscissa
) const
{
    const fluidThermo& flThermo = abscissa.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
    
    volScalarField lambdaf = flThermo.mu()/flThermo.p()*sqrt(Foam::constant::mathematical::pi*Foam::constant::physicoChemical::k*flThermo.T()/(2.0*rhoSoot_*abscissa)); 
    
    return 2.0*lambdaf/dc(abscissa);
}

Foam::volScalarField Foam::populationBalanceSubModels::aggregationKernels::Fuchs::D
(
    const volScalarField& abscissa
) const
{  
    const fluidThermo& flThermo = abscissa.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
    
    volScalarField Kna = Kn(abscissa);
    
    return Foam::constant::physicoChemical::k*flThermo.T()/(3.0*Foam::constant::mathematical::pi*flThermo.mu()*dc(abscissa))*(5.0+4.0*Kna+6.0*sqr(Kna)+8.0*pow3(Kna))/(5.0-Kna+(8.0+Foam::constant::mathematical::pi)*sqr(Kna));
}

Foam::volScalarField Foam::populationBalanceSubModels::aggregationKernels::Fuchs::dc
(
    const volScalarField& abscissa
) const
{       
    return 2.0*pow(3.0/(4.0*Foam::constant::mathematical::pi)*abscissa,1.0/3.0)*pow(abscissa/abscissa0_,1.0/df_);
}

Foam::volScalarField Foam::populationBalanceSubModels::aggregationKernels::Fuchs::velocity
(
    const volScalarField& abscissa
) const
{      
    const fluidThermo& flThermo = abscissa.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
    
    return sqrt(8.0*Foam::constant::physicoChemical::k*flThermo.T()/(Foam::constant::mathematical::pi*rhoSoot_*abscissa));
}

Foam::volScalarField Foam::populationBalanceSubModels::aggregationKernels::Fuchs::l
(
    const volScalarField& abscissa
) const
{  
    return 8.0*D(abscissa)/(Foam::constant::mathematical::pi*velocity(abscissa));
}

Foam::volScalarField Foam::populationBalanceSubModels::aggregationKernels::Fuchs::g
(
    const volScalarField& abscissa
) const
{   
    volScalarField dca = dc(abscissa);
    volScalarField la = l(abscissa);
    
    return 1.0/(3.0*dca*la)*(pow3(dca+la)
    -pow(sqr(dca)+sqr(la),3.0/2.0))-dca;
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{   
    
    volScalarField a1=max(abscissa1*1e-27,abscissa0_);
    volScalarField a2=max(abscissa2*1e-27,abscissa0_);
    
    volScalarField  dc1 = dc(a1);
    volScalarField  dc2 = dc(a2);
    
    volScalarField Da = D(a1) + D(a2);
    
    volScalarField A = (dc1+dc2+2.0*sqrt(sqr(g(a1))+sqr(g(a2))));
    volScalarField B = sqr(dc1+dc2)*sqrt(sqr(velocity(a1))+sqr(velocity(a2)));
    
    //Foam::tmp<Foam::volScalarField>  aggregationKernel = 2.0*Foam::constant::mathematical::pi*(D1+D2)*(dc1+dc2)*pow((dc1+dc2)/(dc1+dc2+2.0*sqrt(sqr(g(a1))+sqr(g(a2))))+8.0*(D1+D2)/((dc1+dc2)*sqrt(sqr(velocity(a1))+sqr(velocity(a2)))),-1.0);
    
    Foam::tmp<Foam::volScalarField>  aggregationKernel = 2.0*Foam::constant::mathematical::pi*Da*A*B/(8.0*Da*A+B);
    
    //Info << "aggregationKernel " << aggregationKernel.ref() << endl;
    
    return aggregationKernel;
}

// ************************************************************************* //
