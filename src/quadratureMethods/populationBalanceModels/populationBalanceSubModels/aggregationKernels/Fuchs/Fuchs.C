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
    df_(readScalar(dict.lookup("df"))),
    lambdaf_(dict.lookup("lambdaf")),
    abscissa0_(dict.lookup("abscissa0"))
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
//     dimensionedScalar abscissa0 
//     (
//         "abscissa0",
//         abscissa.dimensions(),
//         abscissa0_
//     );
//     
    return 2.0*max(abscissa, abscissa0_)*pow(max(abscissa, abscissa0_)/abscissa0_,1.0/df_);
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::velocity
(
    const volScalarField& abscissa
) const
{      
    const fluidThermo& flThermo = abscissa.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
    
    return sqrt(8.0*Foam::constant::physicoChemical::k*flThermo.T()/(Foam::constant::mathematical::pi*flThermo.rho()*pow3(6.0/Foam::constant::mathematical::pi*max(abscissa, abscissa0_))));
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
    -pow(dc(abscissa)*dc(abscissa)+l(abscissa)*l(abscissa),3.0/2.0))-dc(abscissa);
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::aggregationKernels::Fuchs::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{ 
    
//   tmp<volScalarField> aggSource(abscissa1);
    
//     tmp<volScalarField> aggSource
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "aggregationSource",
//                 abscissa1.mesh().time().timeName(),
//                 abscissa1.mesh(),
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE,
//                 false
//             ),
//             abscissa1.mesh(),
//             dimensionedScalar("aggSource", dimensionSet(0,-3,-1,0,0), 0.0)
//         )
//     );   
//     
//     const fluidThermo& flThermo = abscissa1.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
//     
//     forAll(abscissa1,cellI)
//     {
//         if(abscissa1[cellI]>abscissa0_.value() && abscissa2[cellI]>abscissa0_.value())
//         {
//             flThermo.T()[cellI]=1.0;
//             D(abscissa1)[cellI]=1.0;
//             
//             aggSource = 2.0*Foam::constant::mathematical::pi*(D(abscissa1)+D(abscissa2))*(dc(abscissa1)+dc(abscissa2))*pow((dc(abscissa1)+dc(abscissa2))/(dc(abscissa1)+dc(abscissa2)+2.0*sqrt(g(abscissa1)*g(abscissa1)+g(abscissa2)*g(abscissa2)))+8.0*(D(abscissa1)+D(abscissa2))/((dc(abscissa1)+dc(abscissa2))*sqrt(velocity(abscissa1)*velocity(abscissa1)+velocity(abscissa2)*velocity(abscissa2))),-1.0);
//         }
//     }
//     
    
    tmp<Foam::volScalarField> aggregationKernel(2.0*Foam::constant::mathematical::pi*(D(abscissa1)+D(abscissa2))*(dc(abscissa1)+dc(abscissa2))*pow((dc(abscissa1)+dc(abscissa2))/(dc(abscissa1)+dc(abscissa2)+2.0*sqrt(g(abscissa1)*g(abscissa1)+g(abscissa2)*g(abscissa2)))+8.0*(D(abscissa1)+D(abscissa2))/((dc(abscissa1)+dc(abscissa2))*sqrt(velocity(abscissa1)*velocity(abscissa1)+velocity(abscissa2)*velocity(abscissa2))),-1.0));
    
    //Info << "aggregation kernel : " << aggregationKernel.ref().dimensions() << endl;
    
    return aggregationKernel;
}

// ************************************************************************* //
