/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alberto Passalacqua
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

#include "growthReduction.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace convectionModels
{
    defineTypeNameAndDebug(growthReduction, 0);

    addToRunTimeSelectionTable
    (
        convectionModel,
        growthReduction,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::convectionModels::growthReduction::growthReduction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    convectionModel(dict, mesh),
    VCarbon_("CarbonVolume", pow3(dimLength), 1.1e-2),
    Xi_("concentration hydrogen sites", inv(sqr(dimLength)), 17)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::convectionModels::growthReduction::~growthReduction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField Foam::populationBalanceSubModels::convectionModels::growthReduction::Arrhenius
(
    const scalar& A,
    const scalar& n,
    const scalar& E,
    const label& nReactifs //reactif number
)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    volScalarField Ea
    (
        IOobject
        (
            "Ea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("activation energy", T.dimensions(), E*1000)
    );
    
    volScalarField Arrhenius = A*pow(T,n)*exp(-Ea/(8.314*T)); //careful, the dimensions of E in kJ and A divided by 1e6 to have the result by m^3 instead of cm^3
    
    if (nReactifs==2)
    {
        Arrhenius.dimensions().reset(Foam::dimensionSet(0,3,-1,0,-1,0,0));//dimensions 
    }
    else 
    {
        Arrhenius.dimensions().reset(Foam::dimensionSet(0,0,-1,0,0,0,0));
    }
    
    return Arrhenius;
}

Foam::scalar Foam::populationBalanceSubModels::convectionModels::growthReduction::Arrhenius
(
    const scalar& A,
    const scalar& n,
    const scalar& E,
    const scalar& T
)
{
    return A*pow(T,n)*exp(-E*1000/(8.314*T)); 
}

Foam::volScalarField Foam::populationBalanceSubModels::convectionModels::growthReduction::r()
{
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    const volScalarField& concentration_H(mesh_.lookupObject<volScalarField>("H"));
    
    const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& concentration_H2O(mesh_.lookupObject<volScalarField>("H2O"));
    
    const volScalarField& concentration_H2(mesh_.lookupObject<volScalarField>("H2"));
    
    const volScalarField& concentration_C2H2(mesh_.lookupObject<volScalarField>("C2H2"));
    
    /*volScalarField r = Arrhenius(1.00e2,1.80,68.42,2)*concentration_H*rho/0.001+Arrhenius(6.72e-5,3.33,6.09,2)*concentration_OH*rho/0.005+Arrhenius(1.13e16,-0.06,476.05,1)/(Arrhenius(8.68e-2,2.36,25.46,2)*concentration_H2*rho/0.002+Arrhenius(6.44e-7,3.79,27.96,2)*concentration_H2O*rho/0.018+Arrhenius(4.17e6,0.15,0.00,2)*concentration_H*rho/0.001+Arrhenius(2.52e3,1.10,17.13,2)*concentration_C2H2*rho/0.026);*/
    
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    volScalarField r
    (
        IOobject
        (
            "r",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    
    forAll(concentration_H, cellI)
    {
        if (concentration_H[cellI]!=0) 
        {
            r[cellI] = Arrhenius(1e2,1.80,68.42,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001+Arrhenius(8.68e-5,2.36,25.46,T[cellI])*concentration_OH[cellI]*rho[cellI]/0.005+Arrhenius(1.13e16,-0.06,476.05,T[cellI])/(Arrhenius(8.68e-2,2.36,25.46,T[cellI])*concentration_H2[cellI]*rho[cellI]/0.002+Arrhenius(6.44e-7,3.79,27.96,T[cellI])*concentration_H2O[cellI]*rho[cellI]/0.018+Arrhenius(4.17e6,0.15,0.00,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001+Arrhenius(2.52e3,1.10,17.13,T[cellI])*concentration_C2H2[cellI]*rho[cellI]/0.026);
        }
    }
    
    //Info << r << endl;
    
    return r;
}

Foam::volScalarField Foam::populationBalanceSubModels::convectionModels::growthReduction
::cRed()
{
    volScalarField kReaction = Arrhenius(2.20e6,0.0,31.38,2); //reaction : Soot* + O2 -> Soot-H + 2CO + arrhenius law
    
    volScalarField kOH = Arrhenius(1.0e10,0.0,8.4,2); //reaction Soot-H + OH -> Soot-H + C0 + reaction probability 0.13 ???
    
    volScalarField concentration_O2(mesh_.lookupObject<volScalarField>("O2"));
    
    volScalarField concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    volScalarField cRed = (2.0*r()/(r()+1.0)*kReaction*concentration_O2/0.032+kOH*concentration_OH/0.017)*rho*Xi_*VCarbon_*pow(36.0*Foam::constant::mathematical::pi,1.0/3.0);

    /*concentration_OH -= r()/(r()+1.0)*concentration_OH*Xi_*kOH*pow(36.0*Foam::constant::mathematical::pi,1.0/3.0)*pow(abscissa,2.0/3.0)*moment/Foam::constant::physicoChemical::Na*mesh_.time().deltaT();
    
    concentration_O2 -= r()/(r()+1.0)*concentration_O2*Xi_*kReaction*pow(36.0*Foam::constant::mathematical::pi,1.0/3.0)*pow(abscissa,2.0/3.0)*moment/Foam::constant::physicoChemical::Na*mesh_.time().deltaT();*/
    
    return cRed;
}

Foam::volScalarField Foam::populationBalanceSubModels::convectionModels::growthReduction
::cGro()
{
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    volScalarField concentration_C2H2(mesh_.lookupObject<volScalarField>("C2H2"));
    
    volScalarField cGro = 2.0*VCarbon_*r()/(r()+1.0)*concentration_C2H2*Xi_*Arrhenius(2.52e3,1.10,17.13,2)*rho/0.026*pow(36.0*Foam::constant::mathematical::pi,1.0/3.0);
    
    /*concentration_C2H2 -= r()/(r()+1.0)*concentration_C2H2*Xi_*Arrhenius(2.52e3,1.10,17.13,2)*pow(36.0*Foam::constant::mathematical::pi,1.0/3.0)*pow(abscissa,2.0/3.0)*moment/Foam::constant::physicoChemical::Na*mesh_.time().deltaT();*/
    
    return cGro;
}

Foam::scalar
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic
(
    const scalar& cRed,
    const scalar& cGro    
) 
{
    return max(pow3((cRed-cGro)/6.0*mesh_.time().deltaT().value()),0.0);
}

Foam::scalar
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic
(
    const scalar& abscissa, 
    const scalar& cRed,
    const scalar& cGro
) 
{
    return pow3(pow(abscissa,1.0/3.0) + (cGro-cRed)/6.0*mesh_.time().deltaT().value());
}
