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
    VCarbon_(1.9e-30),
    Xi_(readScalar(dict.lookup("Xi")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::convectionModels::growthReduction::~growthReduction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::populationBalanceSubModels::convectionModels::growthReduction::Arrhenius
(
    const scalar& A,
    const scalar& n,
    const scalar& E,
    const scalar& T
) const
{
    return A*pow(T,n)*exp(-E/(0.008314*T)); //careful, the dimensions of E in kJ and A multiply by 1e6 to have the result by m^3 instead of cm^3
}

Foam::scalar Foam::populationBalanceSubModels::convectionModels::growthReduction
::cRed(const label& cellI)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    scalar kReaction = Arrhenius(2.20e18,0.0,31.38,T[cellI]); //reaction : Soot* + O2 -> Soot-H + 2CO + arrhenius law
    
    const volScalarField& concentration_O2(mesh_.lookupObject<volScalarField>("O2"));
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    scalar cRed = 2*VCarbon_*kReaction*Xi_*concentration_O2[cellI]*rho[cellI]/0.016;

    //Info << "cRed : " << cRed << endl;
    //Info << "kReaction : " << kReaction << endl;
    return cRed;
}

Foam::scalar Foam::populationBalanceSubModels::convectionModels::growthReduction
::cGro(const label& cellI)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    const volScalarField& concentration_H(mesh_.lookupObject<volScalarField>("H"));
    
    const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& concentration_H2O(mesh_.lookupObject<volScalarField>("H2O"));
    
    const volScalarField& concentration_H2(mesh_.lookupObject<volScalarField>("H2"));
    
    const volScalarField& concentration_C2H2(mesh_.lookupObject<volScalarField>("C2H2"));
    
    scalar r = 0.0;
    
    if (concentration_H[cellI]!=0)
    {
        r = Arrhenius(1e14,1.80,68.42,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001+Arrhenius(8.68e10,2.36,25.46,T[cellI])*concentration_OH[cellI]*rho[cellI]/0.003+Arrhenius(1.13e25,-0.06,476.05,T[cellI])/(Arrhenius(8.68e10,2.36,25.46,T[cellI])*concentration_H2[cellI]*rho[cellI]/0.002+Arrhenius(6.44e9,3.79,27.96,T[cellI])*concentration_H2O[cellI]*rho[cellI]/0.01+Arrhenius(4.17e19,0.15,0.00,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001+Arrhenius(2.52e15,1.10,17.13,T[cellI])*concentration_C2H2[cellI]*rho[cellI]/0.014);
    }
        
    scalar cGro = 2*VCarbon_*r/(r+1)*concentration_C2H2[cellI]*Arrhenius(2.52e15,1.10,17.13,T[cellI])*rho[cellI]/0.014;
    
    //Info << "cGro : " << cGro << endl;
    
    return cGro;
}

Foam::scalar
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic(const label& cellI) 
{
    return cGro(cellI)/3.0*mesh_.time().deltaT().value();
    //return (cGro(cellI)-cRed(cellI))/3.0*mesh_.time().deltaT().value();
}

Foam::scalar
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic(const scalar& abscissa, const label& cellI) 
{
    return (abscissa - cGro(cellI)/3.0*mesh_.time().deltaT().value());
    //return (abscissa - (cGro(cellI)-cRed(cellI))/3.0*mesh_.time().deltaT().value());
}

// ************************************************************************* //
