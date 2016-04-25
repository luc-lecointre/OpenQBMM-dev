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
    return A*pow(T,n)*exp(-E/(8.314*T));
}

Foam::scalar Foam::populationBalanceSubModels::convectionModels::growthReduction
::cRed(const label& cellI)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    scalar kReaction = Arrhenius(2.20e12,0.0,31.38,T[cellI]); //reaction : Soot* + O2 -> Soot-H + 2CO + arrhenius law
    
    const volScalarField& concentration_O2(mesh_.lookupObject<volScalarField>("O2"));
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    return 2*VCarbon_*kReaction*Xi_*Foam::constant::physicoChemical::NA.value()*concentration_O2[cellI]*rho[cellI]/0.016;
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
    
    scalar r = Arrhenius(1e8,1.80,68.42,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001
        +Arrhenius(8.68e4,2.36,25.46,T[cellI])*concentration_OH[cellI]*rho[cellI]/0.003
        +Arrhenius(1.13e16,-0.06,476.05,T[cellI])/
        (Arrhenius(8.68e4,2.36,25.46,T[cellI])*concentration_H2[cellI]*rho[cellI]/0.002
        +Arrhenius(6.44e-1,3.79,27.96,T[cellI])*concentration_H2O[cellI]*rho[cellI]/0.01
        +Arrhenius(4.17e13,0.15,0.00,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001
        +Arrhenius(2.52e9,1.10,17.13,T[cellI])*concentration_C2H2[cellI]*rho[cellI]/0.014);
        
    return 2*VCarbon_*r/(r+1)*Foam::constant::physicoChemical::NA.value()*concentration_C2H2[cellI]*rho[cellI]/0.014;
}

Foam::scalar
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic(const label& cellI) 
{
    return (cRed(cellI)-cGro(cellI))*mesh_.time().deltaT().value();
}

Foam::scalar
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic(const scalar& abscissa, const label& cellI) 
{
    return (abscissa - (cRed(cellI)-cGro(cellI))*mesh_.time().deltaT().value());
}

// ************************************************************************* //
