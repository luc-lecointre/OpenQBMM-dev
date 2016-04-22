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

#include "reduction.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace oxidationModels
{
    defineTypeNameAndDebug(reduction, 0);

    addToRunTimeSelectionTable
    (
        oxidationModel,
        reduction,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::oxidationModels::reduction::reduction
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    oxidationModel(dict, mesh)
    //cRed_(readScalar(dict.lookup("cRed")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::oxidationModels::reduction::~reduction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::populationBalanceSubModels::oxidationModels::reduction
::cRed(const label& cellI)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    scalar dCarbon = 2.365e-10;
    
    //scalar kReaction = ArrheniusReactionRate(2.20e12,0.0,31.38);
    
    scalar kReaction = 2.20e12*exp(-31.38/(8.314*T[cellI])); //reaction : Soot* + O2 -> Soot-H + 2CO + arrhenius law
    
    return 2*dCarbon*kReaction*1.7e19;
    
    
}

Foam::scalar
Foam::populationBalanceSubModels::oxidationModels::reduction
::characteristic(const label& cellI) 
{
    const volScalarField& concentration_O2(mesh_.lookupObject<volScalarField>("O2"));
    
    //const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    return (cRed(cellI)*concentration_O2[cellI]*mesh_.time().deltaT().value()*rho[cellI]/3);
}

Foam::scalar
Foam::populationBalanceSubModels::oxidationModels::reduction
::characteristic(const scalar& abscissa, const label& cellI) 
{
    const volScalarField& concentration_O2(mesh_.lookupObject<volScalarField>("O2")); //use of pointeurs to avoid to reconstruct each times ???
    
    //const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    return (abscissa - cRed(cellI)*concentration_O2[cellI]*rho[cellI]*mesh_.time().deltaT().value()/3);
    
}

// ************************************************************************* //
