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
    oxidationModel(dict, mesh),
    cRed_(readScalar(dict.lookup("cRed")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::oxidationModels::reduction::~reduction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar
Foam::populationBalanceSubModels::oxidationModels::reduction
::characteristic(const label& cellI) 
{
    const volScalarField& concentration(mesh_.lookupObject<volScalarField>("OH"));
    
    //const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    return (cRed_*concentration[cellI]*mesh_.time().deltaT().value()/3);
}

Foam::scalar
Foam::populationBalanceSubModels::oxidationModels::reduction
::characteristic(const scalar& abscissa, const label& cellI) 
{
    const volScalarField& concentration(mesh_.lookupObject<volScalarField>("OH"));
    
    //const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    return (abscissa - cRed_*concentration[cellI]*mesh_.time().deltaT().value()/3);
    
}

// ************************************************************************* //
