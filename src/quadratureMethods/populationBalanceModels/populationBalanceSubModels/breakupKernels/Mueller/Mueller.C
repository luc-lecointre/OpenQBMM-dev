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

#include "Mueller.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(Mueller, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        Mueller,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::Mueller
::Mueller
(
    const dictionary& dict
)
:
    breakupKernel(dict),
    Xi_("concentration hydrogen sites", inv(sqr(dimLength)), 17.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::Mueller::~Mueller()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::breakupKernels::Mueller::Kb
(
    const volScalarField& abscissa
) const
{   
    Foam::tmp<Foam::volScalarField> breakup
    (
        new volScalarField
        (
            IOobject
            (
                "breakup",
                abscissa.mesh().time().timeName(),
                abscissa.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            abscissa.mesh(),
            dimensionedScalar("zero", inv(dimTime), 0.0)
        )
    );
    
    const volScalarField& T = abscissa.mesh().lookupObject<volScalarField>("T");
    
    const volScalarField& rho = abscissa.mesh().lookupObject<volScalarField>("rho");
    
    const volScalarField& concentration_O2(abscissa.mesh().lookupObject<volScalarField>("O2"));
    
    forAll (abscissa, cellI)
    {
        if (abscissa[cellI] != 0.0)
        {
            breakup.ref()[cellI] = pow(36*Foam::constant::mathematical::pi,1.0/3.0)*2.20e6*exp(-31380/(8.314*T[cellI]))*Xi_*concentration_O2[cellI]*rho[cellI]/(0.032*pow(abscissa[cellI],1.0/3.0));
        }
    }
    
    return breakup;
}

// ************************************************************************* //
