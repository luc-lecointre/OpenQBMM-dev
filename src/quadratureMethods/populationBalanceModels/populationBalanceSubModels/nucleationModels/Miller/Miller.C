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

#include "Miller.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(Miller, 0);

    addToRunTimeSelectionTable
    (
        nucleationModel,
        Miller,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::Miller::Miller
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    nucleationModel(dict, mesh),
    //MCarbon_(dict.lookup("MCarbon")),
    MCarbon_
    (
        "carbonMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.012
    ),
    PAH_(dict.lookup("PAH")),
    nCarbonPAH_(readScalar(dict.lookup("nCarbonPAH"))),
    rhoSoot_(
        "rhoSoot",
        Foam::dimensionSet(1,-3,0,0,0,0,0),
        1800.0
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::nucleationModels::Miller::~Miller()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::dimensionedScalar
Foam::populationBalanceSubModels::nucleationModels::Miller
::volume(const scalar& nC)
{
    return nC*MCarbon_/(rhoSoot_*Foam::constant::physicoChemical::NA); // [m^3]
}

Foam::volScalarField
Foam::populationBalanceSubModels::nucleationModels::Miller
::Kfm(const scalar& nC) 
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    return pow(6.0,2.0/3.0)*8.8*sqrt(Foam::constant::physicoChemical::k
        *flThermo.T()/flThermo.rho())*pow(volume(nC)/Foam::constant::mathematical::pi, 1.0/6.0); // [m^3/s]
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::beta
(
    const volScalarField& abscissa
)
{
    tmp<volScalarField> beta
    (
        new volScalarField
        (
            IOobject
            (
                "beta",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("beta", dimensionSet(0,3,-1,0,0,0,0), 0.0)
        )
    );
    
    dimensionedScalar dDimer = pow(6.0/Foam::constant::mathematical::pi*volume(2*nCarbonPAH_),1.0/3.0);
    
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    forAll(abscissa, cellI)
    {
        if (abscissa[cellI]!=0)
        {
            scalar dSoot = pow(6.0/Foam::constant::mathematical::pi*abscissa[cellI],1.0/3.0);
    
            scalar mu = rhoSoot_.value()*(volume(2*nCarbonPAH_).value()*abscissa[cellI])/(volume(2*nCarbonPAH_).value()+abscissa[cellI]);
    
            beta.ref()[cellI] = sqrt((Foam::constant::mathematical::pi*Foam::constant::physicoChemical::k.value()*T[cellI])/(2.0*mu))*sqr(dDimer.value()+dSoot);
        }
    }
    
    return beta; //(m^3/s)
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::betaN()
{
    return Kfm(2*nCarbonPAH_)*Foam::constant::physicoChemical::NA;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::betaPAH()
{
    const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(PAH_));
    
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    dimensionedScalar MPAH (
        "PAHMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.202
    );
    
    return 1.0/2.0*Kfm(nCarbonPAH_)*Foam::constant::physicoChemical::NA*sqr(pahConcentration*flThermo.rho()/MPAH);
}

Foam::dimensionedScalar
Foam::populationBalanceSubModels::nucleationModels::Miller
::xiNuc()
{
    return volume(4*nCarbonPAH_)*1.0e27;
}


// ************************************************************************* //
