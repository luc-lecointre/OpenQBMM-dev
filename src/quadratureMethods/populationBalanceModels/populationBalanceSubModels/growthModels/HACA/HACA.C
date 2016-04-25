/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Matteo Icardi
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

#include "HACA.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{
    defineTypeNameAndDebug(HACA, 0);

    addToRunTimeSelectionTable
    (
        growthModel,
        HACA,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::HACA
::HACA
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    growthModel(dict, mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModels::HACA
::~HACA()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::HACA::Arrhenius
(
    const scalar& A,
    const scalar& n,
    const scalar& Ea,
    const volScalarField& T
) const
{
    dimensionedScalar E (
        "ActivationEnergy",
        Foam::dimensionSet(0,0,0,1,0,0,0),
        Ea
    );
    
    tmp<volScalarField> Arrhenius =
        A*pow(T,n)*exp(-E/(8.314*T));
    
    Arrhenius.ref().dimensions().reset(dimensionSet(0,3,-1,0,-1,0,0));
    
    return Arrhenius;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::HACA::r
(const volScalarField& T) const
{
    
    const volScalarField& rho(mesh_.lookupObject<volScalarField>("rho"));
    
    const volScalarField& concentration_H(mesh_.lookupObject<volScalarField>("H"));
    
    dimensionedScalar MH (
        "HMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.001
    );
    
    const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    dimensionedScalar MOH (
        "OHMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.017
    );
    
    const volScalarField& concentration_H2O(mesh_.lookupObject<volScalarField>("H2O"));
    
    dimensionedScalar MH2O (
        "H2OMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.018
    );
    
    const volScalarField& concentration_H2(mesh_.lookupObject<volScalarField>("H2"));
    
    dimensionedScalar MH2 (
        "H2MolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.002
    );
    
    const volScalarField& concentration_C2H2(mesh_.lookupObject<volScalarField>("C2H2"));
    
    dimensionedScalar MC2H2 (
        "C2H2MolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.026
    );
    
    Foam::tmp<Foam::volScalarField> r = Arrhenius(1e8,1.80,68.42,T)*concentration_H*rho/MH+Arrhenius(8.68e4,2.36,25.46,T)*concentration_OH*rho/MOH+Arrhenius(1.13e16,-0.06,476.05,T)/(Arrhenius(8.68e4,2.36,25.46,T)*concentration_H2*rho/MH2+Arrhenius(6.44e-1,3.79,27.96,T)*concentration_H2O*rho/MH2O+Arrhenius(4.17e13,0.15,0.00,T)*concentration_H*rho/MH+Arrhenius(2.52e9,1.10,17.13,T)*concentration_C2H2*rho/MC2H2);

    r.ref().dimensions().reset(dimless);
    
    return r;
}


Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::growthModels::HACA::Kg
(
    const volScalarField& abscissa
) const
{   
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    const volScalarField& concentration_C2H2(mesh_.lookupObject<volScalarField>("C2H2"));
    
    dimensionedScalar MC2H2 (
        "C2H2MolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.026
    );
    
    dimensionedScalar Xi
    (
        "Xi", 
        dimensionSet(0,-2,0,0,0,0,0), 
        1.7e19
    );
    dimensionedScalar dCarbon
    (
        "dCarbon",
        dimensionSet(0,1,0,0,0,0,0),
        2.365e-10
    );
    
    return Arrhenius(2.52e9,1.10,17.13,flThermo.T())*concentration_C2H2*flThermo.rho()/MC2H2*2*dCarbon*abscissa*Xi*Foam::constant::physicoChemical::NA;
}

// ************************************************************************* //
