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
    nCarbonPAH_(readScalar(dict.lookup("nCarbonPAH"))),
    rhoSoot_(dict.lookup("rhoSoot"))
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
    
    return 8.8*sqrt(Foam::constant::mathematical::pi
        *Foam::constant::physicoChemical::k
        *flThermo.T()*Foam::constant::physicoChemical::NA
        /(nCarbonPAH_*MCarbon_))
        *pow(6.0*nCarbonPAH_*MCarbon_
        /(Foam::constant::mathematical::pi*rhoSoot_
        *Foam::constant::physicoChemical::NA), 2.0/3.0); // [m^3/s]
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::dimerConcentration(const scalar& nC)
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>("A4"));
    
    volScalarField dimerSource = 0.5*Kfm(nCarbonPAH_)*Foam::constant::physicoChemical::NA
    *sqr(pahConcentration*flThermo.rho()/202.0); //[mol/(m^3*s)]
    
    volScalarField betaN = Kfm(2*nCarbonPAH_)*Foam::constant::physicoChemical::NA; //[m^3/(s*mol)]
    
    dimensionedScalar MDimer (
        "dimerMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.404
    ); //or use A4 Molar Mass 
    
    return sqrt(dimerSource/betaN)*MDimer/flThermo.rho(); 
}
    
Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::nucleationSource(const volUnivariateMoment& moment) 
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    dimensionedScalar abscissaNucleation = pow(6.0/Foam::constant::mathematical::pi*volume(4*nCarbonPAH_),1.0/3.0); //[m]
    
    dimensionedScalar MPAH (
        "PAHMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.202
    );
    
    const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>("A4"));
    
    volScalarField jT = 0.5*Kfm(2*nCarbonPAH_)*sqr(Foam::constant::physicoChemical::NA*pahConcentration*flThermo.rho()/MPAH);
    
    return jT*pow(abscissaNucleation,moment.order());
}


/*Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::nucleationSource(const volUnivariateMoment& moment) 
{
    
    dimensionedScalar abscissaNucleation = pow(6.0/Foam::constant::mathematical::pi*2.0*MCarbon_*2*nCarbonPAH_/(rhoSoot_*Foam::constant::physicoChemical::NA),1.0/3.0);
        
    const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>("A4"));
    
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    dimensionedScalar MPAH (
        "PAHMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.202
    );

    tmp<volScalarField> nucleationSource = 4.4*sqrt(Foam::constant::mathematical::pi
        *Foam::constant::physicoChemical::k
        *flThermo.T()*Foam::constant::physicoChemical::NA
        /(nCarbonPAH_*MCarbon_))
        *pow(6.0*nCarbonPAH_*MCarbon_
        /(Foam::constant::mathematical::pi*rhoSoot_
        *Foam::constant::physicoChemical::NA), 2.0/3.0)
        *sqr(pahConcentration*flThermo.rho()/MPAH)*sqr(Foam::constant::physicoChemical::NA); 
        
    nucleationSource=nucleationSource*pow(abscissaNucleation, moment.order());
    
    nucleationSource.ref().dimensions().reset(moment.dimensions()/dimTime);
    
    //Info << "deltaT : " << mesh_.time().deltaT()<< endl;
    
    //Info << "nucleationSource" << nucleationSource.ref()[0] << endl;
    
    //Info << "moment" << moment.order() << ":" << moment[0] << endl;
    
    return nucleationSource;
}*/

// ************************************************************************* //
