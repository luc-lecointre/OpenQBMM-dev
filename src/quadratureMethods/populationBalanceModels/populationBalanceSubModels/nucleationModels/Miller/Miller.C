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
::dimerConcentration(const scalar& nC)
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(PAH_));
    
    volScalarField dimerSource = 0.5*Kfm(nC)*Foam::constant::physicoChemical::NA
    *sqr(pahConcentration*flThermo.rho()/0.202); //[mol/(m^3*s)]
    
    volScalarField betaN = Kfm(2*nC)*Foam::constant::physicoChemical::NA; //[m^3/(s*mol)]
    
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
    
    volScalarField V
    (
        IOobject
        (
            "volume",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    
    forAll(V, cellI){
        V[cellI]=mesh_.V()[cellI];
    }
    
    dimensionedScalar abscissaNucleation = volume(4*nCarbonPAH_)*1.0e27; //[nm^3]
        
    //Info << "abscissaNucleation" << tokens[specie] << " : " << abscissaNucleation.value() << endl;
        
    dimensionedScalar MPAH (
        "PAHMolarMass",
        Foam::dimensionSet(1,0,0,0,-1,0,0),
        0.202
    );
    
    const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(PAH_));
    
    tmp<volScalarField> jT = 0.5*Kfm(2.0*nCarbonPAH_)*sqr(Foam::constant::physicoChemical::NA*pahConcentration*flThermo.rho()/MPAH)*V*pow(abscissaNucleation,moment.order());
        
    //Info << V << endl;
    //Miller.ref().dimensions().reset(inv(dimTime*pow3(dimLength)));
    
    return jT;
}

// ************************************************************************* //
