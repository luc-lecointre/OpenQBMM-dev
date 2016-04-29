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
    nCarbonPAH_(dict.lookup("nCarbonPAH")),
    rhoSoot_(
        "rhoSoot",
        Foam::dimensionSet(1,-3,0,0,0,0,0),
        1800.0e-27
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
    return nC*MCarbon_/(rhoSoot_*Foam::constant::physicoChemical::NA); // [nm^3]
}

Foam::volScalarField
Foam::populationBalanceSubModels::nucleationModels::Miller
::Kfm(const scalar& nC) 
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    return 8.8*sqrt(Foam::constant::mathematical::pi
        *Foam::constant::physicoChemical::k*1e18
        *flThermo.T()*Foam::constant::physicoChemical::NA
        /(nC*MCarbon_))
        *pow(6.0*nC*MCarbon_
        /(Foam::constant::mathematical::pi*rhoSoot_
        *Foam::constant::physicoChemical::NA), 2.0/3.0); // [nm^3/s]
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::dimerConcentration(const scalar& nC)
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(PAH_));
    
    volScalarField dimerSource = 0.5*Kfm(nC)*Foam::constant::physicoChemical::NA
    *sqr(pahConcentration*1e-27*flThermo.rho()/202.0); //[mol/(nm^3*s)]
    
    volScalarField betaN = Kfm(2*nC)*Foam::constant::physicoChemical::NA; //[nm^3/(s*mol)]
    
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
    
    tmp<volScalarField> Miller
    (
        new volScalarField
        (
            IOobject
            (
                "Miller",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    List<word> tokens;
    List<label> nCarbon;
    
    string::size_type lastPos = PAH_.find_first_not_of(" ", 0);
    string::size_type pos     = PAH_.find_first_of(" ", lastPos);
    


    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.append(PAH_.substr(lastPos, pos - lastPos));
        lastPos = PAH_.find_first_not_of(" ", pos);
        pos = PAH_.find_first_of(" ", lastPos);
    }
    
    lastPos = nCarbonPAH_.find_first_not_of(" ", 0);
    pos     = nCarbonPAH_.find_first_of(" ", lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        nCarbon.append(stoi(nCarbonPAH_.substr(lastPos, pos - lastPos)));
        lastPos = nCarbonPAH_.find_first_not_of(" ", pos);
        pos = nCarbonPAH_.find_first_of(" ", lastPos);
    }
    
    forAll (tokens, specie)
    {
        dimensionedScalar abscissaNucleation = volume(4*nCarbon[specie]); //[nm^3]
        
        //Info << "abscissaNucleation" << tokens[specie] << " : " << abscissaNucleation.value() << endl;
        
        dimensionedScalar MPAH (
            "PAHMolarMass",
            Foam::dimensionSet(1,0,0,0,-1,0,0),
            0.012*nCarbon[specie]+0.010
        );
    
        const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(tokens[specie]));
    
        volScalarField jT = 0.5*Kfm(2*nCarbon[specie])*sqr(Foam::constant::physicoChemical::NA*pahConcentration*flThermo.rho()/MPAH)*pow(abscissaNucleation,moment.order());
        
        Miller.ref().dimensions().reset
        (
            jT.dimensions()
        );

        Miller.ref() = Miller.ref() + jT;
    }
    
    //Miller.ref().dimensions().reset(inv(dimTime*pow3(dimLength)));
    
    return Miller;
}


/*Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::nucleationSource(const volUnivariateMoment& moment) 
{
    
    dimensionedScalar abscissaNucleation = pow(6.0/Foam::constant::mathematical::pi*2.0*MCarbon_*2.0*nCarbonPAH_/(rhoSoot_*Foam::constant::physicoChemical::NA),1.0/3.0);
        
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
        /(2.0*nCarbonPAH_*MCarbon_))
        *pow(12.0*nCarbonPAH_*MCarbon_
        /(Foam::constant::mathematical::pi*rhoSoot_
        *Foam::constant::physicoChemical::NA), 2.0/3.0)
        *sqr(pahConcentration*flThermo.rho()/MPAH*Foam::constant::physicoChemical::NA)*pow(abscissaNucleation, moment.order()); 
    
    //nucleationSource.ref().dimensions().reset(moment.dimensions()/dimTime);
    
    //Info << "deltaT : " << mesh_.time().deltaT()<< endl;
    
    //Info << "nucleationSource" << nucleationSource.ref()[0] << endl;
    
    //Info << "moment" << moment.order() << ":" << moment[0] << endl;
    
    return nucleationSource;
}*/

// ************************************************************************* //
