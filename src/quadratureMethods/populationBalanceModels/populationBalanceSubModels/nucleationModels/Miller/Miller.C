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
    gamma_(dict.lookup("gamma")),
    rhoSoot_(
        "rhoSoot",
        Foam::dimensionSet(1,-3,0,0,0,0,0),
        1800.0
    ),
    V_(
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
    )
{
    
    forAll(V_, cellI){
        V_[cellI]= mesh_.V()[cellI];
    }
    
    string::size_type lastPos = PAH_.find_first_not_of(" ", 0);
    string::size_type pos     = PAH_.find_first_of(" ", lastPos);
    
    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.append(PAH_.substr(lastPos, pos - lastPos));
        lastPos = PAH_.find_first_not_of(" ", pos);
        pos = PAH_.find_first_of(" ", lastPos);
    }
    
    //Info << "PAH : "<<tokens << endl;
    
    lastPos = nCarbonPAH_.find_first_not_of(" ", 0);
    pos     = nCarbonPAH_.find_first_of(" ", lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        nCarbon.append(stoi(nCarbonPAH_.substr(lastPos, pos - lastPos)));
        lastPos = nCarbonPAH_.find_first_not_of(" ", pos);
        pos = nCarbonPAH_.find_first_of(" ", lastPos);
    }
    
    //Info << "nCarbonPAH : "<< nCarbon << endl;
    
    lastPos = gamma_.find_first_not_of(" ", 0);
    pos     = gamma_.find_first_of(" ", lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        gamma.append(stod(gamma_.substr(lastPos, pos - lastPos)));
        lastPos = gamma_.find_first_not_of(" ", pos);
        pos = gamma_.find_first_of(" ", lastPos);
    }
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
::volume(const volScalarField& nCmoy)
{
    return nCmoy*MCarbon_/(rhoSoot_*Foam::constant::physicoChemical::NA); // [m^3]
}

Foam::volScalarField
Foam::populationBalanceSubModels::nucleationModels::Miller
::Kfm(const scalar& nC) 
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    return pow(6.0,2.0/3.0)*8.8*sqrt(Foam::constant::physicoChemical::k
        *flThermo.T()/flThermo.rho())*pow(volume(nC)/Foam::constant::mathematical::pi, 1.0/6.0); // [m^3/s]
}

Foam::volScalarField
Foam::populationBalanceSubModels::nucleationModels::Miller
::Kfm(const volScalarField& nCmoy) 
{
    const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
    
    return pow(6.0,2.0/3.0)*8.8*sqrt(Foam::constant::physicoChemical::k
        *flThermo.T()/flThermo.rho())*pow(volume(nCmoy)/Foam::constant::mathematical::pi, 1.0/6.0); // [m^3/s]
}

//with condensation

Foam::volScalarField
Foam::populationBalanceSubModels::nucleationModels::Miller
::nCmoy()
{
    volScalarField Cmoy
    (
        IOobject
        (
            "Cmoy",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );
        
    forAll (tokens, sI)
    {
    
        const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(tokens[sI]));
        
        //specie PAH = mesh_.lookupObject<specie>(tokens[sI]);
        
        scalar MPAH = nCarbon[sI]*0.012+0.01; 
    
        const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
        
        //Info << betaPAH() << endl;
        
        forAll(Cmoy,mI)
        {
            Cmoy[mI] += 0.5*gamma[sI]*Kfm(nCarbon[sI])[mI]*Foam::constant::physicoChemical::NA.value()*sqr(pahConcentration[mI]*rho[mI]/MPAH)*nCarbon[sI]/betaPAH().ref()[mI];
        }
    }
    
    return Cmoy;
    
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
    
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    volScalarField Cmoy=nCmoy();
    
    forAll(abscissa, cellI)
    {
        if (abscissa[cellI]!=0)
        {
            scalar dDimer = pow(6.0/Foam::constant::mathematical::pi*volume(2*Cmoy[cellI]).value(),1.0/3.0);
            
            scalar dSoot = pow(6.0/Foam::constant::mathematical::pi*abscissa[cellI]*1e-27,1.0/3.0);
    
            scalar mu = rhoSoot_.value()*(volume(2*Cmoy[cellI]).value()*abscissa[cellI]*1e-27)/(volume(2*Cmoy[cellI]).value()+abscissa[cellI]*1e-27);
    
            beta.ref()[cellI] = sqrt((Foam::constant::mathematical::pi*Foam::constant::physicoChemical::k.value()*T[cellI])/(2.0*mu))*sqr(dDimer+dSoot);
        }
    }
    
    return beta; //(m^3/s)
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::betaN()
{
    volScalarField cMoy = nCmoy();
    
    return Kfm(2*cMoy)*Foam::constant::physicoChemical::NA;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::betaPAH()
{

    tmp<volScalarField> betaPAH
    (
        new volScalarField
        (
            IOobject
            (
                "betaPAH",
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
    
    forAll (tokens, sI)
    {
    
        const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(tokens[sI]));
        
        //specie PAH = mesh_.lookupObject<specie>(tokens[sI]);
        
        dimensionedScalar MPAH (
            "PAHMolarMass",
            Foam::dimensionSet(1,0,0,0,-1,0,0),
            nCarbon[sI]*0.012+0.01                
            //PAH.W()*1e-3
        );
    
        const fluidThermo& flThermo = mesh_.lookupObject<fluidThermo>(basicThermo::dictName);
        
        tmp<volScalarField> beta = 1.0/2.0*gamma[sI]*Kfm(nCarbon[sI])*Foam::constant::physicoChemical::NA*sqr(pahConcentration*flThermo.rho()/MPAH);
        
        betaPAH.ref().dimensions().reset(beta().dimensions());
        betaPAH = betaPAH + beta;
        
    }
    
    return betaPAH; 
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::xiNuc()
{

    volScalarField cMoy = nCmoy();
    
    return volume(4*cMoy)*1.0e27;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::nucleationModels::Miller
::xiCond
(
    const volScalarField& abscissa,
    const label& order
)
{
    volScalarField cMoy = nCmoy();
    
    return pow(abscissa+volume(4*cMoy)*1.0e27,order)-pow(abscissa,order);
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
    
//without condensation
    
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
    
    forAll (tokens, sI)
    {
        dimensionedScalar abscissaNucleation = volume(4*nCarbon[sI])*1.0e27; //[m]
        
        //Info << "abscissaNucleation" << tokens[specie] << " : " << abscissaNucleation.value() << endl;
        
        //specie PAH = mesh_.lookupObject<specie>(tokens[sI]);
        
        dimensionedScalar MPAH (
            "PAHMolarMass",
            Foam::dimensionSet(1,0,0,0,-1,0,0),
            nCarbon[sI]*0.012+0.01
            //PAH.W()*1e-3
        );
        
    
        const volScalarField& pahConcentration(mesh_.lookupObject<volScalarField>(tokens[sI]));
    
        volScalarField jT = gamma[sI]*0.5*Kfm(2*nCarbon[sI])*sqr(Foam::constant::physicoChemical::NA*pahConcentration*flThermo.rho()/MPAH)*V_*pow(abscissaNucleation,moment.order());
        
        Miller.ref().dimensions().reset
        (
            jT.dimensions()
        );

        Miller.ref() = Miller.ref() + jT;
        
        //Info << jT << endl;
        //Info << tokens[sI] << " : " << nCarbon[sI] << " gamma : " << gamma[sI] << endl;
    }
    
    //Miller.ref().dimensions().reset(inv(dimTime*pow3(dimLength)));
    
    return Miller;
}

// ************************************************************************* //
