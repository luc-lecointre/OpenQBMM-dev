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
    VCarbon_("CarbonVolume", pow3(dimLength), 1.1e-29),
    Xi_("concentration hydrogen sites", inv(sqr(dimLength)), 1.7e19)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::convectionModels::growthReduction::~growthReduction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::volScalarField Foam::populationBalanceSubModels::convectionModels::growthReduction::Arrhenius
(
    const scalar& A,
    const scalar& n,
    const scalar& E
)
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    volScalarField Ea
    (
        IOobject
        (
            "Ea",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("activation energy", T.dimensions(), E*1000)
    );
    
    volScalarField Arrhenius = A*pow(T,n)*exp(-Ea/(8.314*T)); //careful, the dimensions of E in kJ and A divided by 1e6 to have the result by m^3 instead of cm^3
    
    Arrhenius.dimensions().reset(Foam::dimensionSet(-1,3,-1,0,0,0,0));//dimensions ???
    
    return Arrhenius;
}

Foam::scalar Foam::populationBalanceSubModels::convectionModels::growthReduction::Arrhenius
(
    const scalar& A,
    const scalar& n,
    const scalar& E,
    const scalar& T
)
{
    return A*pow(T,n)*exp(-E*1000/(8.314*T)); 
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::convectionModels::growthReduction
::cRed()
{
    volScalarField kReaction = Arrhenius(2.20e6,0.0,31.38); //reaction : Soot* + O2 -> Soot-H + 2CO + arrhenius law
    
    volScalarField kOH = Arrhenius(1.0e4,0.0,8.4); //reaction Soot-H + OH -> Soot-H + C0 + reaction probability 0.13 ???
    
    const volScalarField& concentration_O2(mesh_.lookupObject<volScalarField>("O2"));
    
    const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    Foam::tmp<Foam::volScalarField> cRed = (2.0*kReaction*concentration_O2/0.032+kOH*concentration_OH/0.017)*rho*Xi_*VCarbon_*pow(36.0*Foam::constant::mathematical::pi,1.0/3.0);

    
    
    //Info << "cRed : " << cRed.ref() << endl;
    //Info << "kReaction : " << kReaction << endl;
    return cRed;
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::convectionModels::growthReduction
::cOH()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    Foam::tmp<Foam::volScalarField> cOH = 1.29e3*sqrt(T)*concentration_OH*rho/0.017*0.13;
    
    cOH.ref().dimensions().reset(dimensionSet(0,-2,-1,0,0,0,0));
    
    return cOH;
}

Foam::tmp<Foam::volScalarField> Foam::populationBalanceSubModels::convectionModels::growthReduction
::cGro()
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    
    const volScalarField& rho = mesh_.lookupObject<volScalarField>("rho");
    
    const volScalarField& concentration_H(mesh_.lookupObject<volScalarField>("H"));
    
    const volScalarField& concentration_OH(mesh_.lookupObject<volScalarField>("OH"));
    
    const volScalarField& concentration_H2O(mesh_.lookupObject<volScalarField>("H2O"));
    
    const volScalarField& concentration_H2(mesh_.lookupObject<volScalarField>("H2"));
    
    const volScalarField& concentration_C2H2(mesh_.lookupObject<volScalarField>("C2H2"));
    
    volScalarField r
    (
        IOobject
        (
            "r",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", dimless , 0.0)
    );
    
    forAll(concentration_H, cellI)
    {
        if (concentration_H[cellI]!=0) 
        {
            r[cellI] = Arrhenius(1e2,1.80,68.42,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001+Arrhenius(8.68e-5,2.36,25.46,T[cellI])*concentration_OH[cellI]*rho[cellI]/0.005+Arrhenius(1.13e16,-0.06,476.05,T[cellI])/(Arrhenius(8.68e-2,2.36,25.46,T[cellI])*concentration_H2[cellI]*rho[cellI]/0.002+Arrhenius(6.44e-7,3.79,27.96,T[cellI])*concentration_H2O[cellI]*rho[cellI]/0.018+Arrhenius(4.17e6,0.15,0.00,T[cellI])*concentration_H[cellI]*rho[cellI]/0.001+Arrhenius(2.52e3,1.10,17.13,T[cellI])*concentration_C2H2[cellI]*rho[cellI]/0.026);
        }
    }
        
    Foam::tmp<Foam::volScalarField> cGro = 2.0*VCarbon_*r/(r+1.0)*concentration_C2H2*Xi_*Arrhenius(2.52e3,1.10,17.13)*rho/0.026*pow(36.0*Foam::constant::mathematical::pi,1.0/3.0);
    
    //Info << "cGro : " << cGro.ref() << endl;
    
    return cGro;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic() 
{
    volScalarField null
    (
        IOobject
        (
            "0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("zero", pow3(dimLength) , 0.0)
    );
    
    //Info << "characteristic in 0" << max(pow3((cRed()-cGro())/3.0*mesh_.time().deltaT()),null) << endl;
    //Info << "difference " << cRed() - cGro() << endl;
    
    return max(pow3((cRed()-cGro())/3.0*mesh_.time().deltaT()),null)*1.0e27;
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::convectionModels::growthReduction
::characteristic(const volScalarField& abscissa) 
{
    //return pow3(pow(abscissa,1.0/3.0) + cGro(cellI)/3.0*mesh_.time().deltaT().value());
    //Info << "abscissa" << abscissa.dimensions() << endl;
    
    return pow3(pow(abscissa*1.0e-27,1.0/3.0) + (cGro()-cRed())/3.0*mesh_.time().deltaT())*1.0e27;
}

// ************************************************************************* //
