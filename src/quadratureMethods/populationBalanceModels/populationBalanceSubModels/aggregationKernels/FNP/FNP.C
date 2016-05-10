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

#include "FNP.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(FNP, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        FNP,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::FNP
::FNP
(
    const dictionary& dict
)
:
    aggregationKernel(dict),
    m0Scaling_(readScalar(dict.lookup("m0Scaling"))),
    nUnimers_(readScalar(dict.lookup("nUnimers"))),
    cellVol_(dict.lookup("cellVol")),
    siteVol_(dict.lookup("siteVol")),
    alpha_(readScalar(dict.lookup("alpha"))),
    na_(readScalar(dict.lookup("na"))),
    nb_(readScalar(dict.lookup("nb"))),
    nuA_(readScalar(dict.lookup("nuA"))),
    nuB_(readScalar(dict.lookup("nuB"))),
    etaS_(dict.lookup("etaS")),
    Xip_(dict.lookup("Xip")),
    Xi0_(dict.lookup("Xi0")),
    T_(readScalar(dict.lookup("T")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::FNP
::~FNP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::FNP::beta
(
    const scalar& p,
    const scalar& i
) const
{
    // Free Coupling
    if (p == 1 && i == 1)
    {
        scalar R = pow(na_*siteVol_.value(), nuA_) + pow(nb_,nuB_)
            *pow(siteVol_.value(), 1.0/3.0);
        
        scalar D = Foam::constant::physicoChemical::k.value()*T_
            /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R);
        
        scalar Rcoll = pow(na_*siteVol_.value(), 1.0/3.0);
        
        return
            16.0*Foam::constant::mathematical::pi*D*Rcoll;
    }
    
    else
    {  
        scalar Rcoll_p = pow((p*na_*siteVol_.value()), 1.0/3.0);
        scalar Rcoll_i = pow((i*na_*siteVol_.value()), 1.0/3.0);
        
        scalar Rcor_p = pow(nb_,nuB_)*pow(p,(1.0 - nuB_)/2.0)
           *pow(siteVol_.value(), 1.0/3.0);
        scalar Rcor_i = pow(nb_,nuB_)*pow(i,(1.0 - nuB_)/2.0)
           *pow(siteVol_.value(), 1.0/3.0);
        
        scalar Ccor_p = 
            (3.0*p*nb_)/(4.0*Foam::constant::mathematical::pi
           *(pow3(Rcoll_p + Rcor_p) - pow3(Rcoll_p)));

        scalar Ccor_i =
            (3.0*i*nb_)/(4.0*Foam::constant::mathematical::pi
           *(pow3(Rcoll_i + Rcor_i) - pow3(Rcoll_i)));
        
        // Unimer Insertion
        if (p > 1 && i == 1)
        {
            scalar R_p = Rcoll_p + Rcor_p;
            scalar R_i = pow(na_*siteVol_.value(), nuA_) + pow(nb_,nuB_)
                *pow(siteVol_.value(), 1.0/3.0);
        
            scalar D_p = Foam::constant::physicoChemical::k.value()*T_
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_p);
            scalar D_i = Foam::constant::physicoChemical::k.value()*T_
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_i);
            
            scalar Cp = 1.0/(siteVol_.value()*pow(nb_, (3.0*nuB_ - 1.0)));
            scalar Dstar = D_i*pow(Cp/Ccor_p, 1.5);
            scalar Ains = exp(-alpha_*sqrt(p));
            
            return
                 4.0*Foam::constant::mathematical::pi
                *Ains*((D_p + D_i)*Dstar*(Rcoll_p + Rcoll_i)
                *(Rcoll_p + Rcoll_i + Rcor_p)/(Dstar*(Rcoll_p + Rcoll_i) 
               + (D_p + D_i)*Rcor_p));
        }
        
        // Unimer Insertion
        else if (p == 1 && i > 1)
        {
            scalar R_p = pow(na_*siteVol_.value(), nuA_) + pow(nb_,nuB_)
                *pow(siteVol_.value(), 1.0/3.0);
            scalar R_i = Rcoll_i + Rcor_i;
        
            scalar D_p = Foam::constant::physicoChemical::k.value()*T_
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_p);
            scalar D_i = Foam::constant::physicoChemical::k.value()*T_
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_i);
            
            scalar Cp = 1.0/(siteVol_.value()*pow(nb_, (3.0*nuB_ - 1.0)));
            scalar Dstar = D_p*pow(Cp/Ccor_i, 1.5);
            scalar Ains = exp(-alpha_*sqrt(i));
            
            return
                 4.0*Foam::constant::mathematical::pi
                *Ains*((D_p + D_i)*Dstar*(Rcoll_p + Rcoll_i)
                *(Rcoll_p + Rcoll_i + Rcor_i)/(Dstar*(Rcoll_p + Rcoll_i) 
               + (D_p + D_i)*Rcor_i));
        }
        
        // Large Aggregate Fusion
        else if (p > 1 && i > 1)
        {
            scalar R_p = Rcoll_p + Rcor_p;
            scalar R_i = Rcoll_i + Rcor_i;
            
            scalar D_p = Foam::constant::physicoChemical::k.value()*T_
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_p);
            scalar D_i = Foam::constant::physicoChemical::k.value()*T_
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_i);
                
            scalar Chi = 
                pow(siteVol_.value(), 1.0/3.0)
               /pow((Ccor_p + Ccor_i)*siteVol_.value(), 3.0/4.0);
                
            scalar L = sqr(Rcor_p + Rcor_i)/Chi;
            scalar N = nb_/pow(Chi/pow(siteVol_.value(), 1.0/3.0), 5.0/3.0);
            
            scalar Dfus =
                Foam::constant::physicoChemical::k.value()*T_
               *sqr(Rcor_p + Rcor_i)/(etaS_.value()*N*Chi*sqr(L));
            
            scalar Afus = exp(-alpha_*min(p, i)*sqrt(max(p,i)));
            
            return
                 4.0*Foam::constant::mathematical::pi*Afus*((D_p + D_i)
                *Dfus*(Rcoll_p + Rcoll_i)*(Rcoll_p + Rcoll_i + Rcor_p + Rcor_i))
                /(Dfus*(Rcoll_p + Rcoll_i) + (D_p + D_i)*(Rcor_p + Rcor_i));
        }
        
        else
        {
            return scalar(0);
        }
    }
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::aggregationKernels::FNP::Ka
(
    const volScalarField& a1,
    const volScalarField& a2
) const
{
    dimensionedScalar smallA("smallA", a1.dimensions(), 0);
    volScalarField abscissa1 = max(a1, smallA);
    volScalarField abscissa2 = max(a2, smallA);
    
    if (min(abscissa1).value() < 0 || min(abscissa2).value() < 0)
    {
        FatalErrorInFunction
            << "Negative abscissa value."
            << abort(FatalError);
    }
    
    /*if (!abscissa1.mesh().foundObject<fluidThermo>(basicThermo::dictName))
    {
        FatalErrorInFunction
            << "No valid thermophysical model found."
            << abort(FatalError);
    }
    
    const fluidThermo& flThermo =
        abscissa1.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
        
    const volScalarField& T = flThermo.T();*/
    volScalarField mixtureFraction //=
        //abscissa1.mesh().lookupObject<volScalarField>("mixtureFraction");
    (
        IOobject
        (
            "mixtureFraction",
            "0",
            a1.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        a1.mesh()
    );
    tmp<volScalarField> betaKernel
    (
        new volScalarField
        (
            IOobject
            (
                "betaKernel",
                a1.mesh().time().timeName(),
                a1.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            a1.mesh(),
            dimensionedScalar("beta", dimless, 0.0)
        )
    );

    forAll(abscissa1, cellI)
    {
        if
        (
            abscissa1[cellI] == 0 
         || abscissa2[cellI] == 0
         || mixtureFraction[cellI] < Xip_.value()
         || mixtureFraction[cellI] > Xi0_.value()
        )
        {
            betaKernel.ref()[cellI] = scalar(0);
        }
        else
        {
            Info<< mixtureFraction[cellI] << endl;
            scalar p = max(1.0, abscissa1[cellI]);
            scalar i = max(1.0, abscissa2[cellI]);
            
            scalar pLow  = floor(p);
            scalar pHigh = ceil(p);
            scalar iLow  = floor(i);
            scalar iHigh = ceil(i);
            
            if(pLow == pHigh && iLow == iHigh)
            {
                betaKernel.ref()[cellI] =
                    beta(p,i);
            }
            else if (pLow == pHigh)
            {
                scalar BpLow = beta(p, iLow);
                scalar BpHigh = beta(p, iHigh);
                
                betaKernel.ref()[cellI] = BpLow + (iHigh - i)
                   *(BpHigh - BpLow)/(iHigh - iLow);
            }
            else if (iLow == iHigh)
            {
                scalar BiLow = beta(pLow, i);
                scalar BiHigh = beta(pHigh, i);
                
                betaKernel.ref()[cellI] = BiLow + (pHigh - p)
                   *(BiHigh - BiLow)/(pHigh - pLow);
            }
            else
            {
                scalar B11 = beta(pLow, iLow);
                scalar B12 = beta(pLow, iHigh);
                scalar B21 = beta(pHigh, iLow);
                scalar B22 = beta(pHigh, iHigh);
                  
                betaKernel.ref()[cellI] = 
                    (B11*(pHigh - p)*(iHigh - i) + B21*(p - pLow)*(iHigh - i)
                  + B12*(pHigh - p)*(i - iLow) + B22*(p - pLow)*(i - iLow))
                   /((pHigh - pLow)*(iHigh - iLow));
            }
        }
    }
    
    betaKernel.ref().dimensions().reset(pow3(abscissa1.dimensions())/dimTime);
    return betaKernel*m0Scaling_;
    
}

// ************************************************************************* //
