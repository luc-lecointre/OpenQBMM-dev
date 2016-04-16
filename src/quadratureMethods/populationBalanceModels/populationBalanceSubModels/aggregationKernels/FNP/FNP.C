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
    Xi0_(dict.lookup("Xi0"))
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
    const scalar& i,
    const scalar& T
) const
{
    dimensionedScalar Xi("Xi", dimless, 1.0);

    // Free Coupling
    if (p == 1 && i == 1)
    {
        scalar R = pow(na_*siteVol_.value(), nuA_) + pow(nb_,nuB_)
            *pow(siteVol_.value(), 1.0/3.0);
        
        scalar D = Foam::constant::physicoChemical::k.value()*T
            /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R);
        
        scalar Rcoll = pow(p*na_*siteVol_.value(), 1.0/3.0);
        
        return
            16.0*Foam::constant::mathematical::pi*pos(Xi - Xip_).value()*D*Rcoll;
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
            scalar R_p = Rcoll_p + Rcor_p
            scalar R_i = pow(na_*siteVol_.value(), nuA_) + pow(nb_,nuB_)
                *pow(siteVol_.value(), 1.0/3.0);
        
            scalar D_p = Foam::constant::physicoChemical::k.value()*T
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_p);
            scalar D_i = Foam::constant::physicoChemical::k.value()*T
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_i);
            
            scalar Cstar = 1/(siteVol_.value()*pow(nb_, (3.0*nuB_ - 1.0)));
            scalar Cratio = Cstar/Ccor_p;
            scalar Dstar = D_i*pow(Cratio, 1.5);
            scalar Ains = Foam::exp(-alpha_*Foam::sqrt(p));
            
            return
                 4.0*Foam::constant::mathematical::pi*pos(Xi - Xip_).value()
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
        
            scalar D_p = Foam::constant::physicoChemical::k.value()*T
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_p);
            scalar D_i = Foam::constant::physicoChemical::k.value()*T
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_i);
            
            scalar Cstar = 1/(siteVol_.value()*pow(nb_, (3.0*nuB_ - 1.0)));
            scalar Cratio = Cstar/Ccor_i;
            scalar Dstar = D_i*pow(Cratio, 1.5);
            scalar Ains = Foam::exp(-alpha_*Foam::sqrt(i));
            
            return
                 4.0*Foam::constant::mathematical::pi*pos(Xi - Xip_).value()
                *Ains*((D_p + D_i)*Dstar*(Rcoll_p + Rcoll_i)
                *(Rcoll_p + Rcoll_i + Rcor_i)/(Dstar*(Rcoll_p + Rcoll_i) 
                + (D_p + D_i)*Rcor_i));
        }
        
        // Large Aggregate Fusion
        else if (p > 1 && i > 1)
        {
            scalar R_p = Rcoll_p + Rcor_p;
            scalar R_i = Rcoll_i + Rcor_i;
            
            scalar D_p = Foam::constant::physicoChemical::k.value()*T
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_p);
            scalar D_i = Foam::constant::physicoChemical::k.value()*T
                /(6.0*Foam::constant::mathematical::pi*etaS_.value()*R_i);
                
            scalar Chi = 
                pow(siteVol_.value(), 1.0/3.0)
               /pow((Ccor_p + Ccor_i)*siteVol_.value(), 3.0/4.0);
                
            scalar L = sqr(Rcor_p + Rcor_i)/Chi;
            scalar N = nb_/(Chi/pow(siteVol_, 1.0/3.0), 5.0/3.0);
            
            scalar Dfus =
                (Foam::constant::physicoChemical::k.value()*T
               *sqr(Rcor_p + Rcor_i)
               /(etaS_.value()*N*Chi*sqr(L));
            
            scalar Afus = Foam::exp(-alpha_*min(p, i)*sqrt(max(p,i)));
            
            return
                 4.0*Foam::constant::mathematical::pi
                *pos(Xi - Xip_).value()*Afus*((D_p + D_i)
                *Dfus*(Rcoll_p + Rcoll_i)*(Rcoll_p + Rcoll_i + Rcor_p + Rcor_i))
                /(Dfus*(Rcoll_p + Rcoll_i) + (D_p + D_i)*(Rcor_p + Rcor_i));
        }
    }
}

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::aggregationKernels::FNP::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{
    if (min(abscissa1).value() < 0 || min(abscissa2).value() < 0)
    {
        FatalErrorIn
        (
            "Foam::populationBalanceSubModels::aggregationKernels::FNP::"
            "Ka\n"
            "(\n"
            "   const volScalarField& p\n"
            "   const volScalarField& i\n"
            ")"
        )   << "zero or negative abscissa value"
            << abort(FatalError);
    }
    
    if (!abscissa1.mesh().foundObject<fluidThermo>(basicThermo::dictName))
    {
        FatalErrorIn
        (
            "Foam::populationBalanceSubModels::aggregationKernels::FNP::"
            "beta\n"
            "(\n"
            "   const volScalarField& p,\n "
            "   const volScalarField& i\n"
            ")"
        )   << "No valid thermophysical model found."
            << abort(FatalError);
    }
    
    const fluidThermo& flThermo =
        abscissa1.mesh().lookupObject<fluidThermo>(basicThermo::dictName);
        
    const volScalarField& T = flThermo.T();
    
    tmp<volScalarField> betaKernel
    (
        new volScalarField
        (
            IOobject
            (
                "betaKernel",
                abscissa1.mesh().time().timeName(),
                abscissa1.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            abscissa1.mesh(),
            dimensionedScalar("betaKernel", pow3(abscissa1.dimensions())/dimTime, 0.0)
        )
    );

    forAll(abscissa1, cellI)
    {
        if (abscissa1[cellI] == 0 || abscissa2[cellI] == 0)
        {
            betaKernel.ref()[cellI] = 0;
        }
        else
        {
            scalar p = max(1.0, abscissa1[cellI]);
            scalar i = max(1.0, abscissa2[cellI]);
            
            scalar pLow  = floor(p);
            scalar pHigh = ceil(p);
            scalar iLow  = floor(i);
            scalar iHigh = ceil(i);
            
            if(pLow == pHigh && iLow == iHigh)
            {
                betaKernel.ref()[cellI] = beta(p,i, T[cellI]);
            }
            else if (pLow == pHigh)
            {
                scalar BpLow = beta(p, iLow, T[cellI]);
                scalar BpHigh = beta(p, iHigh, T[cellI]);
                
                betaKernel.ref()[cellI] = BpLow + (iHigh - i)
                   *(BpHigh - BpLow)/(iHigh - iLow);
            }
            else if (iLow == iHigh)
            {
                scalar BiLow = beta(pLow, i, T[cellI]);
                scalar BiHigh = beta(pHigh, i, T[cellI]);
                
                betaKernel.ref()[cellI] = BiLow + (pHigh - p)
                   *(BiHigh - BiLow)/(pHigh - pLow);
            }
            else
            {
                scalar B11 = beta(pLow, iLow, T[cellI]);
                scalar B12 = beta(pLow, iHigh, T[cellI]);
                scalar B21 = beta(pHigh, iLow, T[cellI]);
                scalar B22 = beta(pHigh, iHigh, T[cellI]);
            
                betaKernel.ref()[cellI] = 
                    (B11*(pHigh - p)*(iHigh - i) + B21*(p - pLow)*(iHigh - i)
                  + B12*(pHigh - p)*(i - iLow) + B22*(p - pLow)*(i - iLow))
                   /((pHigh - pLow)*(iHigh - iLow));
            }
        }
    }
    return betaKernel*9.00e23;
    
}

// ************************************************************************* //
