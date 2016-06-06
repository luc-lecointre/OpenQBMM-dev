/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "univariatePopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(univariatePopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        univariatePopulationBalance,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::univariatePopulationBalance
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    univariatePDFTransportModel(name, dict, U.mesh(), U, phi, "RPlus"),
    populationBalanceModel(name, dict, U, phi),
    name_(name),
    aggregation_(dict.lookup("aggregation")),
    breakup_(dict.lookup("breakup")),
    convection_(dict.lookup("convection")),
    nucleation_(dict.lookup("nucleation")),
    aggregationKernel_
    (
        Foam::populationBalanceSubModels::aggregationKernel::New
        (
            dict.subDict("aggregationKernel")
        )
    ),
    breakupKernel_
    (
        Foam::populationBalanceSubModels::breakupKernel::New
        (
            dict.subDict("breakupKernel")
        )
    ),
    daughterDistribution_
    (
        Foam::populationBalanceSubModels::daughterDistribution::New
        (
            dict.subDict("daughterDistribution")
        )
    ),
    convectionModel_
    (
        Foam::populationBalanceSubModels::convectionModel::New
        (
            dict.subDict("convectionModel"),
            U.mesh()
        )
    ),
    diffusionModel_
    (
        Foam::populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    nucleationModel_
    (
        Foam::populationBalanceSubModels::nucleationModel::New
        (
            dict.subDict("nucleationModel"),
            U.mesh()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::~univariatePopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::aggregationSource
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> aSource
    (
        new volScalarField
        (
            IOobject
            (
                "aSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    if (!aggregation_)
    {
        aSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return aSource;
    }
    label order = moment.order();

    volScalarField& aggregationSource = aSource.ref();

    forAll(quadrature_.nodes(), pNode1I)
    {
        const extendedVolScalarNode& node1 = quadrature_.nodes()[pNode1I];
    
        const volScalarField& pWeight1 = node1.primaryWeight();
        
        //Info << "weight first order : " << pWeight1.dimensions() << endl;

        forAll(node1.secondaryWeights(), sNode1I)
        {

            const volScalarField& sWeight1 = node1.secondaryWeights()[sNode1I];
            
            //Info << "weight second order : " << sWeight1.dimensions() << endl;
                        
            const volScalarField& sAbscissa1
                = node1.secondaryAbscissae()[sNode1I];
            
            //Info << "abscissa" << sAbscissa1 << endl;

            forAll(quadrature_.nodes(), pNode2I)
            {
                const extendedVolScalarNode& node2 = quadrature_.nodes()[pNode2I];

                const volScalarField& pWeight2 = node2.primaryWeight();

                forAll(node2.secondaryWeights(), sNode2I)
                {
                    const volScalarField& sWeight2
                        = node2.secondaryWeights()[sNode2I];

                    const volScalarField& sAbscissa2
                        = node2.secondaryAbscissae()[sNode2I];

                    tmp<volScalarField> aggInnerSum =
                        0.5*pWeight1*sWeight1*
                        (
                            pWeight2*sWeight2*
                            (
                                pow // Birth
                                (
                                    sAbscissa1 + sAbscissa2,
                                    order
                                )
                              - pow(sAbscissa1, order)
                              - pow(sAbscissa2, order)
                            )*aggregationKernel_->Ka(sAbscissa1, sAbscissa2)
                        );
                        
                    //Info << "aggregationKernel" << aggregationKernel_->Ka(sAbscissa1, sAbscissa2).ref().dimensions() << endl;
                    
                    //Info << "sWeight1" << sWeight1.dimensions() << endl;

                    aggregationSource.dimensions().reset
                    (
                        aggInnerSum().dimensions()
                    );

                    aggregationSource == aggregationSource + aggInnerSum();
                }
            }
        }
    }
    
    //Info << "aggregationSource : "<< aggregationSource << endl;
    
    return aSource;
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::breakupSource
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> bSource
    (
        new volScalarField
        (
            IOobject
            (
                "bSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    if (!breakup_)
    {
        bSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return bSource;
    }

    label order = moment.order();

    volScalarField& breakupSource = bSource.ref();

    forAll(nodes_(), pNodeI)
    {
        const extendedVolScalarNode& node = nodes_()[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            tmp<volScalarField> bSrc = node.primaryWeight()
                *node.secondaryWeights()[sNodei]
                *breakupKernel_->Kb(node.secondaryAbscissae()[sNodei])
                *(
                    daughterDistribution_->mD                      //Birth
                    (
                        order,
                        node.secondaryAbscissae()[sNodei]
                    )
                  - pow(node.secondaryAbscissae()[sNodei], order)   //Death
                 );

            breakupSource.dimensions().reset(bSrc().dimensions());
            breakupSource == breakupSource + bSrc;
        }
    }

    return bSource;
}

Foam::tmp<fvScalarMatrix> Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::momentDiffusion
(
    const volUnivariateMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}

void Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance::phaseSpaceConvection()
{
    if (!convection_)
    {
        return;
    }
    else
    {
        forAll(quadrature_.moments(), mI)
        {
            volUnivariateMoment& moment = quadrature_.moments()[mI];
        
            forAll(moment, cellI)
            {
                moment[cellI]=0.0;

                forAll(quadrature_.nodes(), pNodeI)
                {
                    const extendedVolScalarNode& node = quadrature_.nodes()[pNodeI];
                    volScalarField primaryAbscissa = node.primaryAbscissa();
                    //Info << "primaryAbscissa : " << primaryAbscissa << endl;
                    volScalarField sigma = node.sigma();
                    //Info << "sigma : " << sigma << endl;
        
                    scalar characteristic = convectionModel_->characteristic(cellI);
                    //Info << "characteristic in 0 : " << characteristic << endl;
                    scalar primaryAbscissaFinal = convectionModel_->characteristic(primaryAbscissa[cellI], cellI);
        
                    scalar xmax=characteristic+10.0;
            
                    if (sigma[cellI]!=0 && characteristic!=0.0)
                    { 
                        while (quadrature_.momentInverter()->distribution(xmax,primaryAbscissa[cellI],sigma[cellI])>1.0e-12)
                        {
                            xmax+=10.0;
                        }
                        
                        //Formule de quadrature de Gauss_3_points
                        
                        /*moment[cellI] += node.primaryWeight()[cellI]*(xmax -characteristic)/2.0*
                        (8.0/9.0*pow(convectionModel_->characteristic((xmax+characteristic)/2.0,cellI),mI)
                        *quadrature_.momentInverter()->distribution((xmax+characteristic)/2.0,primaryAbscissa[cellI],sigma[cellI])
                        + 5.0/9.0*(pow(convectionModel_->characteristic((xmax-characteristic)/2.0*sqrt(3.0/5.0)+(xmax+characteristic)/2.0,cellI),mI)*quadrature_.momentInverter()->distribution((xmax-characteristic)/2.0*sqrt(3.0/5.0)+(xmax+characteristic)/2.0,primaryAbscissa[cellI],sigma[cellI])
                        +pow(convectionModel_->characteristic((xmax-characteristic)/2.0*(-sqrt(3.0/5.0))+(xmax+characteristic)/2.0,cellI),mI)*quadrature_.momentInverter()->distribution((xmax-characteristic)/2.0*(-sqrt(3.0/5.0))+(xmax+characteristic)/2.0,primaryAbscissa[cellI],sigma[cellI])));*/
                        
                        
                        //trapeze
                        int n = 0;
                        scalar s = 0.0;
                        scalar olds=-1.0e30;
                        
                        do 
                        {
                            olds = s;
                            n+=1;
                            
                            if (n==21)
                            {
                                Info << "Too many steps" << endl;
                            }
                            int it=pow(2,(n-2));
                            scalar del = (xmax-characteristic)/it;
                            scalar x=characteristic+0.5*del;
                            scalar sum=0.0;
                            for (int j=1 ; j<=it; j++)
                            {
                                sum += pow(convectionModel_->characteristic(x,cellI),mI)*quadrature_.momentInverter()->distribution(x,primaryAbscissa[cellI],sigma[cellI]);
                                x += del;
                            }
                            s = 0.5*(s+(xmax-characteristic)*sum/it);
                    
                        }while(fabs(s-olds)>1.e-6*fabs(olds) && (s!=0.0 && olds!=0.0));
                            
                        moment[cellI] += node.primaryWeight()[cellI]*s;
                        
                        //Romberg
                        
                        /*scalar s[21]; //array of lenght 21
                        scalar h[21]; 
                        scalar ss = 0.0;
                        scalar dss = 0.0;
                        int n = 0;
                        
                        s[1]=0.0;

                        do 
                        {
                            n+=1;
                            
                            if (n==21)
                            {
                                Info << "Too many steps in Romberg" << endl;
                            }
                            for (int i=1 ; i<=n ; i++)
                            {
                                int it=pow(2,(i-2));
                                scalar del = (xmax-characteristic)/it;
                                scalar x=characteristic+0.5*del;
                                scalar sum=0.0;
                                for (int j=1 ; j<=it; j++)
                                {
                                    sum += pow(convectionModel_->characteristic(x,cellI),mI)*quadrature_.momentInverter()->distribution(x,primaryAbscissa[cellI],sigma[cellI]);
                                    x += del;
                                }
                                s[n] = 0.5*(s[n]+(xmax-characteristic)*sum/it);
                            }
                            if (n>=5)
                            {
                                scalar c[10];
                                scalar d[10];
                                int ns = 1;
                                scalar dif=h[n-4];
                                
                                for (int i=1; i<=5; i++)
                                {
                                    scalar dift = h[n-4];
                                    if (dift < dif)
                                    {
                                        ns=i;
                                        dif=dift;
                                    }
                                    c[i]=s[n-3];
                                    d[i]=s[n-3];
                                }
                                ss = s[n-4+ns];
                                ns -= 1;
                                for (int m=1; m<5 ; m++)
                                {
                                    for (int i=1; i<=5-m ; i++)
                                    {
                                        scalar ho=h[n-4+i];
                                        scalar hp=h[n-4+i+m];
                                        scalar w = c[i+1]-d[i];
                                        scalar den = w/(ho - hp);
                                        d[i]=hp*den;
                                        c[i]=ho*den;
                                    }
                                    if (2*ns<(5-m))
                                    {
                                        dss=c[ns+1];
                                    }
                                    else
                                    {
                                        dss=d[ns];
                                        ns=ns-1;
                                    }
                                    ss += dss;
                                }
                            }
                            s[n+1]=s[n];
                            h[n+1] = 0.25*h[n];
                        }while(dss<=1.e-6*ss);
                            
                        moment[cellI] += node.primaryWeight()[cellI]*ss;*/
                        
                    }
                    else if (primaryAbscissaFinal>0.0)
                    {
                        forAll(node.secondaryWeights(), sNodeI)
                        {
                            moment[cellI] += node.primaryWeight()[cellI]*node.secondaryWeights()[sNodeI][cellI]*pow(convectionModel_->characteristic(node.secondaryAbscissae()[sNodeI][cellI],cellI),mI);
                        }
                    }
                }
            }
        }
        //Info << quadrature_.moments() << endl;
    
        quadrature_.updateQuadrature();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::momentSource
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> mSource
    (
        new volScalarField
        (
            IOobject
            (
                "mSource",
                moment.mesh().time().timeName(),
                moment.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            moment.mesh(),
            dimensionedScalar
            (
                "mSource",
                moment.dimensions()/dimTime,
                0.0
            )
        )
    );
    mSource.ref() ==
        aggregationSource(moment)
      + breakupSource(moment)
      + nucleationModel_->nucleationSource(moment);

    return mSource;
}

void Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::solve
()
{
    univariatePDFTransportModel::solve();
}

// ************************************************************************* //
