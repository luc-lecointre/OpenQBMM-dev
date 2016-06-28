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
::nucleationSource
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> nSource
    (
        new volScalarField
        (
            IOobject
            (
                "nSource",
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
    
    if (!nucleation_)
    {
        nSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return nSource;
    }
    
    volScalarField V = nSource.ref();
    
    forAll(V, cellI){
        V[cellI]= U_.mesh().V()[cellI];
    }
    
    volScalarField dConcentration = nSource.ref();
    
    forAll(nodes_(), pNodeI)
    {
        const extendedVolScalarNode& node = nodes_()[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            tmp<volScalarField> dSrc = node.primaryWeight()
                *node.secondaryWeights()[sNodei]
                *nucleationModel_->beta(node.secondaryAbscissae()[sNodei]);
                
            dConcentration.dimensions().reset(dSrc().dimensions());
            dConcentration == dConcentration + dSrc;
        }
    }
    
    volScalarField delta = sqr(dConcentration)+4.0*nucleationModel_->betaN()*nucleationModel_->betaPAH();
    
    volScalarField dimerConcentration = -(dConcentration-sqrt(delta))/(2*nucleationModel_->betaN());
    
    label order = moment.order();

    volScalarField& nucleationSource = nSource.ref();
    
    tmp<volScalarField> nSrc = 1.0/2.0*nucleationModel_->betaN()*sqr(dimerConcentration)*pow(nucleationModel_->xiNuc(),order)*Foam::constant::physicoChemical::NA*V;
    
    nucleationSource.dimensions().reset(nSrc().dimensions());
    nucleationSource == nSrc;
    
    forAll(nodes_(), pNodeI)
    {
        const extendedVolScalarNode& node = nodes_()[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            tmp<volScalarField> cSrc = 1.0/2.0*node.primaryWeight()
                *node.secondaryWeights()[sNodei]
                *nucleationModel_->xiCond(node.secondaryAbscissae()[sNodei], order)
                *Foam::constant::physicoChemical::NA*V
                *dimerConcentration
                *nucleationModel_->beta(node.secondaryAbscissae()[sNodei]);
            
            nucleationSource.dimensions().reset(cSrc().dimensions());
            nucleationSource == nucleationSource + cSrc;
        }
    }
    
    return nSource;
}

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
        //Info << quadrature_.moments() << endl;
        
        aSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return aSource;
    }
    
    label order = moment.order();

    volScalarField& aggregationSource = aSource.ref();

    forAll(nodes_(), pNode1I)
    {
        const extendedVolScalarNode& node1 = nodes_()[pNode1I];
    
        const volScalarField& pWeight1 = node1.primaryWeight();
        
        //Info << "weight first order : " << pWeight1.dimensions() << endl;

        forAll(node1.secondaryWeights(), sNode1I)
        {

            const volScalarField& sWeight1 = node1.secondaryWeights()[sNode1I];
            
            //Info << "weight second order : " << sWeight1.dimensions() << endl;
                        
            const volScalarField& sAbscissa1
                = node1.secondaryAbscissae()[sNode1I];
            
            //Info << "abscissa" << sAbscissa1 << endl;

            forAll(nodes_(), pNode2I)
            {
                const extendedVolScalarNode& node2 = nodes_()[pNode2I];

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
    
    //Info << "aggregationSource : "<< max(aggregationSource)*U_.mesh().time().deltaT().value() << endl;
    
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
        //Info << moments_ << endl;
        /*volUnivariateMoment moment = moments_[0];
        volScalarField abscissa = nodes_[0].primaryAbscissa();*/
        
        volScalarField cRed = convectionModel_->cRed();
        volScalarField cGro = convectionModel_->cGro();
        
        //Info << cGro - cRed << endl;
        
        forAll(moments_, mI)
        {
            volUnivariateMoment& moment = moments_[mI];
            
            volScalarField conventionSource
            (
                IOobject
                (
                    "cSource",
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar("zero", moment.dimensions(), 0.0)
            );
            
            //Info << moment << endl;
        

            forAll(nodes_(), pNodeI)
            {
                const extendedVolScalarNode& node = nodes_()[pNodeI];
                volScalarField primaryAbscissa = node.primaryAbscissa();
                volScalarField sigma = node.sigma();
                
                //Info << "sigma : " << max(sigma) << endl;
                //Info << "primaryAbscissa : " << primaryAbscissa << endl;
                
                forAll(conventionSource, cellI)
                {
                    scalar characteristic = convectionModel_->characteristic(cRed[cellI], cGro[cellI]);
                    
                    //Info << "sigma : " << sigma[cellI] << endl;
                    
                    //Info <<"cellI :" << cellI << " -> characteristic : " << characteristic << endl;
                    
                    scalar primaryAbscissaFinal = convectionModel_->characteristic(primaryAbscissa[cellI], cRed[cellI], cGro[cellI]);
                    
                    //Info << cellI << " : " << primaryAbscissa[cellI] << " -> " << primaryAbscissaFinal << endl; 
                    
                    scalar xmax=characteristic+10.0;
            
                    if (sigma[cellI]!=0 && characteristic!=0.0)
                    {
                        while (momentInverter_->distribution(xmax,primaryAbscissa[cellI],sigma[cellI])>1.0e-8)
                        {
                            xmax+=10.0;
                        }
                        
                        //Info << cellI << " -> xmax : " << xmax << endl;
                        
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
                                sum += pow(convectionModel_->characteristic(x,cRed[cellI], cGro[cellI]),mI)*momentInverter_->distribution(x,primaryAbscissa[cellI],sigma[cellI]);
                                x += del;
                            }
                            s = 0.5*(s+(xmax-characteristic)*sum/it);
                    
                        }while(fabs(s-olds)>1.e-6*fabs(olds) && (s!=0.0 && olds!=0.0));
                            
                        conventionSource[cellI] += node.primaryWeight()[cellI]*s;
                        
                    }
                    else if (primaryAbscissaFinal>0.0)
                    {
                        //Info << "point 2 : "<< cellI << endl;
                        forAll(node.secondaryWeights(), sNodeI)
                        {
                            conventionSource[cellI] += node.primaryWeight()[cellI]*node.secondaryWeights()[sNodeI][cellI]*pow(convectionModel_->characteristic(node.secondaryAbscissae()[sNodeI][cellI],cRed[cellI], cGro[cellI]),mI);
                        }
                    }
                }
            }
            
            forAll(moment, cellI)
            {
                //moments_[mI][cellI]=conventionSource[cellI];
                quadrature_.moments()[mI][cellI]=conventionSource[cellI];
            }
        }
    
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
      //+ nucleationSource(moment);
      
    //Info << "aggregation : " << aggregationSource(moment) << endl;
    //Info << "breakup : " << max(breakupSource(moment))*U_.mesh().time().deltaT().value() << endl;
    //Info << "nucleation : " << nucleationModel_->nucleationSource(moment) << endl;

    return mSource;
}

void Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::solve
()
{
    univariatePDFTransportModel::solve();
}

// ************************************************************************* //
