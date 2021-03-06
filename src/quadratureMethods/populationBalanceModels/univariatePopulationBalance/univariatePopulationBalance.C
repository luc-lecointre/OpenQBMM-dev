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
    univariatePDFTransportModel(name, dict, U.mesh(), U, "RPlus"),
    populationBalanceModel(name, dict, U, phi),
    name_(name),
    aggregation_(dict.lookup("aggregation")),
    breakup_(dict.lookup("breakup")),
    growth_(dict.lookup("growth")),
    oxidation_(dict.lookup("oxidation")),
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
    growthModel_
    (
        Foam::populationBalanceSubModels::growthModel::New
        (
            dict.subDict("growthModel"),
            U.mesh()
        )
    ),
    oxidationModel_
    (
        Foam::populationBalanceSubModels::oxidationModel::New
        (
            dict.subDict("oxidationModel"),
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
                        pWeight1*sWeight1*
                        (
                            pWeight2*sWeight2*
                            (
                                0.5*pow // Birth
                                (
                                    pow3(sAbscissa1) + pow3(sAbscissa2),
                                    order/3.0
                                )
                              - pow(sAbscissa1, order)
                            )*aggregationKernel_->Ka(sAbscissa1, sAbscissa2)
                        );

                    aggregationSource.dimensions().reset
                    (
                        aggInnerSum().dimensions()
                    );

                    aggregationSource == aggregationSource + aggInnerSum();
                }
            }
        }
    }
    
    Info << "aggregationSource : "<< aggregationSource << endl;
    
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

    forAll(quadrature_.nodes(), pNodeI)
    {
        const extendedVolScalarNode& node = quadrature_.nodes()[pNodeI];

        forAll(node.secondaryWeights(), sNodeI)
        {
            tmp<volScalarField> bSrc = node.primaryWeight()
                *node.secondaryWeights()[sNodeI]
                *breakupKernel_->Kb(node.secondaryAbscissae()[sNodeI])
                *(
                    daughterDistribution_->mD                      //Birth
                    (
                        order,
                        node.secondaryAbscissae()[sNodeI]
                    )
                  - pow(node.secondaryAbscissae()[sNodeI], order)   //Death
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

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::phaseSpaceConvection
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> gSource
    (
        new volScalarField
        (
            IOobject
            (
                "gSource",
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

    if (!growth_)
    {
        gSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return gSource;
    }

    label order = moment.order();

    if (order < 1)
    {
        gSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return gSource;
    }

    volScalarField& growthSource = gSource.ref();

    forAll(quadrature_.nodes(), pNodeI)
    {
        const extendedVolScalarNode& node = quadrature_.nodes()[pNodeI];

        forAll(node.secondaryWeights(), sNodeI)
        {
            tmp<volScalarField> gSrc = node.primaryWeight()
                *node.secondaryWeights()[sNodeI]
                *growthModel_->Kg(node.secondaryAbscissae()[sNodeI])
                *order*pow(node.secondaryAbscissae()[sNodeI],order-1);

            growthSource.dimensions().reset(gSrc().dimensions());
            growthSource == growthSource + gSrc;
        }
    }

    return gSource;
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance::oxidationSpaceConvection
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> oxidationMoment
    (
        new volScalarField
        (
            IOobject
            (
                "oxiMoment",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless , 0.0)
        )
    );
    
    if (!oxidation_)
    {
        oxidationMoment.ref().dimensions().reset(moment.dimensions()/dimTime);

        return oxidationMoment;
    }

    label order = moment.order();

    if (order < 1)
    {
        oxidationMoment.ref().dimensions().reset(moment.dimensions()/dimTime);

        return oxidationMoment;
    }
    
    volScalarField& oxiMoment = oxidationMoment.ref();
    
    label N2 = quadrature_.nMoments()+1;
    
    forAll(moment,cellI)
    {
        univariateMomentSet mSet(N2,0,"RPlus");
        
        forAll(mSet,mI)
        {
            forAll(quadrature_.nodes(), pNodeI)
            {
                const extendedVolScalarNode& node = quadrature_.nodes()[pNodeI];

                forAll(node.secondaryWeights(), sNodeI)
                {
                    scalar m = node.primaryWeight()[cellI]
                        *node.secondaryWeights()[sNodeI][cellI]
                        *pow(node.secondaryAbscissae()[sNodeI][cellI]+oxidationModel_->characteristic(cellI),mI);
                    
                    mSet[mI] = mSet[mI] + m;
                }
            }
        }
        mSet.invert();
        
        forAll (mSet.weights(),I)
        {
            scalar oxi = mSet.weights()[I]
                *pow(oxidationModel_->characteristic(mSet.abscissae()[I],cellI),order)/U_.mesh().time().deltaT().value();
        
            oxiMoment[cellI] = oxiMoment[cellI] + oxi;
        }
    }
    
    oxidationMoment.ref().dimensions().reset(moment.dimensions()/dimTime);
    
    Info << "oxidation dimension :" << oxidationMoment.ref().dimensions() << endl;
    
    return oxidationMoment;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::momentSource
(
    const volUnivariateMoment& moment
)
{
    tmp<fvScalarMatrix> mSource
    (
        new fvScalarMatrix
        (
            moment,
            moment.dimensions()*dimVol/dimTime
        )
    );
    
    //Info << "ordre " << moment.order() << endl;
    
    mSource.ref() +=
        aggregationSource(moment) + breakupSource(moment)
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
