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

#include "univariatePDFTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::univariatePDFTransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const volVectorField& U,
    const word support
)
:
    PDFTransportModel(name, dict, mesh),
    name_(name),
    support_(support),
    ode_(dict.lookup("ode")),
    quadrature_(name, mesh, support),
    momentInverter_(),
    deltaM_(quadrature_.nMoments()),
    k1_(quadrature_.nMoments()),
    k2_(quadrature_.nMoments()),
    k3_(quadrature_.nMoments()),
    ATol_(readScalar(dict.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(dict.subDict("odeCoeffs").lookup("RTol"))),
    fac_(readScalar(dict.subDict("odeCoeffs").lookup("fac"))),
    facMin_(readScalar(dict.subDict("odeCoeffs").lookup("facMin"))),
    facMax_(1.0),
    h_(facMin_*U.mesh().time().deltaT()),
    maxDeltaT_(false),
    U_(U),
    nodes_(),
    moments_
    (
        name_,
        quadrature_.nMoments(),
        nodes_,
        quadrature_.nDimensions(),
        quadrature_.moments().momentMap()
    )
{
    nodes_ = autoPtr<PtrList<extendedVolScalarNode> >
    (
        new PtrList<extendedVolScalarNode> 
        (
            quadrature_.lookup("nodes"),
            Foam::extendedVolScalarNode::iNew
            (
                name_,
                mesh_,
                quadrature_.moments()[0].dimensions(),
                quadrature_.moments()[1].dimensions()
                   /quadrature_.moments()[0].dimensions(),
                quadrature_.moments()[0].boundaryField().types()
            )
                
        )
    );
    
    forAll(moments_, mI)
    {
        moments_.set
        (
            mI,
            new volUnivariateMoment
            (
                name_,
                quadrature_.moments()[mI].cmptOrders(),
                nodes_,
                quadrature_.moments()[mI]
            )
        );
        
        deltaM_.set
        (
            mI,
            new volScalarField
            (
                IOobject
                (
                    "deltaM." + moments_[mI].name(),
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U.mesh(),
                dimensionedScalar
                (
                    "deltaM",
                    moments_[mI].dimensions()/dimTime,
                    0.0
                )
            )
        );
        
        k1_.set
        (
            mI,
            new volScalarField
            (
                IOobject
                (
                    "k1." + moments_[mI].name(),
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar("k1", moments_[mI].dimensions(), 0.0)
            )
        );
        
        k2_.set
        (
            mI,
            new volScalarField
            (
                IOobject
                (
                    "k2." + moments_[mI].name(),
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar("k2", moments_[mI].dimensions(), 0.0)
            )
        );
        
        k3_.set
        (
            mI,
            new volScalarField
            (
                IOobject
                (
                    "k3." + moments_[mI].name(),
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar("k3", moments_[mI].dimensions(), 0.0)
            )
        );
    }
    
    momentInverter_ = autoPtr<Foam::extendedMomentInversion>
    (
        Foam::extendedMomentInversion::New
        (
            quadrature_.subDict("extendedMomentInversionCoeff"),
            quadrature_.nMoments(),
            quadrature_.nodes()[0].nSecondaryNodes()
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::~univariatePDFTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::univariatePDFTransportModel
::updatePhysicalSpaceConvection
(
    surfaceScalarField& phiOwn,
    surfaceScalarField& phiNei
)
{
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("nei", dimless, -1.0)
    );

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("own", dimless, 1.0)
    );

    phiOwn = fvc::interpolate(U_, own, "reconstruct(U)") & mesh_.Sf();
    phiNei = fvc::interpolate(U_, nei, "reconstruct(U)") & mesh_.Sf();

    // Update interpolated nodes
    quadrature_.interpolateNodes();

    // Updated reconstructed moments
    quadrature_.momentsNei().update();
    quadrature_.momentsOwn().update();
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::univariatePDFTransportModel::physicalSpaceConvection
(
    const volUnivariateMoment& moment,
    const surfaceScalarField& phiOwn,
    const surfaceScalarField& phiNei
)
{
    dimensionedScalar zeroPhi("zero", phiNei.dimensions(), 0.0);

    tmp<volScalarField> divMoment
    (
        new volScalarField
        (
            IOobject
            (
                "divMoment",
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

    label order = moment.order();

    surfaceScalarField mFlux
    (
        quadrature_.momentsNei()[order]*min(phiNei, zeroPhi)
      + quadrature_.momentsOwn()[order]*max(phiOwn, zeroPhi)
    );

    fvc::surfaceIntegrate(divMoment.ref(), mFlux);
    divMoment.ref().dimensions().reset(moment.dimensions()/dimTime);

    return divMoment;
}


void Foam::PDFTransportModels::univariatePDFTransportModel::updateQuadrature()
{
    const volScalarField& m0(moments_[0]);

    PtrList<extendedVolScalarNode>& nodes(nodes_());

    forAll(m0, cellI)
    {
        univariateMomentSet momentsToInvert(quadrature_.nMoments(), 0.0, support_);

        // Copying moment set from a cell to univariateMomentSet
        forAll(momentsToInvert, mI)
        {
            momentsToInvert[mI] = moments_[mI][cellI];
        }

        // Inverting moments and updating secondary quadrature
        momentInverter_->invert(momentsToInvert);

        // Recovering primary weights and abscissae from moment inverter
        const scalarDiagonalMatrix& pWeights(momentInverter_->primaryWeights());

        const scalarDiagonalMatrix& pAbscissae
        (
            momentInverter_->primaryAbscissae()
        );

        // Copying to fields
        for (label pNodeI = 0; pNodeI < quadrature_.nodes().size(); pNodeI++)
        {
            extendedVolScalarNode& node(nodes[pNodeI]);

            // Copy primary node
            node.primaryWeight()[cellI] = pWeights[pNodeI];
            node.primaryAbscissa()[cellI] = pAbscissae[pNodeI];

            // Copy secondary nodes
            PtrList<volScalarField>& sWeightFields(node.secondaryWeights());
            PtrList<volScalarField>& sAbscissaFields(node.secondaryAbscissae());

            const scalarRectangularMatrix& sWeights
            (
                momentInverter_->secondaryWeights()
            );

            const scalarRectangularMatrix& sAbscissae
            (
                momentInverter_->secondaryAbscissae()
            );

            for (label sNodeI = 0; sNodeI < quadrature_.nodes()[0].nSecondaryNodes(); sNodeI++)
            {
                sWeightFields[sNodeI][cellI] = sWeights[pNodeI][sNodeI];
                sAbscissaFields[sNodeI][cellI] = sAbscissae[pNodeI][sNodeI];
            }

            // Copy sigma
            node.sigma()[cellI] = momentInverter_->sigma();
        }
    }

    // Updating boundary conditions
    forAll(nodes, pNodeI)
    {
        extendedVolScalarNode& pNode(nodes[pNodeI]);

        pNode.primaryWeight().correctBoundaryConditions();
        pNode.primaryAbscissa().correctBoundaryConditions();
        pNode.sigma().correctBoundaryConditions();

        for (label sNodeI = 0; sNodeI < quadrature_.nodes()[0].nSecondaryNodes(); sNodeI++)
        {
            pNode.secondaryWeights()[sNodeI].correctBoundaryConditions();
            pNode.secondaryAbscissae()[sNodeI].correctBoundaryConditions();
        }
    }

    moments_.update();
}

void Foam::PDFTransportModels::univariatePDFTransportModel::solveMomentSource()
{
    Info<< quadrature_nMoments() << endl;
    if (!ode_)
    {
        forAll(moments_, mI)
        {
            moments_[mI] == quadrature_.moments()[mI];
        }
        
        updateQuadrature();
        
        forAll(moments_, mI)
        {
            deltaM_[mI] = momentSource(quadrature_.moments()[mI]);
        }
        return;
    }
    
    // create a volScalarField copy of moments, updates before each itteration
    PtrList<volScalarField> momentsOld(quadrature_.nMoments());
    
    forAll(moments_, mI)
    {
        momentsOld.set
        (
            mI,
            new volScalarField(quadrature_.moments()[mI])
        );
        
        moments_[mI] == quadrature_.moments()[mI];
    }
    
    updateQuadrature();
    
    // Read current deltaT
    dimensionedScalar dt0 = U_.mesh().time().deltaT();
    
    if (h_ > dt0)
    {
        maxDeltaT_ = true;
    }
    
    if (maxDeltaT_)
    {
        h_ = dt0;
    }
    
    dimensionedScalar dTime("dTime", dimTime, 0.0);
    
    label nItt = 0;
    bool timeComplete = false;
    
    while (!timeComplete)
    {
        if (dTime + h_ > dt0)
        {
            h_ = dt0 - dTime;
        }
        
        dTime += h_;
        
        // set original moments for current itteration
        forAll(moments_, mI)
        {
            momentsOld[mI] == moments_[mI];
        }

        nItt++;
            
        // Calculate k1 for all moments
        forAll(moments_, mI)
        {
            k1_[mI] = h_*momentSource(moments_[mI]);
            volScalarField mZero = momentsOld[mI] + k1_[mI];
            
            moments_[mI] == mZero;
        }
        
        updateQuadrature();
        
        // Calculate k2 for all moments
        forAll(moments_, mI)
        {
            k2_[mI] = h_*momentSource(moments_[mI]);
            volScalarField mWhole = momentsOld[mI] + (k1_[mI] + k2_[mI])/4.0;
            
            moments_[mI] == mWhole;
        }
    
        updateQuadrature();
        
        // calculate k3 and new moments for all moments
        forAll(moments_, mI)
        {
            k3_[mI] = h_*momentSource(moments_[mI]);
            
            moments_[mI] == momentsOld[mI]
                + (k1_[mI] + k2_[mI] + 4.0*k3_[mI])/6.0;
        }
        
        updateQuadrature();        
        
        // Calculate error
        scalar err = 0.0;
        
        forAll(moments_, mI)
        {
            scalar sc = 
                ATol_ + min(max(mag(momentsOld[mI]), mag(moments_[mI]))).value()*RTol_;
            err += sqr(max((k1_[mI] + k2_[mI] - 2.0*k3_[mI])/(3.0*sc)).value());
        }
        
        err = sqrt(err/quadrature_.nMoments());
        
        if (err == 0.0)
        {
            h_ = dt0 - dTime;
            maxDeltaT_ = true;
        }
        
        else
        {
            h_ = dt0*min(facMax_, max(facMin_, fac_/pow(err, 1.0/3.0)));
        }
        
        if (dTime.value() >= dt0.value())
        {
            timeComplete = true;
        }
    }

    if (h_.value() == dt0.value())
    {
        maxDeltaT_ = true;
    }
    
    forAll(moments_, mI)
    {
        deltaM_[mI] == (moments_[mI] - quadrature_.moments()[mI])/dTime;
    }

    Info<< "Solving for moment source terms, No itterations " << nItt << endl;
    return;
}


void Foam::PDFTransportModels::univariatePDFTransportModel::solve()
{
    surfaceScalarField phiOwn("phiOwn", fvc::interpolate(U_) & mesh_.Sf());
    surfaceScalarField phiNei("phiNei", phiOwn);
    updatePhysicalSpaceConvection(phiOwn, phiNei);
    
    solveMomentSource();
    
    // Solve moment transport equations
    forAll(quadrature_.moments(), mI)
    {
        volUnivariateMoment& m = quadrature_.moments()[mI];

        fvScalarMatrix momentEqn
        (
            fvm::ddt(m)
          + physicalSpaceConvection(m, phiOwn, phiNei)
          - momentDiffusion(m)
          ==
            deltaM_[mI]
          + phaseSpaceConvection(m)
        );

        momentEqn.relax();
        momentEqn.solve();

    }

    quadrature_.updateQuadrature();
}


// ************************************************************************* //
