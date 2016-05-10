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
    const surfaceScalarField& phi,
    const word& support
)
:
    PDFTransportModel(name, dict, mesh),
    name_(name),
    support_(support),
    interfaceModel_(dict.lookup("interfaceModel")),
    ode_(dict.lookup("ode")),
    momentInverter_(),
    ATol_(readScalar(dict.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(dict.subDict("odeCoeffs").lookup("RTol"))),
    fac_(readScalar(dict.subDict("odeCoeffs").lookup("fac"))),
    facMin_(readScalar(dict.subDict("odeCoeffs").lookup("facMin"))),
    facMax_(readScalar(dict.subDict("odeCoeffs").lookup("facMax"))),
    h_(facMin_*U.mesh().time().deltaT()),
    maxDeltaT_(false),
    quadrature_(name, mesh, support),
    U_(U),
    phi_(phi),
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
    
    forAll(moments_, mi)
    {
        moments_.set
        (
            mi,
            new volUnivariateMoment
            (
                name_,
                quadrature_.moments()[mi].cmptOrders(),
                nodes_,
                quadrature_.moments()[mi]
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
::updatePhysicalSpaceConvection()
{
    // Update interpolated nodes
    quadrature_.interpolateNodes();

    // Updated reconstructed moments
    quadrature_.momentsNei().update();
    quadrature_.momentsOwn().update();
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::univariatePDFTransportModel::physicalSpaceConvection
(
    const volUnivariateMoment& moment
)
{
    dimensionedScalar zeroPhi("zero", phi_.dimensions(), 0.0);

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
        quadrature_.momentsNei()[order]*min(phi_, zeroPhi)
      + quadrature_.momentsOwn()[order]*max(phi_, zeroPhi)
    );

    fvc::surfaceIntegrate(divMoment.ref(), mFlux);
    divMoment.ref().dimensions().reset(moment.dimensions()/dimTime);
    
    return divMoment;
}


void Foam::PDFTransportModels::univariatePDFTransportModel::updateQuadrature()
{
    const volScalarField& m0(moments_[0]);

    PtrList<extendedVolScalarNode>& nodes(nodes_());

    forAll(m0, celli)
    {
        univariateMomentSet 
            momentsToInvert(quadrature_.nMoments(), 0.0, "Gauss", support_);

        // Copying moment set from a cell to univariateMomentSet
        forAll(momentsToInvert, mi)
        {
            momentsToInvert[mi] = moments_[mi][celli];
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
            node.primaryWeight()[celli] = pWeights[pNodeI];
            node.primaryAbscissa()[celli] = pAbscissae[pNodeI];

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
                sWeightFields[sNodeI][celli] = sWeights[pNodeI][sNodeI];
                sAbscissaFields[sNodeI][celli] = sAbscissae[pNodeI][sNodeI];
            }

            // Copy sigma
            node.sigma()[celli] = momentInverter_->sigma();
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

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::univariatePDFTransportModel::interfaceSource
(
    const volUnivariateMoment& m
)
{
    tmp<volScalarField> iSource
    (
        new volScalarField
        (
            IOobject
            (
                "iSource",
                m.mesh().time().timeName(),
                m.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            m.mesh(),
            dimensionedScalar("iSource", m.dimensions()/dimTime, 0.0)
        )
    );
    
    if (!interfaceModel_)
    {
        return iSource;
    }
    
    volScalarField& interSrc = iSource.ref();
    
    word order;
    
    forAll(m.cmptOrders(), dimi)
    {
        order += Foam::name(m.cmptOrders()[dimi]);
    }
    
    PtrList<volScalarField> P(2);
    
    forAll(P, nodei)
    {
        P.set
        (
            nodei,
            volScalarField
            (
                m.mesh().lookupObject<volScalarField>
                (
                    "mixingWeights."
                  + Foam::name(nodei)
                )
            )
        );
    }            
        
    typedef incompressible::turbulenceModel icoTurbModel;

    const incompressible::turbulenceModel& turb =
        m.mesh().lookupObject<icoTurbModel>
        (
            icoTurbModel::propertiesName
        );

    volScalarField gamma = 2.0*turb.epsilon()/turb.k(); // Cphi = 2.0
    
    dimensionedScalar minDiff("min", m.dimensions(), 1e-4);
    
    if (name_ == "populationBalance2")
    {
        const volUnivariateMoment& m1 =
            m.mesh().lookupObject<volUnivariateMoment>
            (
                "moment."
                + order
                + ".populationBalance1"
            );
        
        const volUnivariateMoment& m2 = m;
        
        interSrc == gamma*P[0]*P[1]*(m1 - m2) 
          + turb.nut()/max(minDiff, m2 - m1)
           *(P[0]*(fvc::grad(m1) & fvc::grad(m1))
          + P[1]*(fvc::grad(m2) & fvc::grad(m2)));
    }
    else
    {
        const volUnivariateMoment& m1 = m;
        const volUnivariateMoment& m2 =
            m.mesh().lookupObject<volUnivariateMoment>
            (
                "moment."
                + order
                + ".populationBalance2"
            );
        interSrc == gamma*P[0]*P[1]*(m2 - m1)
          - turb.nut()/max(minDiff, m1 - m2)
           *(P[0]*(fvc::grad(m1) & fvc::grad(m1))
          + P[1]*(fvc::grad(m2) & fvc::grad(m2)));
    }
    iSource().write();
    return iSource;
}

void Foam::PDFTransportModels::univariatePDFTransportModel::solveMomentSource()
{
    Info<< "RK23-SSP: Solving for moment source terms" << endl;
    
    // Read current deltaT
    dimensionedScalar dt0 = U_.mesh().time().deltaT();
    
    //- Initialize rate change PtrLists
    PtrList<volScalarField> k1(quadrature_.nMoments());
    PtrList<volScalarField> k2(quadrature_.nMoments());
    PtrList<volScalarField> k3(quadrature_.nMoments());
    
    forAll(moments_, mi)
    {
        k1.set
        (
            mi,
            new volScalarField
            (
                IOobject
                (
                    "k1",
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar("k1", moments_[mi].dimensions(), 0.0)
            )
        );
        
        k2.set
        (
            mi,
            new volScalarField
            (
                IOobject
                (
                    "k2",
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar("k2", moments_[mi].dimensions(), 0.0)
            )
        );
        
        k3.set
        (
            mi,
            new volScalarField
            (
                IOobject
                (
                    "k3",
                    U_.mesh().time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                U_.mesh(),
                dimensionedScalar("k3", moments_[mi].dimensions(), 0.0)
            )
        );
        
        moments_[mi] == quadrature_.moments()[mi];
    }
    
    updateQuadrature();
    
    // create a volScalarField copy of moments, updates before each itteration
    PtrList<volScalarField> momentsOld(quadrature_.nMoments());
    
    forAll(moments_, mi)
    {
        momentsOld.set
        (
            mi,
            new volScalarField(quadrature_.moments()[mi])
        );
        
        moments_[mi] == quadrature_.moments()[mi];
    }
    
    updateQuadrature();
    
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
        forAll(moments_, mi)
        {
            momentsOld[mi] == moments_[mi];
        }

        nItt++;
            
        // Calculate k1 for all moments
        forAll(moments_, mi)
        {
            k1[mi] = h_*momentSourceODE(moments_[mi]);
            moments_[mi] == momentsOld[mi] + k1[mi];
        }
        
        updateQuadrature();
        
        // Calculate k2 for all moments
        forAll(moments_, mi)
        {
            k2[mi] = h_*momentSourceODE(moments_[mi]);
            moments_[mi] == momentsOld[mi] + (k1[mi] + k2[mi])/4.0;
        }
    
        updateQuadrature();
        
        // calculate k3 and new moments for all moments
        forAll(moments_, mi)
        {
            k3[mi] = h_*momentSourceODE(moments_[mi]);
            
            // Second order accurate, k3 only used for error estimation
            moments_[mi] == momentsOld[mi] + (k1[mi] + k2[mi] + 4.0*k3[mi])/6.0;
        }
        updateQuadrature();
        
        // Calculate error
        scalar sc = 0.0;
        scalar error = 0.0;
        
        forAll(moments_, mi)
        {
            sc = min(sc, ATol_
              + min(max(mag(momentsOld[mi]), mag(moments_[mi]))).value()*RTol_);
            
            error = max(error, max(mag(k1[mi] + k2[mi] - 2.0*k3[mi])).value());
        }
        
        scalar err = error/(3.0*sc);
        
        if (err == 0.0)
        {
            h_ = dt0 - dTime;
            maxDeltaT_ = true;
        }
        
        else
        {
            h_ = max(dt0*facMin_, min(h_*facMax_, dt0*fac_/pow(err, 1.0/3.0)));
        }
        
        if (dTime.value() >= dt0.value())
        {
            timeComplete = true;
        }
        
        if 
        (
            Foam::name(U_.mesh().time().deltaT().value()) == U_.mesh().time().timeName()
         && timeComplete == true
        )
        {
            maxDeltaT_ = false;
            h_ = dt0*facMin_;
        }
        
        // Write some stuff
        Info<< "Iteration " << nItt 
            << ", "
            << "Residual = " << err
            << ","
            << "Local time = " << dTime.value() << "s"
            << endl;
            
    }

    if (h_.value() == dt0.value())
    {
        maxDeltaT_ = true;
    }
    
    return;
}


void Foam::PDFTransportModels::univariatePDFTransportModel::solve()
{
    updatePhysicalSpaceConvection();
    
    if (ode_)
    {
        solveMomentSource();
    }
    
    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());

    // Solve moment transport equations
    forAll(quadrature_.moments(), momenti)
    {
        volUnivariateMoment& m = quadrature_.moments()[momenti];

        momentEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              + physicalSpaceConvection(m)
              - momentDiffusion(m)
              ==
                (moments_[momenti] - m)/U_.mesh().time().deltaT()
              + momentSource(m)
              + interfaceSource(m)
            )
        );
    }

    forAll(momentEqns, mEqni)
    {
        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }
        
    quadrature_.updateQuadrature();
    
}


// ************************************************************************* //
