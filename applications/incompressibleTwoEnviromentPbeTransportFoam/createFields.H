Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info<< "Reading populationBalanceProperties\n" << endl;

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

#include "invertMixtureFraction.H"

autoPtr<populationBalanceModel> populationBalance1
(
    populationBalanceModel::New
    (
        "populationBalance1", populationBalanceProperties, U, phi
    )
);

autoPtr<populationBalanceModel> populationBalance2
(
    populationBalanceModel::New
    (
        "populationBalance2", populationBalanceProperties, U, phi
    )
);

PtrList<volScalarField> moments(5);

forAll(moments, mi)
{
    const volUnivariateMoment& m = 
        mesh.lookupObject<volUnivariateMoment>
        (
            "moment."
          + Foam::name(mi)
          + ".populationBalance1"
        );
        
    moments.set
    (
        mi,
        new volScalarField
        (
            IOobject
            (
                "moment." + Foam::name(mi) + ".global",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "moments",
                m.dimensions(),
                0.0
            )
        )
    );
}
