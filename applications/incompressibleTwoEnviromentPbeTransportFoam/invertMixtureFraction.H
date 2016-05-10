PtrList<volScalarField> mixingWeights(2);

forAll(mixingWeights, nodei)
{
    mixingWeights.set
    (
        nodei,
        new volScalarField
        (
            IOobject
            (
                "mixingWeights." + Foam::name(nodei),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("mixingWeight", dimless, 0.0)
        )
    );
}

{
    PtrList<volScalarField> mixingMoments(4);

    forAll(mixingMoments, mi)
    {
        mixingMoments.set
        (
            mi,
            new volScalarField
            (
                IOobject
                (
                    "moment." + Foam::name(mi) + ".mixing",
                    "0",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            )
        );
    }

    forAll(mixingMoments[0], celli)
    {
        univariateMomentSet momentsToInvert(4, 0, "RPlus");
        
        forAll(mixingMoments, mi)
        {
            momentsToInvert[mi] = mixingMoments[mi][celli];
        }
        
        momentsToInvert.invert();
        
        forAll(momentsToInvert.weights(), nodei)
        {
            mixingWeights[nodei][celli] = momentsToInvert.weights()[nodei];
        }
    }
}

volScalarField mixtureFraction
(
    IOobject
    (
        "mixtureFraction",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
