/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
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

Application
    Test-UnivariateMomentSet

Description
    Test univariateMomentSet class and methods.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "scalarMatrices.H"
#include "univariateMomentSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info << "Testing univariateMomentSet\n" << endl;

<<<<<<< HEAD
    label nMoments = 2;
=======
    label nMoments = 5;
>>>>>>> 732c5e01c23d931aeea9dd093ce6d7f6b0e3c5bf

    scalarDiagonalMatrix m(nMoments, 0.0);
    // Computing integer moments of a log-normal function
//     scalarDiagonalMatrix m(nMoments, 1.0);
//     scalar mu = 0.0;
//     scalar sigma = 0.25;
//
<<<<<<< HEAD
//     for (label momentI = 0; momentI < nMoments; momentI++)
//     {
//         m[momentI] = Foam::exp(momentI*mu + Foam::sqr(momentI*sigma)/2.0);
=======
//     for (label momenti = 0; momenti < nMoments; momenti++)
//     {
//         m[momenti] = Foam::exp(momenti*mu + Foam::sqr(momenti*sigma)/2.0);
>>>>>>> 732c5e01c23d931aeea9dd093ce6d7f6b0e3c5bf
//     }

    // Computing integer moments of a gaussian function
//     scalar mu = 0.0;
//     scalar sigma = 0.25;


//     m[0] = 1.0;
//
<<<<<<< HEAD
//     for (label momentI = 2; momentI < nMoments; momentI = momentI + 2)
//     {
//         m[momentI] = pow(sigma, momentI)*Foam::factorial(momentI-1);
=======
//     for (label momenti = 2; momenti < nMoments; momenti = momenti + 2)
//     {
//         m[momenti] = pow(sigma, momenti)*Foam::factorial(momenti-1);
>>>>>>> 732c5e01c23d931aeea9dd093ce6d7f6b0e3c5bf
//     }

    // Computing integer moments of sum of two beta function
//     scalar mu1 = 0.5;
//     scalar sigma1 = 0.3;
//     scalar mu2 = 0.5;
//     scalar sigma2 = 0.3;
//
//     m[0] = 1.0;
//
<<<<<<< HEAD
//     for (label momentI = 1; momentI < nMoments; momentI++)
//     {
//         m[momentI] = (((mu1 + (momentI - 1.0)*sigma1)*m[momentI-1]
//             /(1.0 + (momentI - 1.0)*sigma1))
//             +((mu2 + (momentI - 1.0)*sigma2)*m[momentI-1]
//             /(1.0 + (momentI - 1.0)*sigma2)))/2;
//     }

    // Computing integer moments of a beta function
//     scalar mu = 0.5;
//     scalar sigma = 0.3;
//
//     m[0] = 1.0;
//
//     for (label momentI = 1; momentI < nMoments; momentI++)
//     {
//         m[momentI] = (mu + (momentI - 1.0)*sigma)*m[momentI-1]
//             /(1.0 + (momentI - 1.0)*sigma);
//     }

    m[0] = 1.0;
    m[1] = 0.0;

    word support = "R";

    Info << "Support: " << support << endl;
=======
//     for (label momenti = 1; momenti < nMoments; momenti++)
//     {
//         m[momenti] = (((mu1 + (momenti - 1.0)*sigma1)*m[momenti-1]
//             /(1.0 + (momenti - 1.0)*sigma1))
//             +((mu2 + (momenti - 1.0)*sigma2)*m[momenti-1]
//             /(1.0 + (momenti - 1.0)*sigma2)))/2;
//     }

// Computing integer moments of a beta function
    scalar mu = 0.5;
    scalar sigma = 0.3;

    m[0] = 1.0;

    for (label momenti = 1; momenti < nMoments; momenti++)
    {
        m[momenti] = (mu + (momenti - 1.0)*sigma)*m[momenti-1]
            /(1.0 + (momenti - 1.0)*sigma);
    }



//    m[0] = 1.0;
//    m[1] = 0.0;

//     m[0] = 1;
//     m[1] = 0.000205634192732;
//     m[2] = 4.25189233395e-08;
//     m[3] = 3.63331214177e-10;
    //m[4] = 5.0e-11;
    //m[5] = 1.0e-16;
    //m[6] = 2.0e-22;

    word support = "01";
    word quadratureType = "GaussRadau";

    Info << "Support: " << support << endl;
    Info << "Quadrature type: " << quadratureType << endl;
>>>>>>> 732c5e01c23d931aeea9dd093ce6d7f6b0e3c5bf

    Info << setprecision(16);
    Info << "\nInput moments\n" << endl;

<<<<<<< HEAD
    for (label momentI = 0; momentI < nMoments; momentI++)
    {
        Info << "Moment " << momentI << " = " << m[momentI] << endl;
    }

    univariateMomentSet moments(m, support);

    Info << "\nStored moments\n" << endl;

    forAll(moments, momentI)
    {
        Info << "Moment " << momentI << " = " << moments[momentI] << endl;
=======
    for (label momenti = 0; momenti < nMoments; momenti++)
    {
        Info << "Moment " << momenti << " = " << m[momenti] << endl;
    }

    univariateMomentSet moments(m, quadratureType, support);

    Info << "\nStored moments\n" << endl;

    forAll(moments, momenti)
    {
        Info << "Moment " << momenti << " = " << moments[momenti] << endl;
>>>>>>> 732c5e01c23d931aeea9dd093ce6d7f6b0e3c5bf
    }

    moments.invert();

    if (moments.isFullyRealizable())
    {
	Info << "\nThe full set of moments is realizable.\n" << endl ;
    }
    else if (moments.isSubsetRealizable())
    {
        Info << "\nThe full set of moments is not realizable.\n" << endl
             << "The number of realizable moments is "
             << moments.nRealizableMoments() << "\n" << endl;
    }
    else
    {
        Info << "\nThe moment set is not realizable.\n" << endl;
    }

    Info << "The number of invertible moments is "
	 << moments.nInvertibleMoments() << "\n" << endl;

    scalarDiagonalMatrix weights(moments.weights());
    scalarDiagonalMatrix abscissae(moments.abscissae());

    Info << "Weights and abscissae:\n" << endl;

<<<<<<< HEAD
    for (label nodeI = 0; nodeI < moments.nNodes(); nodeI++)
    {
	Info << "Node " << nodeI
             << " Weight: " << weights[nodeI]
             << " Abscissa: " << abscissae[nodeI] << endl;
=======
    for (label nodei = 0; nodei < moments.nNodes(); nodei++)
    {
	Info << "Node " << nodei
             << " Weight: " << weights[nodei]
             << " Abscissa: " << abscissae[nodei] << endl;
>>>>>>> 732c5e01c23d931aeea9dd093ce6d7f6b0e3c5bf
    }

    moments.update();

    Info << "\nMoments computed from quadrature\n" << endl;

<<<<<<< HEAD
    for (label momentI = 0; momentI < nMoments; momentI++)
    {
        Info << "Moment " << momentI << " = " << moments[momentI] << endl;
=======
    for (label momenti = 0; momenti < nMoments; momenti++)
    {
        Info << "Moment " << momenti << " = " << moments[momenti] << endl;
>>>>>>> 732c5e01c23d931aeea9dd093ce6d7f6b0e3c5bf
    }

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
