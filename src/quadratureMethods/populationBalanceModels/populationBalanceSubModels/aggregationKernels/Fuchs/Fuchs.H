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

Class
    Foam::populationBalanceSubModels::aggregationKernels::Fuchs

Description
    Fuchs aggregation kernel.
    
    Reference
    \verbatim
    ???
    \endverbatim    

SourceFiles
    Fuchs.C

\*---------------------------------------------------------------------------*/

#ifndef Fuchs_H
#define Fuchs_H

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class Fuchs Declaration
\*---------------------------------------------------------------------------*/

class Fuchs
:
    public aggregationKernel
{
private :
    
    //- Molecular mass of one carbon atom
    dimensionedScalar MCarbon_;
    
    //- fractal dimension depending of flame configuration
    scalar df_;
        
    //- mean free path of molecule
    dimensionedScalar lambdaf_;

    //- Number of carbon atoms in PAM species
    label nCarbonSoot_;
    
    //- density of soot
    dimensionedScalar rhoSoot_;
    
    //- abscissa nucleation
    dimensionedScalar abscissa0_;
    
    //functions
    tmp<volScalarField> Kn
    (
        const volScalarField& abscissa
    ) const;
    
    tmp<volScalarField> D
    (
        const volScalarField& abscissa
    ) const;
    
    tmp<volScalarField> dc
    (
        const volScalarField& abscissa
    ) const;
    
    tmp<volScalarField> velocity
    (
        const volScalarField& abscissa
    ) const;
    
    tmp<volScalarField> l
    (
        const volScalarField& abscissa
    ) const;
    
    tmp<volScalarField> g
    (
        const volScalarField& abscissa
    ) const;
    
public:

        //- Runtime type information
        TypeName("Fuchs");


    // Constructors

        //- Construct from components
        Fuchs(const dictionary& dict);


    //- Destructor
    virtual ~Fuchs();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> Ka
        (
            const volScalarField& abscissa1,
            const volScalarField& abscissa2
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

