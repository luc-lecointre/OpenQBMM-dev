/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::convectionModels::growthReduction

Description
    growthReduction convection model for soot formation.

SourceFiles
    growthReduction.C

\*---------------------------------------------------------------------------*/

#ifndef growthReduction_H
#define growthReduction_H

#include "convectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace convectionModels
{

/*---------------------------------------------------------------------------*\
                        Class growthReduction Declaration
\*---------------------------------------------------------------------------*/

class growthReduction
:
    public convectionModel
{
private:
    // Constants
    
        //- Carbon atom Volume
        dimensionedScalar VCarbon_;
        
        //- concentration of active sites per units of surface
        dimensionedScalar Xi_;
        
    // Functions
        
        //- Arrhenius law
        volScalarField Arrhenius
        (
            const scalar& A,
            const scalar& n,
            const scalar& E,
            const label& nReactifs
        );
        
        scalar Arrhenius
        (
            const scalar& A,
            const scalar& n,
            const scalar& E,
            const scalar& T
        );
        
        volScalarField r();
        

public:

        //- Runtime type information
        TypeName("growthReduction");


    // Constructors

        //- Construct from components
        growthReduction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~growthReduction();


    // Member Functions
    
        //-constant growthReduction
        virtual volScalarField cRed();
        virtual volScalarField cGro();

        //- Caracteristic convection equation
        virtual scalar characteristic
        (
            const scalar& abscissa,
            const scalar& cRed,
            const scalar& cGro
        );
        
        //- Caracteristic convection equation
        virtual scalar characteristic
        (
            const scalar& cRed,
            const scalar& cGro
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace oxidationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
