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
    Foam::populationBalanceSubModels::oxidationModels::reduction

Description
    reduction oxidation model for soot formation.

SourceFiles
    reduction.C

\*---------------------------------------------------------------------------*/

#ifndef reduction_H
#define reduction_H

#include "oxidationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace oxidationModels
{

/*---------------------------------------------------------------------------*\
                        Class reduction Declaration
\*---------------------------------------------------------------------------*/

class reduction
:
    public oxidationModel
{
private:
        
    // Functions
        
        //-constant reduction
        scalar cRed( const label& cellI );
        

public:

        //- Runtime type information
        TypeName("reduction");


    // Constructors

        //- Construct from components
        reduction
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~reduction();


    // Member Functions

        //- Caracteristic oxidation equation
        virtual scalar characteristic
        (
            const scalar& abcissa,
            const label& cellI
        );
        
                //- Caracteristic oxidation equation
        virtual scalar characteristic
        (
            const label& cellI
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace oxidationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
