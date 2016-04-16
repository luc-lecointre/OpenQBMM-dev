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
    Foam::populationBalanceModels::aggregationKernels::FNP

Description
    Flash nano percipitation aggregation kernl.
    
    References
    \verbatim
        "Versuch Einer Mathematischen Theorie Der Koagulationskinetic Kolloider 
        Losunger"
        M X Smoluchowski
        Zeitschrift fur Physikalische Chemie
        Volume 92, Pages 129-142, 1917
    \endverbatim 
    
    \verbatim
        "A competitive aggregation model for Flash NanoPercipitation"
        Janine Chungyin Cheng, R.D. Vigil, R.O. Fox
        Journal of Colloid and Interface Science
        Volume 351, Pages 330-341, 2010
    \endverbatim 

SourceFiles
    FNP.C

\*---------------------------------------------------------------------------*/

#ifndef FNP_H
#define FNP_H

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class FNP Declaration
\*---------------------------------------------------------------------------*/

class FNP
:
    public aggregationKernel
{
//- Private variables

    //- Constant refrence to fluid thermo model
    const fluidThermo& flThermo_;

    //- Boltzmann's constant
    dimensionedScalar kb_;
    
    //- Avagradro's Number
    dimensionedScalar avagadro_;
    
    //- Number of unimers
    scalar nUnimers_;
    
    //- Volume of one polymer site
    dimensionedScalar polySite_;
    
    //- Length of polymers
    dimensionedScalar na_;
    dimensionedScalar nb_;
    
    //- Flory exponents
    dimensionedScalar nuA_;
    dimensionedScalar nuB_;
    
public:

        //- Runtime type information
        TypeName("FNP");


    // Constructors

        //- Construct from components
        FNP(const dictionary& dict);


    //- Destructor
    virtual ~FNP();


    // Member Functions

        //- Calcuation of beta kernel
        volScalarField beta(scalar p, scalar i);
        
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
