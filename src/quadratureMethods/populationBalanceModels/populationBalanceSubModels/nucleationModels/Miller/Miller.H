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
    Foam::populationBalanceSubModels::nucleationModels::Miller

Description
    Miller nucleation model for soot formation.

SourceFiles
    Miller.C

\*---------------------------------------------------------------------------*/

#ifndef Miller_H
#define Miller_H

#include "nucleationModel.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                        Class Miller Declaration
\*---------------------------------------------------------------------------*/

class Miller
:
    public nucleationModel
{
private:

    // Private data

        //- Molecular mass of one carbon atom
        dimensionedScalar MCarbon_;
        
        //name of the PAH
        string PAH_;

        //- Number of carbon atoms in PAH species
        string nCarbonPAH_;
        
        //sticking coefficient
        string gamma_;

        //- Soot density
        dimensionedScalar rhoSoot_;
        
        //-volume of the mesh
        volScalarField V_;
        
        List<word> tokens;
        List<label> nCarbon;
        List<scalar> gamma;
        
        dimensionedScalar volume(const scalar& nC);
        volScalarField volume(const volScalarField& nC);
        volScalarField Kfm(const scalar& nC);
        volScalarField Kfm(const volScalarField& nC);
        tmp<volScalarField> dimerConcentration(const scalar& nC);
        volScalarField nCmoy();

public:

        //- Runtime type information
        TypeName("Miller");


    // Constructors

        //- Construct from components
        Miller
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~Miller();


    // Member Functions

        //collision models
        virtual tmp<volScalarField> beta
        (
            const volScalarField& abscissa
        );
        
        virtual tmp<volScalarField> betaN();
        
        virtual tmp<volScalarField> betaPAH();
        
        virtual tmp<volScalarField> xiNuc();
        
        virtual tmp<volScalarField> xiCond
        (
            const volScalarField& abscissa,
            const label& order
        );
        
        //- Nucleation rate
        virtual tmp<volScalarField> nucleationSource
        (
            const volUnivariateMoment& moment
        ) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
