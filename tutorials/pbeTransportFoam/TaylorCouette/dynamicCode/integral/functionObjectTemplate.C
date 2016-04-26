/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "functionObjectTemplate.H"
#include "Time.H"
#include "fvCFD.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(integralFunctionObject, 0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const objectRegistry& integralFunctionObject::obr() const
{
    return obr_;
}


const fvMesh& integralFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

integralFunctionObject::integralFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool
)
:
    name_(name),
    obr_(obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

integralFunctionObject::~integralFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void integralFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read integral sha1: 3811ea68234a05845d6b0899299560d4ea54b3ad\n";
    }

//{{{ begin code
    
//}}} end code
}


void integralFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute integral sha1: 3811ea68234a05845d6b0899299560d4ea54b3ad\n";
    }

//{{{ begin code
    
//}}} end code
}


void integralFunctionObject::end()
{
    if (false)
    {
        Info<<"end integral sha1: 3811ea68234a05845d6b0899299560d4ea54b3ad\n";
    }

//{{{ begin code
    
//}}} end code
}


void integralFunctionObject::timeSet()
{
    if (false)
    {
        Info<<"timeSet integral sha1: 3811ea68234a05845d6b0899299560d4ea54b3ad\n";
    }

//{{{ begin codeTime
    
//}}} end code
}


void integralFunctionObject::write()
{
    if (false)
    {
        Info<<"write integral sha1: 3811ea68234a05845d6b0899299560d4ea54b3ad\n";
    }

//{{{ begin code
    #line 98 "/home/luc/OpenFOAM/OpenQBMM-dev/tutorials/pbeTransportFoam/TaylorCouette/system/controlDict.functions.d43Average"
const volScalarField& m3 
                = mesh().lookupObject<volScalarField>("moment.3.populationBalance");
                
        const volScalarField& m4
                = mesh().lookupObject<volScalarField>("moment.4.populationBalance");
                
        volScalarField d43(m4/m3);

        scalar volAverage = 0;
        scalar totalVolume = 0;

        forAll (d43, cellI)
        {
            totalVolume += mesh().V()[cellI];
        }

        forAll (d43, cellI)
        {
            volAverage += d43[cellI]*mesh().V()[cellI]/totalVolume;
        }

        Info<<"Volume averaged normalized d43: " << volAverage/(2.0e-6) << endl;
//}}} end code
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

