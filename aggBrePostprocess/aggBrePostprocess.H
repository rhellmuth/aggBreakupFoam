/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

Class
    Foam::aggBrePostprocess

Description

SourceFiles
    aggBrePostprocessI.H
    aggBrePostprocess.C
    aggBrePostprocessIO.C

\*---------------------------------------------------------------------------*/

#ifndef aggBrePostprocess_H
#define aggBrePostprocess_H

#include "volFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//// Forward declaration of classes
//class Istream;
//class Ostream;

//// Forward declaration of friend functions and operators
//class aggBrePostprocess;
//Istream& operator>>(Istream&, aggBrePostprocess&);
//Ostream& operator<<(Ostream&, const aggBrePostprocess&);


/*---------------------------------------------------------------------------*\
                         Class aggBrePostprocess Declaration
\*---------------------------------------------------------------------------*/

class aggBrePostprocess
{
    // Private data

        //- Dictionary
        dictionary aggBreakupDict_;

        //-Colloid cluster fields
        const PtrList<volScalarField>& CMD_;

        const Time& runTime_;
        const fvMesh& mesh_;

        const PtrList<scalar>& vList_;
        const PtrList<dimensionedScalar>& rList_;

        //- Swiches
        bool isThetaFieldOn_;
        bool isAggCharTimeOn_;
        bool isBreCharTimeOn_;
        bool isAggPecletOn_;
        bool isAvgConcentrationOn_;
        bool isZerothMomentOn_;
        bool isFirstMomentOn_;
        bool isSecondMomentOn_;
        bool isVmeanOn_;
        bool isRmeanOn_;
        bool isAggFractionOn_;
        bool isCellFlowCharTimeOn_;
        bool isAggAdvRatioOn_;

        //-Aggregation characteristic time field
        autoPtr<volScalarField> ta_;

        //-Breakup characteristic time field
        autoPtr<volScalarField> tb_;

        //-Aggregation-breakup ratio field
        autoPtr<volScalarField> theta_;

        //-Aggregation Péclet number field
        autoPtr<volScalarField> aggPe_;

        //-Average monomer concentration in the whole system
        //scalar ta_;

        //-Zeroth moment field
        autoPtr<volScalarField> M0_;

        //-First moment field
        autoPtr<volScalarField> M1_;

        //-Second moment field
        autoPtr<volScalarField> M2_;

        //-Mean number of monomers per aggregate
        autoPtr<volScalarField> vMean_;

        //-Aggregates mean radius of gyration
        autoPtr<volScalarField> rMean_;

        //-Aggregation fraction in %
        autoPtr<volScalarField> PA_;

        //-Convection characteristic time
        autoPtr<volScalarField> tf_;

        //-Damköhler number for aggregation rate / advection rate
        autoPtr<volScalarField> Da_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        aggBrePostprocess(const aggBrePostprocess&);

        //- Disallow default bitwise assignment
        void operator=(const aggBrePostprocess&);


public:

    ClassName("aggBrePostprocess");

    // Constructors

        //- Construct from dictionary
        aggBrePostprocess
        (
                const dictionary& aggBreakupDict,
                const PtrList<volScalarField>& CMD,
                const PtrList<scalar>& vList,
                const PtrList<dimensionedScalar>& rList
        );

    // Selectors

        //- Select null constructed
        static autoPtr<aggBrePostprocess> New();


    //- Destructor
    ~aggBrePostprocess();


    // Member Functions

        // Access

        // Check

        // Edit
        tmp<volScalarField> G() const;
        tmp<volScalarField> M_0() const;
        tmp<volScalarField> M_1() const;
        tmp<volScalarField> M_2() const;
        tmp<volScalarField> vMean() const;
        tmp<volScalarField> rMean() const;
        tmp<volScalarField> t_a() const;
        tmp<volScalarField> t_b() const;
        tmp<volScalarField> theta() const;
        tmp<volScalarField> PA() const;
        tmp<volScalarField> t_f() const;
        tmp<volScalarField> Da() const;
//        tmp<volScalarField> t_d() const;
//        tmp<volScalarField> Pe() const;

        void update();

        // Write



    // IOstream Operators

//        friend Istream& operator>>(Istream&, aggBrePostprocess&);
        friend Ostream& operator<<(Ostream&, const aggBrePostprocess&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "aggBrePostprocessI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
