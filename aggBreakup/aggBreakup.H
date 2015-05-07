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
    Foam::aggBreakup

Description

SourceFiles
    aggBreakupI.H
    aggBreakup.C
    aggBreakupIO.C

\*---------------------------------------------------------------------------*/

#ifndef aggBreakup_H
#define aggBreakup_H

#include "volFields.H"
#include "PtrList.H"
#include "fvCFD.H"
#include "hashedWordList.H"
#include "ODESystem.H"
#include "ODESolver.H"

#include "../aggBrePostprocess/aggBrePostprocess.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <typename T>
std::string to_string(T value)
{
    std::ostringstream os ;
    os << value ;
    return os.str() ;
}

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class aggBreakup Declaration
\*---------------------------------------------------------------------------*/

class aggBreakup
:
    private ODESystem
{
    // Private data

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        //-Viscosity
        const volScalarField& mu_;

        const fvMesh& mesh_;
        const Time& runTime_;

        //-Colloid cluster fields
        PtrList<volScalarField> CMD_;

        PtrList<volScalarField> CMDold_;

        //-Absolute velocity gradient
        autoPtr<volScalarField> GA_;

        // CMD discretization scheme data

            bool isPBEcoupled_;

            autoPtr<IOdictionary> aggBreakupDict_;

            //-Specify the bin structure you want to use
            // * True => use a linear bin structure and solve the PBE in
            // standard form.
            // * False => use a non-linear (2nd order geometrically spaced) bin
            // structure and solve the PBE using a interpolation operator.
            bool isGridUniform_;

            //-Number of discretization points of the CMD
            label nBins_;

            //-Number of monomers in the cluster 'i'
            PtrList<scalar> vList_;

        // Cluster structure data

            //-Monomer mean radius (m)
            dimensionedScalar Rmono_;

            //-Radius of gyration the cluster 'i'
            PtrList<dimensionedScalar> Rlist_;

            //-Monomer volume (m^3)
            dimensionedScalar Vmono_;

            //-Cluster fractal dimension
            scalar DF_;

        // PBE data

            //-Array of terms obtained from interpolation function of
            // Kumar & Ramkrishna (1997), which also presented by
            // Garrick, Lehtinen & Zachariah (2006)
            scalar*** chi_;

            // Aggregation

                bool isBrowninanAggOn_;
                bool isShearAggOn_;

                //-Aggregation efficiency
                scalar eta_;

                //-Shear collision kernel shape factor
                tmp<volScalarField> alpha_;

                //-Base Temperature
                dimensionedScalar T_;

                //-Base aggregation kernel (shear excluded)
                autoPtr<scalarSquareMatrix> brownianAggKernel_;
                autoPtr<scalarSquareMatrix> shearAggKernel_;
                autoPtr<scalarSquareMatrix> aggKernel_;

            // Breakup data

                bool isBreakupOn_;

                dimensionedScalar a_;
                dimensionedScalar Gstar_;
                scalar b_;
                scalar c_;

                //-Base breakup kernel (just the Radii part, shear part
                // excluded)
                autoPtr<scalarField> Rc_;

                //-Breakup kernel
                autoPtr<scalarField> breKernel_;

                //-Fragment mass distribution
                autoPtr<scalarSquareMatrix> fragMassDistr_;

            // Activation rule

                bool isActivationOn_;
                scalar activThreshold_; //shear rate
                scalar activRate_;
                autoPtr<volScalarField> Crp_;

        // ODE solver
        autoPtr<ODESolver> solver_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        aggBreakup(const aggBreakup&);

        //- Disallow default bitwise assignment
        void operator=(const aggBreakup&);

        // ODE parameters

            //- Return number of equations
            virtual label nEqns() const;
//            {
//                return label(nBins_);
//            }

            //- Return derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Return Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;

        scalarField smoluchowski(const scalarField& y);

        void updateAggKernel(const label& celli);
        void updateBreKernel(const label& celli);

        tmp<volScalarField> S(label clusterI) const;

        tmp<volScalarField> kag(label i, label j) const;
        tmp<volScalarField> kbr(label i) const;

public:

    ClassName("aggBreakup");

    // Public data
    autoPtr<aggBrePostprocess> postProc;


    // Constructors

        //- Construct from components
        aggBreakup
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& mu
        );


//        //- Construct as copy
//        aggBreakup(const aggBreakup&);


    // Selectors

        //- Select null constructed
        static autoPtr<aggBreakup> New();


    //- Destructor
    ~aggBreakup();


    // Member Functions

        // Access

        void readDict(string dictName);

        void setCMD();

        void setPBE();

        // Check

        // Edit

            void update();

            void solveAggBreakup();

            void solveTransport();

        // Write

            //- WriteData member function required by regIOobject
            bool writeData(Ostream&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const aggBreakup&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "aggBreakupI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
