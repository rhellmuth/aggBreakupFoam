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

\*---------------------------------------------------------------------------*/

#include "aggBreakup.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//const Foam::label Foam::aggBreakup::nEqns() const
//{
//    return label(nBins_);
//}

void Foam::aggBreakup::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    // Set the derivatives for displacement
    dydx[0] = 1;
    dydx[1] = 2;
    dydx[2] = 3;
}


void Foam::aggBreakup::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    notImplemented("Oscilator::jacobian(...) const");
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::aggBreakup::aggBreakup
(
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& mu
)
:
    ODESystem(),
    U_(U),
    phi_(phi),
    mu_(mu),
    mesh_(U_.mesh()),
    runTime_(U_.time()),
    Rmono_(dimensionedScalar("Rmono", dimensionSet(0,1,0,0,0,0,0), 0.0)),
    T_(dimensionedScalar("T", dimTemperature, 300.0)),
    a_(dimensionedScalar("a", dimensionSet(0,0,-1,0,0,0,0), 1.0)),
    Gstar_(dimensionedScalar("Gstar", dimensionSet(0,0,-1,0,0,0,0), 1000.0))
{
    readDict("aggBreakupProperties");
    setCMD();
    setPBE();

    writeData(Info);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::aggBreakup>
//Foam::aggBreakup::New()
//{
//    return autoPtr<aggBreakup>(new aggBreakup);
//}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aggBreakup::~aggBreakup()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::aggBreakup::readDict(string dictName)
{
    Info << "Reading " << dictName << endl << endl;

    aggBreakupDict_.set
    (
        new IOdictionary
        (
            IOobject
            (
                dictName,
                runTime_.constant(),
                runTime_.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    isGridUniform_ = readBool(aggBreakupDict_->lookup("isGridUniform"));
    nBins_ = readLabel( aggBreakupDict_->lookup("nBins") );
    Rmono_ = aggBreakupDict_->lookup("Rmono");
    DF_ = readScalar(aggBreakupDict_->lookup("DF"));

    isBrowninanAggOn_ = readBool(aggBreakupDict_->lookup("isBrowninanAggOn"));
    isShearAggOn_ = readBool(aggBreakupDict_->lookup("isShearAggOn"));
    eta_ = readScalar( aggBreakupDict_->lookup("eta") );
    T_ = aggBreakupDict_->lookup("T");

    isBreakupOn_ = readBool(aggBreakupDict_->lookup("isBreakupOn"));
//    a_ = aggBreakupDict_->lookup("a");
    Gstar_ = aggBreakupDict_->lookup("Gstar");
    b_ = readScalar(aggBreakupDict_->lookup("b"));
    c_ = readScalar(aggBreakupDict_->lookup("c"));
}


void Foam::aggBreakup::setCMD()
{
    Info << "Setting CMD up..." << endl << endl;

    vList_.setSize(nBins_);
    Rlist_.setSize(nBins_);
    CMD_.setSize(nBins_);

    forAll(vList_, i)
    {
        if (isGridUniform_)
        {
            // Linear bin grid
            vList_.set(i, new scalar(i + 1));

            // Cluster radius of gyration in bin i
            Rlist_.set
            (
                i,
                new dimensionedScalar(pow(vList_[i], 1.0/DF_) * Rmono_)
            );
        }
        else
        {
            // Geometric bin grid
            vList_.set
            (
                i,
                new scalar(pow(2.0, i))
            );

            // Cluster radius of gyration in bin i
            Rlist_.set
            (
                i,
                new dimensionedScalar(pow(vList_[i], 1.0/DF_) * Rmono_)
            );
        }

        // Rename radii of gyration acording to number of monomers
        // in the cluster
        Rlist_[i].name() = "R_" + to_string(vList_[i]);
    }

    forAll(vList_, i)
    {
        const word name = to_string(vList_[i]);

        IOobject header
        (
            "C_" + name,
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        // check if dictionary exists and can be read
        if (header.headerOk())
        {
            Info << "Reading field C_" + name << endl << endl;

            CMD_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "C_" + name,
                        runTime_.timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            );
        }
        else
        {
            Info << "Warning: Field C_" + name + " not found!"
                 << " Reading Cdefault instead" << endl;

            volScalarField Cdefault
            (
                IOobject
                (
                    "Cdefault",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh_
            );

            CMD_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "C_" + name,
                        runTime_.timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    Cdefault
                )
            );

            Info << "Writing field C_" + name + " = Cdefault" << endl << endl;
            CMD_[i].write();
        }
    }
}

void Foam::aggBreakup::setPBE()
{
    Info << "Setting the PBE up..." << endl << endl;

    // Calculate the interpolation operator on the non-uniform grid
//            chi = interpolate(v_list);
    // Generate the kernel matrices:
    if((isBrowninanAggOn_ == true) && (isShearAggOn_ == true))
    {
//        aggKernel = Brownian_kernel(R_list) + shear_kernel(R_list);
    }
    else if((isBrowninanAggOn_ == true) && (isShearAggOn_ == false))
    {
//        aggKernel = Brownian_kernel(R_list);
    }
    else if((isBrowninanAggOn_ == false) && (isShearAggOn_ == true))
    {
//        aggKernel = shear_kernel(R_list);
    }
    else
    {
//        print "Error: Aggregation OFF!";
    }

    if(isBreakupOn_)
    {
//        breKernel = breakup_kernel(R_list);

        if(isGridUniform_)
        {
//            fragDistr = binarySplitFrag();
        }
        else
        {
//            fragDistr = nonunifBinaryFrag();
        }
    }
    else
    {
//        fragDistr = np.zeros((n_bin, n_bin));
    }

}

void Foam::aggBreakup::solveAggBreakup()
{}


bool Foam::aggBreakup::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const aggBreakup& sds)
{
    os.writeKeyword("isGridUniform") << sds.isGridUniform_
        << token::END_STATEMENT << endl;
    os.writeKeyword("isBrowninanAggOn") << sds.isBrowninanAggOn_
        << token::END_STATEMENT << endl;
    os.writeKeyword("isShearAggOn") << sds.isShearAggOn_
        << token::END_STATEMENT << endl;
    os.writeKeyword("isBreakupOn") << sds.isBreakupOn_
        << token::END_STATEMENT << endl;

    os.writeKeyword("nBins") << sds.nBins_
        << token::END_STATEMENT << endl;
    os.writeKeyword("Rmono") << sds.Rmono_
        << token::END_STATEMENT << endl;
    os.writeKeyword("DF") << sds.DF_
        << token::END_STATEMENT << endl;
    os.writeKeyword("eta") << sds.eta_
        << token::END_STATEMENT << endl;
    os.writeKeyword("T") << sds.T_
        << token::END_STATEMENT << endl;

    os.writeKeyword("a") << sds.a_
        << token::END_STATEMENT << endl;
    os.writeKeyword("Gstar") << sds.Gstar_
        << token::END_STATEMENT << endl;
    os.writeKeyword("b") << sds.b_
        << token::END_STATEMENT << endl;
    os.writeKeyword("c") << sds.c_
        << token::END_STATEMENT << endl;

    os << endl << "CMD:" << endl;
    os << "v" << tab
       << "R" << endl;
    forAll(sds.vList_,i)
    {
        os << sds.vList_[i] << tab
           << sds.Rlist_[i] << endl;
    }

    return os;
}

// ************************************************************************* //
