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
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(aggBreakup, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label Foam::aggBreakup::nEqns() const
{
    return nBins_;
}


void Foam::aggBreakup::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    // This is the actual PBE
    scalarField aggregation(nBins_, 0.0);
    scalarField breakup(nBins_, 0.0);
    scalarField creation(nBins_, 0.0);
    scalarField destruction(nBins_, 0.0);

    // aggregation PBE by Smoluchowski (1917)
    if(isBrowninanAggOn_ || isShearAggOn_)
    {
        creation = 0.0;
        destruction = 0.0;

        if(isGridUniform_)
        {
            for(int i = 0; i < nBins_; ++i)
            {
                for(int j = 0; j < i; ++j)
                {
                    creation[i] += 0.5 * aggKernel_()[i-j-1][j] * y[i-j-1] * y[j];
                }
                for(int j = 0; j < nBins_; ++j)
                {
                    destruction[i] += aggKernel_()[i][j] * y[i] * y[j];
                }
            }
        }
        else
        {
            for(int i = 0; i < nBins_; ++i)
            {
                for(int j = 0; j < nBins_; ++j)
                {
                    for(int k = 0; k < nBins_; ++k)
                    {
                        creation[i] += 0.5 * chi_[j][k][i] * aggKernel_()[j][k] *
                                y[j] * y[k];
                    }
                }
                for(int j = 0; j < nBins_; ++j)
                {
                    destruction[i] += aggKernel_()[i][j] * y[i] * y[j];
                }
            }
        }

        aggregation = eta_ * ( creation - destruction);
    }

    // breakup PBE by Pandya & Spielman (1982)
    if(isBreakupOn_)
    {
        creation = 0.0;
        destruction = 0.0;

        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = i + 1; j < nBins_; ++j)
            {
                creation[i] += fragMassDistr_()[i][j] * breKernel_()[j] * y[j];
            }

            destruction[i] = breKernel_()[i] * y[i];
        }

        breakup = creation - destruction;
    }

    dydx = aggregation + breakup;
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

scalarField Foam::aggBreakup::smoluchowski(const scalarField& y)
{
    scalarField creation(nBins_, 0.0);
    scalarField destruction(nBins_, 0.0);

    if(isGridUniform_)
    {
        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < i; ++j)
            {
                creation[i] += 0.5 * aggKernel_()[i-j-1][j] * y[i-j-1] * y[j];
            }
            for(int j = 0; j < nBins_; ++j)
            {
                destruction[i] += aggKernel_()[i][j] * y[i] * y[j];
            }
        }
    }
    else
    {
        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                for(int k = 0; k < nBins_; ++k)
                {
                    creation[i] += 0.5 * chi_[j][k][i] * aggKernel_()[j][k] *
                                    y[j] * y[k];
                }
            }
            for(int j = 0; j < nBins_; ++j)
            {
                destruction[i] += aggKernel_()[i][j] * y[i] * y[j];
            }
        }
    }

    return eta_ * ( creation - destruction);

}

void Foam::aggBreakup::updateAggKernel(const label& celli)
{
    // clear aggKernel
    for(int i = 0; i < nBins_; ++i)
    {
        for(int j = 0; j < nBins_; ++j)
        {
            aggKernel_()[i][j] = 0.0;
        }
    }

    if(isShearAggOn_)
    {
       scalar G(2.0 * mag(strainRate_()[celli]));

       scalar alpha(4./3.);
       // NOTE: implement general alpha later

       for(int i = 0; i < nBins_; ++i)
       {
           for(int j = 0; j < nBins_; ++j)
           {
               aggKernel_()[i][j] += alpha * G * shearAggKernel_()[i][j];
           }
       }
    }
    if(isBrowninanAggOn_)
    {
        // NOTICE: review this!
        // I should rather define diffusivity elsewhere, so it would
        // be easier to add diffusivity generated by RBC collision.

        //Boltzmann constant [J K−1]
        scalar kB(1.3806488E-23);

        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                const scalar& mu = mu_[celli];
                if (mu <= 0.0)
                {
                    FatalError
                           << "Nonpositive dynamic viscosity in cell "
                           << celli << nl << exit(FatalError);
                }
                aggKernel_()[i][j] += 2./3. * kB * T_.value()
                        / mu_[celli] * brownianAggKernel_()[i][j];
            }
        }
    }
}

void Foam::aggBreakup::updateBreKernel(const label& celli)
{
    // clear breKernel
    for(int i = 0; i < nBins_; ++i)
    {
        breKernel_()[i] =  0.0;
    }

    if(isBreakupOn_)
    {
        for(int i = 1; i < nBins_; ++i) // begins at i=0 because monomer does not break
        {
            scalar G(2.0 * mag(strainRate_()[celli]));
            breKernel_()[i] = pow(G / Gstar_.value(), b_) * Rc_()[i];
        }
    }
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
    Rmono_(dimensionedScalar("Rmono", dimLength, 0.0)),
    T_(dimensionedScalar("T", dimTemperature, 300.0)),
    a_(dimensionedScalar("a", dimless/dimTime, 1.0)),
    Gstar_(dimensionedScalar("Gstar", dimless/dimTime, 1000.0))
{
    readDict("aggBreakupProperties");
    setCMD();
    setPBE();

    postProc.set
    (
        new aggBrePostprocess
        (
            aggBreakupDict_,
            CMD_,
            vList_,
            Rlist_

        )
    );

    if(debug)
    {
        writeData(Info);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::aggBreakup>
//Foam::aggBreakup::New()
//{
//    return autoPtr<aggBreakup>(new aggBreakup);
//}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aggBreakup::~aggBreakup()
{
    // De-Allocate memory to prevent memory leak

    // Calculate the interpolation operator on the non-uniform grid
    if(!isGridUniform_)
    {
        Info << "Deallocating the interpolation array" << nl << endl;

        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                delete [] chi_[i][j];
            }
            delete [] chi_[i];
        }
        delete [] chi_;
    }

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::aggBreakup::readDict(string dictName)
{
    Info << "Reading " << dictName << nl << endl;

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

    // Create the selected ODE system solver
    solver_ =  ODESolver::New(
                                *this,
                                aggBreakupDict_->subDict("odeSolver")
                              );
}


void Foam::aggBreakup::setCMD()
{
    Info << "Setting CMD up..." << nl << endl;

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
            Info << "Reading field " << "C_" + name << nl << endl;

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

            Info << "Writing field C_" + name + " = Cdefault" << nl << endl;
            CMD_[i].write();
        }
    }


}

void Foam::aggBreakup::setPBE()
{
    Info << "Setting the PBE up..." << nl << endl;

    // Calculate the interpolation operator on the non-uniform grid
    if(!isGridUniform_)
    {
        Info << "Allocating the interpolation (chi) array" << nl << endl;

        chi_ = new scalar**[nBins_];
        for(int i = 0; i < nBins_; ++i)
        {
            chi_[i] = new scalar*[nBins_];
            for(int j = 0; j < nBins_; ++j)
            {
                //allocates and initializes = 0
                chi_[i][j] = new scalar[nBins_]();
            }
        }

        scalar vSum[nBins_][nBins_];
        forAll(vList_, i)
        {
            forAll(vList_, j)
            {
                vSum[i][j] = vList_[i] + vList_[j];
            }
        }

        //NEEDS DEBUG!!!

        //no aggregation into the last bin
        for(int k = 0; k < nBins_ - 1; ++k)
        {
            for(int i = 0; i < nBins_; ++i)
            {
                for(int j = 0; j < nBins_; ++j)
                {
                    if
                    (
                        (vSum[i][j] <= vList_[k+1]) &&
                        (vSum[i][j] >= vList_[k])
                    )
                    {
                        chi_[i][j][k] = (vList_[k+1] - vSum[i][j])
                                        / (vList_[k+1] - vList_[k]);
                    }
                    if
                    (
                        (vSum[i][j] <= vList_[k]) &&
                        (vSum[i][j] >= vList_[k-1])
                    )
                    {
                        chi_[i][j][k] = (vSum[i][j] - vList_[k-1])
                                        / (vList_[k] - vList_[k-1]);
                    }
                }
            }

        }
    }

    // Generate the kernel matrices:
    if(isBrowninanAggOn_ || isShearAggOn_ )
    {
        Info << "Allocating the aggregation kernel array" << nl << endl;
        aggKernel_.set(new scalarSquareMatrix(nBins_));

        if(isBrowninanAggOn_)
        {
            // allocating brownian aggregation kernel
            brownianAggKernel_.set(new scalarSquareMatrix(nBins_));

            // filling the 2nd order array
            for(int i = 0; i < nBins_; ++i)
            {
                for(int j = 0; j < nBins_; ++j)
                {
                    brownianAggKernel_()[i][j] = pow(Rlist_[i].value()
                                                     + Rlist_[j].value(),
                                                     2.0)
                                                 / (Rlist_[i].value()
                                                    * Rlist_[j].value());
                }
            }

            if(debug)
            {
                //output file
                OFstream ofBrownKernel("kd");
                for(int i = 0; i < nBins_; ++i)
                {
                    for(int j = 0; j < nBins_; ++j)
                    {
                        ofBrownKernel << brownianAggKernel_()[i][j]
                                      << tab;
                    }
                    ofBrownKernel << nl;
                }
                ofBrownKernel << endl;
            }
        }

        if(isShearAggOn_)
        {
            // allocating array
            shearAggKernel_.set(new scalarSquareMatrix(nBins_));

            // filling array
            for(int i = 0; i < nBins_; ++i)
            {
                for(int j = 0; j < nBins_; ++j)
                {
                    shearAggKernel_()[i][j] = pow(Rlist_[i].value()
                                                  + Rlist_[j].value(),
                                                  3.0);
                }
            }

            if(debug)
            {
                //output file
                OFstream ofShearKernel("kf");
                for(int i = 0; i < nBins_; ++i)
                {
                    for(int j = 0; j < nBins_; ++j)
                    {
                        ofShearKernel << shearAggKernel_()[i][j]
                                      << tab;
                    }
                    ofShearKernel << nl;
                }
            }

            volTensorField gradU = fvc::grad(U_);

            strainRate_.set
            (
                new volTensorField
                (
                    IOobject
                    (
                        "E",
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    0.5 * (gradU + gradU.T())
                )
            );

            Info << "Writing field E" << nl << endl;
            strainRate_->write();

            symmGradU_.set
            (
                new volSymmTensorField
                (
                    IOobject
                    (
                        "symmGradU",
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    symm(gradU)
                )
            );

            Info << "Writing field symmGradU" << nl << endl;
            symmGradU_->write();
        }
    }
    else
    {
        WarningIn("Foam::aggBreakup::setPBE()")
            << "Aggregation OFF!"
            << nl << endl;
    }

    if(isBreakupOn_)
    {
        breKernel_.set(new scalarField(nBins_, 0.0));
        Rc_.set(new scalarField(nBins_, 0.0));

        // fill the radii part of the breakup kernel. The shear part depends
        // on the cell shear rate, therefore the breKernel is calculated
        // elsewhere.
        for(int i = 0; i < nBins_; ++i)
        {
            Rc_()[i] = pow(Rlist_[i].value() / Rlist_[0].value(), c_);
        }

        if(debug)
        {
            Info << "Writing (R_i/R_p)^c in file \"kb\"" << nl << endl;
            //output file
            OFstream ofBreakupKernelDebug("kb");
            for(int i = 0; i < nBins_; ++i)
            {
                ofBreakupKernelDebug << Rc_()[i]
                                     << tab;
            }
            ofBreakupKernelDebug << endl;
        }

        if(a_.value() != 1.0)
        {
            Info << "In the dictionary " << aggBreakupDict_->name()
                 << " the breakup constant \'a != 1 s^-1\' was given."
                 << " Calculating new Gstar from a, b, c, and R_1." << endl;
            Gstar_.value() = pow(a_.value() * pow(Rlist_[0].value(), c_),
                                 - 1.0  / b_);
        }

        // Set binary fragment mass distribution
        // You can implement other fragment mass distributions if you like
        // (e.g., normal, erosion, uniform...)
        fragMassDistr_.set(new scalarSquareMatrix(nBins_));

        if(isGridUniform_)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                scalarField p(nBins_, 0.0);
                for(int i = 0; i < (j + 1) / 2; ++i)
                {
                    p[j/2] = 1.0;
                }
                for(int i = 0; i < nBins_; ++i)
                {
                    fragMassDistr_()[i][j] = p[i];
                    if(j-i-1 >= 0)
                    {
                        fragMassDistr_()[i][j] += p[j-i-1];
                    }
                }
            }
        }
        else
        {
            for(int j = 0; j < nBins_; ++j)
            {
                for(int i = 0; i < nBins_; ++i)
                {
                    if(i == j - 1)
                    {
                        fragMassDistr_()[i][j] = 2.0;
                    }
                    else
                    {
                        fragMassDistr_()[i][j] = 0.0;
                    }
                }
            }
        }

        if(debug)
        {
            // write the fragment mass distribution matrix in a file named "g"
            OFstream ofFragMassDistr("g");
            for(int i = 0; i < nBins_; i++)
            {
                for(int j = 0; j < nBins_; j++)
                {
                    ofFragMassDistr << fragMassDistr_()[i][j]
                                    << tab;
                }
                ofFragMassDistr << nl;
            }
            ofFragMassDistr << endl;
        }
    }
    else
    {
        WarningIn("Foam::aggBreakup::setPBE()")
            << "Breakup OFF!"
            << nl << endl;
    }

}


void Foam::aggBreakup::solveAggBreakup()
{
    scalar eps = readScalar
                (
                    aggBreakupDict_->subDict("odeSolver").lookup("eps")
                );

    Info << "Solving PBE" << tab
         << "T = " << runTime_.value() << " s"
         << endl;

    forAll(mesh_.C(), celli)
    {
        scalarField y(nBins_, 0.0);

        // set values to vector y
        forAll(y, dofi)
        {
            y[dofi] = CMD_[dofi][celli];
        }

        if(isBrowninanAggOn_ || isShearAggOn_)
        {
            updateAggKernel(celli);
        }

        if(isBreakupOn_)
        {
            updateBreKernel(celli);
        }
        solver_->solve
        (
            runTime_.value(),
            runTime_.value() + runTime_.deltaT().value(),
            y,
            eps
        );

        // update CMD in celli
        forAll(y, dofi)
        {
            CMD_[dofi][celli] = y[dofi];
        }
    }

    postProc->update();
}


void Foam::aggBreakup::solveTransport()
{
    forAll(CMD_, i)
    {
        volScalarField& Ci = CMD_[i];

        tmp<fvScalarMatrix> massTransport
        (
            fvm::ddt(Ci)
          + fvm::div(phi_, Ci, "div(phi,C_*)")
//          - fvm::laplacian(D, Ci, "laplacian(D,C_*)")
        );

       massTransport().solve(mesh_.solver("C_*"));
    }

    if(isShearAggOn_)
    {
        // update strainRate field
        volTensorField gradU = fvc::grad(U_);
        strainRate_() = 0.5*(gradU + gradU.T());

        symmGradU_() = symm(gradU);

        // NOTICE:
        // strainRate_() = symm(gradU); could also be used.
        // the problem is that the symm operator generates a symmTensor instead
        // of a simple tensor. The symmTensor stores 6 scalars instead of 9 of
        // tensor. The strainRate field would have to be defined as symmTensor.
    }
}


bool Foam::aggBreakup::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const aggBreakup& sds)
{
    os.writeKeyword("isGridUniform") << sds.isGridUniform_
                                     << token::END_STATEMENT
                                     << endl;
    os.writeKeyword("isBrowninanAggOn") << sds.isBrowninanAggOn_
                                        << token::END_STATEMENT
                                        << endl;
    os.writeKeyword("isShearAggOn") << sds.isShearAggOn_
                                    << token::END_STATEMENT
                                    << endl;
    os.writeKeyword("isBreakupOn") << sds.isBreakupOn_
                                   << token::END_STATEMENT
                                   << endl;

    os.writeKeyword("nBins") << sds.nBins_
                             << token::END_STATEMENT
                             << endl;
    os.writeKeyword("Rmono") << sds.Rmono_
                             << token::END_STATEMENT
                             << endl;
    os.writeKeyword("DF") << sds.DF_
                          << token::END_STATEMENT
                          << endl;
    os.writeKeyword("eta") << sds.eta_
                           << token::END_STATEMENT
                           << endl;
    os.writeKeyword("T") << sds.T_
                         << token::END_STATEMENT
                         << endl;

    os.writeKeyword("a") << sds.a_
                         << token::END_STATEMENT
                         << endl;
    os.writeKeyword("Gstar") << sds.Gstar_
                             << token::END_STATEMENT
                             << endl;
    os.writeKeyword("b") << sds.b_
                         << token::END_STATEMENT
                         << endl;
    os.writeKeyword("c") << sds.c_
                         << token::END_STATEMENT
                         << endl;

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
