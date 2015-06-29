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
    return 12;
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
    if(isBrowninanAggOn_ || isShearAggOn_ || isSorensenianAggOn_)
    {
        creation = 0.0;
        destruction = 0.0;

        if(isGridUniform_)
        {
            for(int i = 0; i < nBins_; ++i)
            {
                for(int j = 0; j < i; ++j)
                {
                    creation[i] += 0.5 * kc_()[i-j-1][j]
                                    * y[i-j-1] * y[j];
                }
                for(int j = 0; j < nBins_; ++j)
                {
                    destruction[i] += kc_()[i][j] * y[i] * y[j];
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
                        creation[i] += 0.5 * chi_[j][k][i] * kc_()[j][k]
                                        * y[j] * y[k];
                    }

                    destruction[i] += kc_()[i][j] * y[i] * y[j];
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
                creation[i] += fragMassDistr_()[i][j] * kb_()[j] * y[j];
            }

            destruction[i] = kb_()[i] * y[i];
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
    notImplemented("aggBreakup::jacobian(...) const");
}


void Foam::aggBreakup::updateKernelsInCell(const label& celli)
{
    const scalar& G = GA_()[celli];

    if(isShearAggOn_ || isBrowninanAggOn_ || isSorensenianAggOn_)
    {
        // clear aggKernel
        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                kc_()[i][j] = 0.0;

                if(isShearAggOn_)
                {
                    kc_()[i][j] += G * ksBase_()[i][j];
                }
                if(isBrowninanAggOn_ || isSorensenianAggOn_)
                {
                    kc_()[i][j] += kdBase_()[i][j];
                    if(isSorensenianAggOn_)
                    {
                        kc_()[i][j] += 1.05 * G * kdBase_()[i][j];
                    }
                }
            }
        }
    }

    if(isBreakupOn_)
    {
        for(int i = 0; i < nBins_; ++i)
        {
            kb_()[i] = pow(G / Gstar_.value(), b_) * kbBase_()[i];
        }
    }
}

Foam::dimensionedScalar Foam::aggBreakup::einsteinStokes
(
    const dimensionedScalar& R,
    const dimensionedScalar& T,
    const dimensionedScalar& mu
)
{
    //Defines brownian diffusivity constant using to Einstein-Stokes equation.

    //Boltzmann constant [J K−1]
    dimensionedScalar kB
    (
        "k_B",
        dimEnergy / dimTemperature,
        1.3806488E-23
    );
    return kB * T / (6.0 * M_PI * mu * R);
}

Foam::dimensionedScalar Foam::aggBreakup::radiusOfGyration
(
    const scalar& v,
    const scalar& D_F,
    const dimensionedScalar& R_p
)
{
    //Radius of gyration of a fractal aggregate ball.
    return pow(v, 1.0/D_F) * R_p;
}

void Foam::aggBreakup::allocChi()
{
    //Interpolation function of Kumar & Ramkrishna (1997),
    //also presented by Garrick, Lehtinen & Zachariah (2006).

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

    for(int k = 0; k < nBins_; ++k)
    {
        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                if(k < nBins_ - 1) //all but the last bin
                {
                    if
                    (
                           (vSum[i][j] <= vList_[k+1])
                        && (vSum[i][j] >= vList_[k])
                    )
                    {
                        chi_[i][j][k] =   (vList_[k+1] - vSum[i][j])
                                        / (vList_[k+1] - vList_[k]);
                    }
                }
                if
                (
                       (vSum[i][j] <= vList_[k])
                    && (vSum[i][j] >= vList_[k-1])
                )
                {
                    chi_[i][j][k] =   (vSum[i][j] - vList_[k-1])
                                    / (vList_[k] - vList_[k-1]);
                }
            }
        }
        if(debug)
        {
            //output file
            OFstream ofChi("chi_k" + to_string(k));
            for(int i = 0; i < nBins_; ++i)
            {
                for(int j = 0; j < nBins_; ++j)
                {
                    ofChi << chi_[i][j][k] << tab;
                }
                ofChi << nl;
            }
        }
    }
}


Foam::autoPtr<scalarSquareMatrix> Foam::aggBreakup::Brownian_kernel()
{
    autoPtr<scalarSquareMatrix> k_d(new scalarSquareMatrix(nBins_));
    for(int i = 0; i < nBins_; ++i)
    {
        for(int j = 0; j < nBins_; ++j)
        {
            //Remove i,j = n from Smoluchowski collision kernel. Aggregation
            // with the largest class of aggregate would cause mass loss.
            if((i == nBins_-1) || (j == nBins_-1))
            {
                k_d()[i][j] = 0.0;
            }
            else
            {
                k_d()[i][j] =
                (
                4.0 * M_PI
                * (Rlist_[i].value() + Rlist_[j].value())
                * (Dlist_[i].value() + Dlist_[j].value())
                );
            }
        }
    }

    if(debug)
    {
        //output file
        OFstream ofBrownKernel("k_d");
        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                ofBrownKernel << k_d()[i][j]
                              << tab;
            }
            ofBrownKernel << nl;
        }
        ofBrownKernel << endl;
    }
    return k_d;
}


Foam::autoPtr<scalarSquareMatrix> Foam::aggBreakup::shear_kernelBase()
{
    // allocating shear collision kernel (excluding shear)
    autoPtr<scalarSquareMatrix> ksBase(new scalarSquareMatrix(nBins_));

    for(int i = 0; i < nBins_; ++i)
    {
        for(int j = 0; j < nBins_; ++j)
        {
            //Remove i,j = n from Smoluchowski collision kernel. Aggregation
            // with the largest class of aggregate would cause mass loss.
            if((i == nBins_-1) || (j == nBins_-1))
            {
                ksBase()[i][j] = 0.0;
            }
            else
            {
                ksBase()[i][j] =
                (
                    (32.0/3.0) * pow
                    (
                        Rlist_[i].value() + Rlist_[j].value(),
                        3.0
                    )
                );
            }
        }
    }

    if(debug)
    {
        //output file
        OFstream ofShearKernel("k_s");
        for(int i = 0; i < nBins_; ++i)
        {
            for(int j = 0; j < nBins_; ++j)
            {
                ofShearKernel << ksBase()[i][j]
                              << tab;
            }
            ofShearKernel << nl;
        }
    }

    return ksBase;
}

Foam::autoPtr<scalarField> Foam::aggBreakup::breakup_kernelBase()
{
    autoPtr<scalarField> kbBase(new scalarField(nBins_, 0.0));
    // fill the radii part of the breakup kernel. The shear part depends
    // on the cell shear rate, therefore the breKernel is calculated
    // elsewhere.
    for(int i = 0; i < nBins_; ++i)
    {
        //monomer doesn't break
        if(i==0)
        {
             kbBase()[i] = 0.0;
        }
        else
        {
            kbBase()[i] = pow(Rlist_[i].value() / Rlist_[0].value(), c_);
        }
    }

    if(debug)
    {
        Info << "Writing (R_i/R_p)^c in file \"k_b\"" << nl << endl;
        //output file
        OFstream ofBreakupKernelDebug("k_b");
        for(int i = 0; i < nBins_; ++i)
        {
            ofBreakupKernelDebug << kbBase()[i]
                                 << tab;
        }
        ofBreakupKernelDebug << endl;
    }

    if(a_.value() != 1.0)
    {
        Info << "In the dictionary " << aggBreakupDict_->name()
             << " the breakup constant \'a != 1 s^-1\' was given."
             << " Calculating new Gstar from a, b, c, and R_p." << endl;
        Gstar_.value() = pow
        (
            a_.value() * pow(Rlist_[0].value(), c_),
            -1.0/b_
        );
    }

    return kbBase;
}

Foam::autoPtr<scalarSquareMatrix> Foam::aggBreakup::fragMassDistr()
{
    // allocating fragment mass distribution
    autoPtr<scalarSquareMatrix> g(new scalarSquareMatrix(nBins_));

    // Set binary fragment mass distribution for uniform and geometric grids
    // You can implement other fragment mass distributions if you like
    // (e.g., normal, erosion, uniform...)
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
                g()[i][j] = p[i];
                if(j-i-1 >= 0)
                {
                    g()[i][j] += p[j-i-1];
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
                    g()[i][j] = 2.0;
                }
                else
                {
                    g()[i][j] = 0.0;
                }
            }
        }
    }

    if(debug)
    {
        // write the fragment mass distribution matrix in a file
        // named "fragMassDistr"
        OFstream ofFragMassDistr("fragMassDistr");
        for(int i = 0; i < nBins_; i++)
        {
            for(int j = 0; j < nBins_; j++)
            {
                ofFragMassDistr << g()[i][j]
                                << tab;
            }
            ofFragMassDistr << nl;
        }
        ofFragMassDistr << endl;
    }

    return g;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::aggBreakup::aggBreakup
(
        const volVectorField& U,
        const surfaceScalarField& phi
)
:
    ODESystem(),
    U_(U),
    phi_(phi),
    mesh_(U_.mesh()),
    runTime_(U_.time()),
    Rp_(dimensionedScalar("Rmono", dimLength, 0.0)),
    T_(dimensionedScalar("T", dimTemperature, 300.0)),
    muPlasma_(dimensionedScalar("mu", dimMass/dimLength/dimTime, 3.5e-03)),
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
    // Deallocate memory to prevent memory leak

    // Deallocate the interpolation operator on the non-uniform grid
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

    isPBEcoupled_ = readBool(aggBreakupDict_->lookup("isPBEcoupled"));

    if(!isPBEcoupled_)
    {
        // Create the selected ODE system solver
        solver_ =  ODESolver::New
        (
            *this,
            aggBreakupDict_->subDict("odeSolver")
        );
    }

    isGridUniform_ = readBool(aggBreakupDict_->lookup("isGridUniform"));
    nBins_ = readLabel( aggBreakupDict_->lookup("nBins") );
    Rp_ = aggBreakupDict_->lookup("Rmono");
    DF_ = readScalar(aggBreakupDict_->lookup("DF"));

    isBrowninanAggOn_ = readBool(aggBreakupDict_->lookup("isBrowninanAggOn"));
    isShearAggOn_ = readBool(aggBreakupDict_->lookup("isShearAggOn"));
    isSorensenianAggOn_ = readBool(aggBreakupDict_->lookup("isSorensenianAggOn"));
    eta_ = readScalar( aggBreakupDict_->lookup("eta") );
    T_ = aggBreakupDict_->lookup("T");
    muPlasma_ = aggBreakupDict_->lookup("muPlasma");

    isBreakupOn_ = readBool(aggBreakupDict_->lookup("isBreakupOn"));
//    a_ = aggBreakupDict_->lookup("a");
    Gstar_ = aggBreakupDict_->lookup("Gstar");
    b_ = readScalar(aggBreakupDict_->lookup("b"));
    c_ = readScalar(aggBreakupDict_->lookup("c"));

    isActivationOn_ = readBool(aggBreakupDict_->lookup("isActivationOn"));
    activThreshold_ = readScalar(aggBreakupDict_->lookup("activThreshold"));

}


void Foam::aggBreakup::setCMD()
{
    Info << "Setting CMD up..." << nl << endl;

    vList_.setSize(nBins_);
    Rlist_.setSize(nBins_);
    CMD_.setSize(nBins_);
    Dlist_.setSize(nBins_);
    if(isPBEcoupled_)
    {
        CMDold_.setSize(nBins_);
    }

    forAll(vList_, i)
    {
        // Number of monomers per cluster of class i
        if (isGridUniform_)
        {
            // Linear bin grid
            vList_.set(i, new scalar(i + 1));
        }
        else
        {
            // Geometric bin grid
            vList_.set(i,  new scalar(pow(2.0, i)));
        }

        const word name = to_string(vList_[i]);

        // Cluster radius of gyration in class of cluster i
        Rlist_.set
        (
            i,
            new dimensionedScalar
            (
                "R_" + to_string(vList_[i]),
                radiusOfGyration(vList_[i], DF_, Rp_)
            )
        );

        // Brownian diffusivity in class of cluster i
        Dlist_.set
        (
            i,
            new dimensionedScalar
            (
                "D_" + to_string(vList_[i]),
                einsteinStokes(Rlist_[i], T_, muPlasma_)
            )
        );

        IOobject header
        (
            "C_" + name,
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        // Check if dictionary exists and can be read
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

        if(isPBEcoupled_)
        {
            CMDold_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Cold_" + name,
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    CMD_[i]
                )
            );
        }
    }

    if(isActivationOn_)
    {
        Info << "Reading field " << "C_RP"<< nl << endl;

        Crp_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "C_RP",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        vWF_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "vWF",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("vWF", dimless/dimTime, 0.0)
            )
        );
    }
}

void Foam::aggBreakup::setPBE()
{
    Info << "Setting the PBE up..." << nl << endl;

    // According to Pedocchi & Piedra-Cueva (2005) the aggregation kernel in general flow
    // depends on the absolute velocity gradient, as had been previously formulated by
    // Camp & Stein (1943).
    GA_.set
    (
        new volScalarField
        (
            IOobject
            (
                "G_A",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt(2.0) * mag(symm(fvc::grad(U_)()))
        )
    );

    Info << "Writing field G_A" << nl << endl;
    GA_->write();

    // Calculate the interpolation operator on the non-uniform grid
    if(!isGridUniform_)
    {
        allocChi();
    }

    // Generate the kernel matrices:
    if(isBrowninanAggOn_ || isShearAggOn_ )
    {
        Info << "Allocating the collision kernel array" << nl << endl;
        kc_.set(new scalarSquareMatrix(nBins_));

        if(isBrowninanAggOn_ || isSorensenianAggOn_)
        {
            // allocating brownian collision kernel
            kdBase_ = Brownian_kernel();
        }
        if(isShearAggOn_)
        {
            // allocating shear collision kernel (excluding shear)
            ksBase_ = shear_kernelBase();
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
        Info << "Allocating the breakup kernel array" << nl << endl;
        kb_.set(new scalarField(nBins_, 0.0));

        kbBase_ = breakup_kernelBase();

        fragMassDistr_ = fragMassDistr();

    }
    else
    {
        WarningIn("Foam::aggBreakup::setPBE()")
            << "Breakup OFF!"
            << nl << endl;
    }
}

void Foam::aggBreakup::update()
{
    if(isShearAggOn_)
    {
        // update the absolute velocity gradient field
        GA_() = sqrt(2.0) * mag(symm( fvc::grad(U_)() ));
    }

    if(isActivationOn_)
    {
        //solve transport of resting platelets
        if(isSorensenianAggOn_)
        {
            //activate platelets
            forAll(GA_->internalField(),i)
            {
                if(GA_->internalField()[i] >= activThreshold_)
                {
                    // C_RP -> C_1
//                    CMD_[0].internalField()[i] += max(Crp_->internalField()[i], 0.0);
//                    // C_RP = 0
//                    Crp_->internalField()[i] = 0.0;
                    vWF_->internalField()[i] = 2e4;
                }
            }

            volScalarField D
            (
                IOobject
                (
                    "D",
                    runTime_.timeName(),
                    mesh_
                ),
                mesh_,
                Dlist_[0]
            );
            D.internalField() += 1.05 * D.internalField() * GA_->internalField();
            fvScalarMatrix massTransport
            (
                fvm::ddt(Crp_())
              + fvm::div(phi_, Crp_(), "div(phi,C_*)")
              - fvm::laplacian
                (
                    D,
                    Crp_(),
                    "laplacian(D,C_*)"
                )
              ==
              - fvc::Sp(vWF_(),Crp_())
            );
            massTransport.solve(mesh_.solver("C_*"));
        }
        else
        {
            fvScalarMatrix massTransport
            (
                fvm::ddt(Crp_())
              + fvm::div(phi_, Crp_(), "div(phi,C_*)")
              - fvm::laplacian(Dlist_[0], Crp_(), "laplacian(D,C_*)")
            );
            massTransport.solve(mesh_.solver("C_*"));
        }
    }

    if(isPBEcoupled_)
    {
        forAll(CMD_, i)
        {
            // CMDold stores CMD that is used in the source term S().
            // Otherwise the source for one cluster would be computed
            // with values partially of step t, and partially of t+dt.
            CMDold_[i] = CMD_[i];
        }

        forAll(CMD_, i)
        {
            volScalarField& Ci = CMD_[i];
            dimensionedScalar& Di = Dlist_[i];
            fvScalarMatrix massTransport
            (
                fvm::ddt(Ci)
              + fvm::div(phi_, Ci, "div(phi,C_*)")
              - fvm::laplacian(Di, Ci, "laplacian(D,C_*)")
              ==
                S(i)
            );
            massTransport.relax();
            massTransport.solve(mesh_.solver("C_*"));
        }
    }
    else
    {
        //solve transport
        forAll(CMD_, i)
        {
            volScalarField& Ci = CMD_[i];
            dimensionedScalar& Di = Dlist_[i];

            //solve transport of resting platelets
            if(isSorensenianAggOn_)
            {
                volScalarField D
                (
                    IOobject
                    (
                        "D",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    Di
                );
                D.internalField() += 1.05 * D.internalField() * GA_->internalField();
                if(i==0)
                {
                    fvScalarMatrix massTransport
                    (
                        fvm::ddt(Ci)
                      + fvm::div(phi_, Ci, "div(phi,C_*)")
                      - fvm::laplacian(D, Ci, "laplacian(D,C_*)")
                    ==
                        fvc::Sp(vWF_(),Crp_().oldTime())
                    );
                    massTransport.solve(mesh_.solver("C_*"));

                }
                else
                {
                    fvScalarMatrix massTransport
                    (
                        fvm::ddt(Ci)
                      + fvm::div(phi_, Ci, "div(phi,C_*)")
                      - fvm::laplacian(D, Ci, "laplacian(D,C_*)")
                    );
                    massTransport.solve(mesh_.solver("C_*"));
                }

            }
            else
            {
                fvScalarMatrix massTransport
                (
                    fvm::ddt(Ci)
                  + fvm::div(phi_, Ci, "div(phi,C_*)")
                  - fvm::laplacian(Di, Ci, "laplacian(D,C_*)")
                );
                massTransport.solve(mesh_.solver("C_*"));
            }

        }

        solveAggBreakup();
    }

    //Eliminate spurious negative values.
    forAll(CMD_, i)
    {
        CMD_[i] = max
                  (
                    CMD_[i],
                    dimensionedScalar("zero", dimMoles/dimVol, 0.0)
                  );
    }

    postProc->update();
}

void Foam::aggBreakup::solveAggBreakup()
{
    scalar subDeltaT = readScalar
    (
        aggBreakupDict_->subDict("odeSolver").lookup("subDeltaT")
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

        updateKernelsInCell(celli);

        solver_->solve
        (
            0.0,
            runTime_.deltaT().value(),
            y,
            subDeltaT
        );

        // update CMD in celli
        forAll(y, dofi)
        {
            CMD_[dofi][celli] = y[dofi];
        }
    }
}


//-This is the actual PBE
Foam::tmp<volScalarField> Foam::aggBreakup::S(label clusterI) const
{
    label& i = clusterI;

    // assign auxiliary fields
    volScalarField aggregation
    (
        IOobject
        (
            "aggregation",
            runTime_.timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("dC_i/dt", dimMoles/dimVolume/dimTime, 0.0)
    );
    volScalarField breakup
    (
        IOobject
        (
            "breakup",
            runTime_.timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("dC_i/dt", dimMoles/dimVolume/dimTime, 0.0)
    );

    // aggregation PBE by Smoluchowski (1917)
    if(isBrowninanAggOn_ || isShearAggOn_)
    {
        volScalarField creation
        (
            IOobject
            (
                "creation",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("dC_i/dt", dimMoles/dimVolume/dimTime, 0.0)
        );
        volScalarField destruction
        (
            IOobject
            (
                "destruction",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("dC_i/dt", dimMoles/dimVolume/dimTime, 0.0)
        );

        if(isGridUniform_)
        {
            for(int j = 0; j < i; ++j)
            {
                creation += 0.5 * kag(i-j-1, j) * CMDold_[i-j-1] * CMDold_[j];
            }
            for(int j = 0; j < nBins_; ++j)
            {
                destruction += kag(i, j) * CMDold_[i] * CMDold_[j];
            }
        }
        else
        {
            for(int j = 0; j < nBins_; ++j)
            {
                for(int k = 0; k < nBins_; ++k)
                {
                    creation += 0.5 * chi_[j][k][i] * kag(j, k) * CMDold_[j] * CMDold_[k];
                }
            }
            for(int j = 0; j < nBins_; ++j)
            {
                destruction += kag(i, j) * CMDold_[i] * CMDold_[j];
            }
        }

        aggregation = eta_ * ( creation - destruction);
    }

    // breakup PBE by Pandya & Spielman (1982)
    if(isBreakupOn_)
    {
        volScalarField creation
        (
            IOobject
            (
                "creation",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("dC_i/dt", dimMoles/dimVolume/dimTime, 0.0)
        );
        volScalarField destruction
        (
            IOobject
            (
                "destruction",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("dC_i/dt", dimMoles/dimVolume/dimTime, 0.0)
        );

        for(int j = i + 1; j < nBins_; ++j)
        {
            creation += fragMassDistr_()[i][j] * kbr(j) * CMDold_[j];
        }

        destruction = kbr(i) * CMDold_[i];

        breakup = creation - destruction;
    }

    tmp<volScalarField> output
    (
        new volScalarField
        (
            IOobject
            (
                "S_" + to_string(vList_[i]),
                runTime_.timeName(),
                mesh_
            ),
            aggregation + breakup
        )
    );

    return output;
}

Foam::tmp<volScalarField> Foam::aggBreakup::kag(label i, label j) const
{
    tmp<volScalarField> k
    (
        new volScalarField
        (
            IOobject
            (
                "kag",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("dC_i/dt", dimVolume/dimMoles/dimTime, 0.0)
        )
    );

    if(isBrowninanAggOn_)
    {
        k->internalField() += kdBase_()[i][j];
        if(isSorensenianAggOn_)
        {
            k->internalField() += 1.05 * GA_() * kdBase_()[i][j];
        }
    }
    if(isShearAggOn_)
    {
        k->internalField() += GA_() * ksBase_()[i][j];
    }

    return k;
}

Foam::tmp<volScalarField> Foam::aggBreakup::kbr(label i) const
{
    dimensionedScalar alin("a\'", dimless/dimTime, 1.0);

    tmp<volScalarField> k
    (
        new volScalarField
        (
            IOobject
            (
                "kbr",
                runTime_.timeName(),
                mesh_
            ),
            alin * pow(GA_() / Gstar_, b_) * kbBase_()[i]
        )
    );

    return k;
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
    os.writeKeyword("Rmono") << sds.Rp_
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

    os << "v" << tab
       << "D" << endl;
    forAll(sds.vList_,i)
    {
        os << sds.vList_[i] << tab
           << sds.Dlist_[i] << endl;
    }

    return os;
}

// ************************************************************************* //
