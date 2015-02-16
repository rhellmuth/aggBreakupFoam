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

#include "aggBrePostprocess.H"

namespace Foam
{
    defineTypeNameAndDebug(aggBrePostprocess, 0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aggBrePostprocess::aggBrePostprocess
(
    const dictionary& aggBreakupDict,
    const PtrList<volScalarField>& CMD,
    const PtrList<scalar>& vList,
    const PtrList<dimensionedScalar>& rList
)
:
  aggBreakupDict_(aggBreakupDict),
  CMD_(CMD),
  runTime_(CMD_[0].time()),
  mesh_(CMD_[0].mesh()),
  vList_(vList),
  rList_(rList),
  isThetaFieldOn_(false),
  isAggCharTimeOn_(false),
  isBreCharTimeOn_(false),
  isAggPecletOn_(false),
  isAvgConcentrationOn_(false),
  isZerothMomentOn_(false),
  isFirstMomentOn_(false),
  isSecondMomentOn_(false),
  isVmeanOn_(false),
  isRmeanOn_(false),
  isAggFractionOn_(false)
{

    // Post processing switches
    dictionary& postProcDict(aggBreakupDict_.subDict("postprocessing"));

    isThetaFieldOn_  = postProcDict.lookupOrDefault("isThetaFieldOn", false);
    isAggCharTimeOn_ = postProcDict.lookupOrDefault("isAggCharTimeOn", false);
    isBreCharTimeOn_ = postProcDict.lookupOrDefault("isBreCharTimeOn", false);
    isAggPecletOn_   = postProcDict.lookupOrDefault("isAggPecletOn", false);
    isAvgConcentrationOn_ = postProcDict.lookupOrDefault
                          (
                              "isAvgConcentrationOn",
                              false
                          );
    isZerothMomentOn_ = postProcDict.lookupOrDefault("isZerothMomentOn", false);
    isFirstMomentOn_ = postProcDict.lookupOrDefault("isFirstMomentOn", false);
    isSecondMomentOn_ = postProcDict.lookupOrDefault("isSecondMomentOn", false);
    isVmeanOn_ = postProcDict.lookupOrDefault("isVmeanOn", false);
    isRmeanOn_ = postProcDict.lookupOrDefault("isRmeanOn", false);
    isAggFractionOn_ = postProcDict.lookupOrDefault("isAggFractionOn", false);

    if(debug)
    {
        Info << *this;
    }


    if(isZerothMomentOn_)
    {
        zerothMoment_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "M_0",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                zerothMoment()
            )
        );

        Info << "Writing field M_0" << nl << endl;
        zerothMoment_->write();
    }
    if(isFirstMomentOn_)
    {
        firstMoment_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "M_1",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                firstMoment()
            )
        );

        Info << "Writing field M_1" << nl << endl;
        firstMoment_->write();
    }
    if(isSecondMomentOn_)
    {
        secondMoment_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "M_2",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                secondMoment()
            )
        );

        Info << "Writing field M_2" << nl << endl;
        secondMoment_->write();
    }

    if(isVmeanOn_)
    {
        vMean_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "vMean",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vMean()
            )
        );

        Info << "Writing field vMean" << nl << endl;
        vMean_->write();
    }

    if(isRmeanOn_)
    {
        rMean_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "rMean",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rMean()
            )
        );

        Info << "Writing field rMean" << nl << endl;
        rMean_->write();
    }

    if(isAggFractionOn_)
    {
        PA_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "PA",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                PA()
            )
        );

        Info << "Writing field PA" << nl << endl;
        PA_->write();
    }

    if(isAggCharTimeOn_ || isBreCharTimeOn_ || isThetaFieldOn_)
    {
        const volScalarField& C_0(CMD[0]);

        // Check whether strain rate tensor is not on objectRegistry
        if(C_0.db().foundObject<volTensorField>(word("E")))
        {
            const volTensorField& E
            (
                C_0.db().lookupObject<volTensorField>(word("E"))
            );
            G_.set
            (
                new volScalarField
                (
                    IOobject
                    (
                        "G",
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mag(E)/ sqrt(2.0)
                )
            );

            G_->write();

            if(isAggCharTimeOn_)
            {
                ta_.set
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "t_a",
                            runTime_.timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        t_a()
                    )
                );
                Info << "Writing field t_a" << endl;
                ta_->write();
            }

            if(isBreCharTimeOn_)
            {
                tb_.set
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "t_b",
                            runTime_.timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        t_b()
                    )
                );

                Info << "Writing field t_b" << nl << endl;
                tb_->write();
            }

            if(isThetaFieldOn_)
            {
                theta_.set
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "theta",
                            runTime_.timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        theta()
                    )
                );

                Info << "Writing field theta" << nl << endl;
                theta_->write();
            }

        }
        else
        {
            WarningIn("Foam::aggBreakupPostprocess::aggBreakupPostprocess")
                << "Characteristic times dependent on strain rate (E) "
                << "cannot be calculated because E is not defined on the "
                << "objectRegistry"
                << nl << endl;
        }


    }

    if(isAggPecletOn_)
    {

    }
    if(isAvgConcentrationOn_)
    {

    }


}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::aggBrePostprocess>
//Foam::aggBrePostprocess::New()
//{
//    return autoPtr<aggBrePostprocess>(new aggBrePostprocess);
//}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::aggBrePostprocess::~aggBrePostprocess()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<volScalarField> Foam::aggBrePostprocess::zerothMoment() const
{
    tmp<volScalarField> M_0
    (
        new volScalarField
        (
            IOobject
            (
                "M_0",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("M_0", dimMoles/dimVolume, 0.0)
        )
    );

    forAll(CMD_, i)
    {
        M_0() += CMD_[i];
    }
    return M_0;
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::firstMoment() const
{
    tmp<volScalarField> M_1
    (
        new volScalarField
        (
            IOobject
            (
                "M_1",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("M_1", dimMoles/dimVolume, 0.0)
        )
    );

    forAll(CMD_, i)
    {
        M_1() += CMD_[i] * vList_[i];
    }
    return M_1;
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::secondMoment() const
{
    tmp<volScalarField> M_2
    (
        new volScalarField
        (
            IOobject
            (
                "M_2",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("M_2", dimMoles/dimVolume, 0.0)
        )
    );

    forAll(CMD_, i)
    {
        M_2() += CMD_[i] * pow(vList_[i], 2.0);
    }
    return M_2;
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::vMean() const
{
    return firstMoment() / zerothMoment();
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::rMean() const
{
    tmp<volScalarField> RiCi
    (
        new volScalarField
        (
            IOobject
            (
                "RiCi",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("RiCi", dimLength * dimMoles / dimVolume, 0.0)
        )
    );

    forAll(CMD_, i)
    {
        RiCi() += CMD_[i] * rList_[i];
    }

    return RiCi() / zerothMoment();
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::PA() const
{
    return 100.0 * (1.0 - CMD_[0] / (firstMoment() +
                dimensionedScalar("VSMALL", dimMoles / dimVolume, VSMALL)));
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::t_a() const
{
    scalar eta(readScalar(aggBreakupDict_.lookup("eta")));
    dimensionedScalar alpha("alpha", dimless/dimMoles, 4./3.);
    const dimensionedScalar& R_mono(rList_[0]);

    return 1.0 / (eta * alpha * G_() * pow(R_mono, 3.) * firstMoment()
                  + dimensionedScalar("VSMALL", dimless/dimTime, VSMALL));
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::t_b() const
{
    scalar b(readScalar(aggBreakupDict_.lookup("b")));
    dimensionedScalar Gstar("Gstar", dimless/dimTime, 1000.0);

    if(aggBreakupDict_.found("Gstar"))
    {
        Gstar = aggBreakupDict_.lookup("Gstar");
    }
    else if(aggBreakupDict_.found("a") && aggBreakupDict_.found("c"))
    {
        const dimensionedScalar& R_mono(rList_[0]);
        scalar a(readScalar(aggBreakupDict_.lookup("a")));
        scalar c(readScalar(aggBreakupDict_.lookup("c")));
        dimensionedScalar alin("alin", dimless/dimTime, 1.0);

        Gstar = alin / (a * pow(R_mono, c) + VSMALL);
    }

    return 1.0 / (pow( G_()/Gstar + VSMALL, b) + VSMALL);
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::theta() const
{
    return t_a() / (t_b() + VSMALL);
}


void Foam::aggBrePostprocess::update()
{
    if(zerothMoment_.valid())
    {
        zerothMoment_() = zerothMoment();
    }

    if(firstMoment_.valid())
    {
        firstMoment_() = firstMoment();
    }

    if(secondMoment_.valid())
    {
        secondMoment_() = secondMoment();
    }

    if(vMean_.valid())
    {
        vMean_() = vMean();
    }

    if(rMean_.valid())
    {
        rMean_() = rMean();
    }

    if(G_.valid())
    {
        // update shear rate
        const volScalarField& C_0(CMD_[0]);
        const volTensorField& E
        (
            C_0.db().lookupObject<volTensorField>(word("E"))
        );
        G_() = mag(E) / sqrt(2.0);

        if(ta_.valid())
        {
            ta_() = t_a();
        }

        if(tb_.valid())
        {
            tb_() = t_b();
        }

        if(theta_.valid())
        {
            theta_() = theta();
        }
    }

}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::aggBrePostprocess::operator=(const aggBrePostprocess& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::aggBrePostprocess::operator=(const Foam::aggBrePostprocess&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const aggBrePostprocess& sds)
{
    os << endl << "Post-processing switches:" << endl;
    os.writeKeyword("isThetaFieldOn") << sds.isThetaFieldOn_
                                      << token::END_STATEMENT
                                      << endl;
    os.writeKeyword("isAggCharTimeOn") << sds.isAggCharTimeOn_
                                       << token::END_STATEMENT
                                       << endl;
    os.writeKeyword("isBreCharTimeOn") << sds.isBreCharTimeOn_
                                       << token::END_STATEMENT
                                       << endl;
    os.writeKeyword("isAggPecletOn") << sds.isAggPecletOn_
                                     << token::END_STATEMENT
                                     << endl;
    os.writeKeyword("isAvgConcentrationOn") << sds.isAvgConcentrationOn_
                                            << token::END_STATEMENT
                                            << endl;
    os.writeKeyword("isZerothMomentOn") << sds.isZerothMomentOn_
                                        << token::END_STATEMENT
                                        << endl;
    os.writeKeyword("isFirstMomentOn") << sds.isFirstMomentOn_
                                       << token::END_STATEMENT
                                       << endl;
    os.writeKeyword("isSecondMomentOn") << sds.isSecondMomentOn_
                                        << token::END_STATEMENT
                                        << endl;
    os.writeKeyword("isVmeanOn") << sds.isVmeanOn_
                                 << token::END_STATEMENT
                                 << endl;
    os.writeKeyword("isRmeanOn") << sds.isRmeanOn_
                                 << token::END_STATEMENT
                                 << endl;

    return os;
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
