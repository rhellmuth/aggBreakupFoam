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
  isAggCharTimeOn_(false),
  isDifAggCharTime_(false),
  isBreCharTimeOn_(false),
  isAggPecletOn_(false),
  isKabFieldOn_(false),
  isAdvCharTimeOn_(false),
  isAggDamkoehlerOn_(false),
  isZerothMomentOn_(false),
  isFirstMomentOn_(false),
  isSecondMomentOn_(false),
  isVmeanOn_(false),
  isRmeanOn_(false),
  isI0On_(false),
  isPAon_(false),
  isAvgConcentrationOn_(false)
{

    // Post processing switches
    dictionary& postProcDict(aggBreakupDict_.subDict("postprocessing"));

    isAggCharTimeOn_ = postProcDict.lookupOrDefault("isAggCharTimeOn", false);
    isBreCharTimeOn_ = postProcDict.lookupOrDefault("isBreCharTimeOn", false);
    isAdvCharTimeOn_ = postProcDict.lookupOrDefault
    (
        "isAdvCharTimeOn",
        false
    );
    isDifAggCharTime_ = postProcDict.lookupOrDefault
    (
        "isDifAggCharTime",
        false
    );
    isAggPecletOn_   = postProcDict.lookupOrDefault("isAggPecletOn", false);
    isAggDamkoehlerOn_ = postProcDict.lookupOrDefault
    (
        "isAggDamkoehlerOn",
        false
    );
    isKabFieldOn_  = postProcDict.lookupOrDefault("isKabFieldOn", false);
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
    isI0On_ = postProcDict.lookupOrDefault("isI0On", false);
    isPAon_ = postProcDict.lookupOrDefault("isAggFractionOn", false);

    if(debug)
    {
        Info << *this;
    }

    if(isZerothMomentOn_)
    {
        M0_.set
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
                M_0()
            )
        );

        Info << "Writing field M_0" << nl << endl;
        M0_->write();
    }
    if(isFirstMomentOn_)
    {
        M1_.set
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
                M_1()
            )
        );

        Info << "Writing field M_1" << nl << endl;
        M1_->write();
    }
    if(isSecondMomentOn_)
    {
        M2_.set
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
                M_2()
            )
        );

        Info << "Writing field M_2" << nl << endl;
        M2_->write();
    }
    if(isVmeanOn_)
    {
        vMean_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "v_mean",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vMean()
            )
        );

        Info << "Writing field v_mean" << nl << endl;
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
                    "R_mean",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rMean()
            )
        );

        Info << "Writing field R_mean" << nl << endl;
        rMean_->write();
    }
    if(isI0On_)
    {
        I0_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "I0",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                I0()
            )
        );

        Info << "Writing field R_mean" << nl << endl;
        I0_->write();
    }
    if(isPAon_)
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

    if(isDifAggCharTime_)
    {
        td_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "t_d",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                t_d()
            )
        );

        Info << "Writing field t_d" << nl << endl;
        td_->write();
    }

    if(isAggCharTimeOn_ || isBreCharTimeOn_ || isKabFieldOn_)
    {
        const volScalarField& C_0(CMD[0]);

        // Check whether strain rate tensor is not on objectRegistry
        if(C_0.db().foundObject<volScalarField>(word("G")))
        {            
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
                        t_s()
                    )
                );

                Info << "Writing field t_a" << nl << endl;
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

            if(isKabFieldOn_)
            {
                Kab_.set
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "K_ab",
                            runTime_.timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        K_ab()
                    )
                );

                Info << "Writing field K_ab" << nl << endl;
                Kab_->write();
            }

        }
        else
        {
            WarningIn("Foam::aggBreakupPostprocess::aggBreakupPostprocess")
                << "Characteristic times dependent on absolute shear rate (G) "
                << "cannot be calculated because G is not defined on the "
                << "objectRegistry"
                << nl << endl;
        }


    }

    if(isAggPecletOn_)
    {
        Pe_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "Pe",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Pe()
            )
        );

        Info << "Writing field Pe"
             << nl << endl;
        Pe_->write();
    }
    if(isAvgConcentrationOn_)
    {

    }

    if(isAdvCharTimeOn_)
    {
        tadv_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "t_adv",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                t_adv()
            )
        );

        Info << "Writing field t_adv (advection characteristic time)"
             << nl << endl;
        tadv_->write();
    }

    if(isAggDamkoehlerOn_)
    {
        Da_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "Da",
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                Da()
            )
        );

        Info << "Writing field Da"
             << nl << endl;
        Da_->write();
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

Foam::tmp<volScalarField> Foam::aggBrePostprocess::G() const
{
    // update shear rate
    const volScalarField& C_0(CMD_[0]);

    const volScalarField& G(C_0.db().lookupObject<volScalarField>(word("G")));
    return G;
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::M_0() const
{
    tmp<volScalarField> output
    (
        new volScalarField
        (
            IOobject
            (
                "summation",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("summation", dimMoles/dimVolume, 0.0)
        )
    );

    bool isActivationOn = readBool(aggBreakupDict_.lookup("isActivationOn"));
    if(isActivationOn)
    {
        const volScalarField& C_0(CMD_[0]);
        const volScalarField& Crp(C_0.db().lookupObject<volScalarField>(word("C_RP")));
        output() += Crp;
    }

    forAll(CMD_, i)
    {
        output() += CMD_[i];
    }

    return output;
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::M_1() const
{
    tmp<volScalarField> output
    (
        new volScalarField
        (
            IOobject
            (
                "summation",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("summation", dimMoles/dimVolume, 0.0)
        )
    );

    bool isActivationOn = readBool(aggBreakupDict_.lookup("isActivationOn"));
    if(isActivationOn)
    {
        const volScalarField& C_0(CMD_[0]);
        const volScalarField& Crp(C_0.db().lookupObject<volScalarField>(word("C_RP")));
        output() += Crp;
    }

    forAll(CMD_, i)
    {
        output() += CMD_[i] * vList_[i];
    }
    return output;
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::M_2() const
{
    tmp<volScalarField> output
    (
        new volScalarField
        (
            IOobject
            (
                "summation",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("summation", dimMoles/dimVolume, 0.0)
        )
    );

    bool isActivationOn = readBool(aggBreakupDict_.lookup("isActivationOn"));
    if(isActivationOn)
    {
        const volScalarField& C_0(CMD_[0]);
        const volScalarField& Crp(C_0.db().lookupObject<volScalarField>(word("C_RP")));
        output() += Crp;
    }

    forAll(CMD_, i)
    {
        output() += CMD_[i] * pow(vList_[i], 2.0);
    }

    return output;
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::vMean() const
{
    return
    (
        M_1() / (M_0() + dimensionedScalar("VSMALL", dimMoles/dimVolume, VSMALL))
    );
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

    bool isActivationOn = readBool(aggBreakupDict_.lookup("isActivationOn"));
    if(isActivationOn)
    {
        const volScalarField& C_0(CMD_[0]);
        const volScalarField& Crp(C_0.db().lookupObject<volScalarField>(word("C_RP")));
        RiCi() += Crp * rList_[0];
    }

    forAll(CMD_, i)
    {
        RiCi() += CMD_[i] * rList_[i];
    }

    return
    (
        RiCi() / (M_0() + dimensionedScalar("VSMALL", dimMoles/dimVolume, VSMALL))
    );
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::I0() const
{
    return
    (
        M_2() / (M_1() + dimensionedScalar("VSMALL", dimMoles/dimVolume, VSMALL))
    );
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::PA() const
{
    bool isActivationOn = readBool(aggBreakupDict_.lookup("isActivationOn"));
    if(isActivationOn)
    {
        const volScalarField& C_0(CMD_[0]);
        const volScalarField& Crp(C_0.db().lookupObject<volScalarField>(word("C_RP")));

        return
        (
            100.0 *
            (
                1.0 -
                (
                    (CMD_[0] + Crp) /
                    (
                        M_1() +
                        dimensionedScalar("VSMALL", dimMoles / dimVolume, VSMALL)
                    )
                )
            )
        );
    }
    else
    {
        return
        (
            100.0 *
            (
                1.0 -
                (
                    CMD_[0] /
                    (
                        M_1() +
                        dimensionedScalar("VSMALL", dimMoles / dimVolume, VSMALL)
                    )
                )
            )
        );
    }
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::t_s() const
{
    scalar eta(readScalar(aggBreakupDict_.lookup("eta")));
    //alpha is actually dimless but we are using dimMoles to account platelets
    dimensionedScalar alpha("alpha", dimless/dimMoles, 32./3.);
    const dimensionedScalar& R_p(rList_[0]);

    return
    (
        min
        (
            1.0 /
            (
                eta * alpha * G() * pow(R_p, 3.) * M_1() +
                dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
            ),
            dimensionedScalar("GREAT", dimTime, GREAT)
        )
    );
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::t_b() const
{
    dimensionedScalar alin("alin", dimless/dimTime, 1.0);
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

        Gstar = alin / (a * pow(R_mono, c) + VSMALL);
    }

    return
    (
        1.0 /
        (
            alin * pow(G()/Gstar, b) +
            dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
        )
    );
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::K_ab() const
{
    return
    (
        t_b() / (t_s() + dimensionedScalar("VSMALL", dimTime, VSMALL))
    );
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::t_adv() const
{
    tmp<volScalarField> charTime
    (
        new volScalarField
        (
            IOobject
            (
                "charTime",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("charTime", dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>("phi");

    charTime->internalField() =
              2.0 * mesh_.V().field()
            / (fvc::surfaceSum(mag(phi))().internalField() + VSMALL);

    charTime->correctBoundaryConditions();

    return charTime;
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::Da() const
{
    return
    (
        t_adv() / (t_s() + dimensionedScalar("VSMALL", dimTime, VSMALL))
    );
}


Foam::tmp<volScalarField> Foam::aggBrePostprocess::t_d() const
{
    dimensionedScalar T = aggBreakupDict_.lookup("T");
    dimensionedScalar muPlasma = aggBreakupDict_.lookup("muPlasma");
    //Boltzmann constant [J K−1]
    dimensionedScalar kB("k_B", dimEnergy / dimTemperature, 1.3806488E-23);
    const dimensionedScalar& R_mono(rList_[0]);

    dimensionedScalar D_mono = kB * T / (6.0 * M_PI * muPlasma * R_mono);

    return
    (
        min
        (
            // one is actually dimless but we are using dimMoles (amount of
            // substance) to account platelets
            dimensionedScalar("one", dimMoles, 1.0) /
            (
                4.0 * M_PI * (2.0 * R_mono) * (2.0 * D_mono) * M_1() +
                dimensionedScalar("VSMALL", dimMoles/dimTime, VSMALL)
            ),
            dimensionedScalar("GREAT", dimTime, GREAT)
        )
    );
}

Foam::tmp<volScalarField> Foam::aggBrePostprocess::Pe() const
{
    return
    (
        t_d() / (t_s() + dimensionedScalar("VSMALL", dimTime, VSMALL))
    );
}

void Foam::aggBrePostprocess::update()
{
    if(runTime_.outputTime()) //calculate only at outputTime
    {
        if(M0_.valid())
        {
            M0_() = M_0();
        }

        if(M1_.valid())
        {
            M1_() = M_1();
        }

        if(M2_.valid())
        {
            M2_() = M_2();
        }

        if(vMean_.valid())
        {
            vMean_() = vMean();
        }

        if(rMean_.valid())
        {
            rMean_() = rMean();
        }

        if(I0_.valid())
        {
            I0_() = I0();
        }

        if(PA_.valid())
        {
            PA_() = PA();
        }

        if(td_.valid())
        {
            td_() = t_d();
        }

        const volScalarField& C_0(CMD_[0]);

        // Check whether strain rate tensor is not on objectRegistry
        if(C_0.db().foundObject<volScalarField>(word("G")))
        {
            if(ta_.valid())
            {
                ta_() = t_s();
            }

            if(tb_.valid())
            {
                tb_() = t_b();
            }

            if(Kab_.valid())
            {
                Kab_() = K_ab();
            }

            if(tadv_.valid())
            {
                tadv_() = t_adv();
            }

            if(Da_.valid())
            {
                Da_() = Da();
            }
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
    os.writeKeyword("isKabFieldOn") << sds.isKabFieldOn_
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
