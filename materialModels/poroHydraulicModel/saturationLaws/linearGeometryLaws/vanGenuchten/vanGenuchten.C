/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "vanGenuchten.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {
        defineTypeNameAndDebug(vanGenuchten, 0);

        addToRunTimeSelectionTable(
            saturationLaw,
            vanGenuchten,
            dictionary);

        // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //


        scalar vanGenuchten::SFunc(const scalar HS, const scalar S0, const scalar Sr)
        {
            return HS * (S0 - Sr) + Sr;
        }

        scalar vanGenuchten::krFunc(const scalar p, const scalar Hk, const scalar mk)
        {
            return max
            (
                neg(p) * sqrt(Hk) * pow(pow(1 - pow(Hk, 1.0 / mk), mk) - 1, 2) + pos(p),
                ROOTVSMALL
            );
        }

        scalar vanGenuchten::H
        (
            const scalar p,
            const scalar nH,
            const scalar mH,
            const scalar alphaH
        )
        {
            scalar pTilda = max(mag(p),SMALL);
            return neg(p) * pow(1.0 + pow(alphaH * mag(pTilda), nH), -mH) + pos(p);
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        vanGenuchten::vanGenuchten(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField)
            : saturationLaw(name, poroHydraulicProperties, pField),
              vanGenuchtenCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
              n_(
                  IOobject(
                      "n_swcc",
                      mesh().time().timeName(),
                      db(),
                      IOobject::READ_IF_PRESENT,
                      IOobject::NO_WRITE),
                  pField.mesh(),
                  dimensionedScalar(vanGenuchtenCoeffs_.lookup("n"))),
              m_(
                  IOobject(
                      "m_",
                      mesh().time().timeName(),
                      db(),
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  1.0 - 1.0 / n_),
              S_r(
                  IOobject(
                      "S_r",
                      mesh().time().timeName(),
                      db(),
                      IOobject::READ_IF_PRESENT,
                      IOobject::NO_WRITE),
                  pField.mesh(),
                  dimensionedScalar(vanGenuchtenCoeffs_.lookup("S_r"))),
              S_0(
                  IOobject(
                      "S_0",
                      mesh().time().timeName(),
                      db(),
                      IOobject::READ_IF_PRESENT,
                      IOobject::NO_WRITE),
                  pField.mesh(),
                  dimensionedScalar(vanGenuchtenCoeffs_.lookup("S_0"))),
              alpha_(
                  IOobject(
                      "alpha",
                      mesh().time().timeName(),
                      db(),
                      IOobject::READ_IF_PRESENT,
                      IOobject::NO_WRITE),
                  pField.mesh(),
                  vanGenuchtenCoeffs_.get<dimensionedScalar>("alpha"))
        {
           if ((dimless/alpha_.dimensions())!=pField.dimensions())
            {
              FatalErrorIn("vanGenuchten::vanGenuchten") << "Pressure dimensions are: " << pField.dimensions() << endl
                                                    << "alpha has wrong dimensions " << alpha_.dimensions() << endl;
            }

            if (debug)
            {
                n_.write();
                m_.write();
                S_r.write();
                S_0.write();
                alpha_.write();
            }
        }

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        scalar vanGenuchten::pStar(const label cellI) const
        {
            // pStar = 1/alpha_*(n-1)/pow(m*n+1,1/n);
            scalar alphaP = alpha_.internalField()[cellI];
            scalar nP = n_.internalField()[cellI];    

            return 
            (
                -1 / alphaP
                * pow
                    (
                    (nP - 1) / nP,
                    1 / nP
                    )
            );
        }

        scalar vanGenuchten::C(const scalar p, const label cellI)
        {
            const scalar nC = n_.internalField()[cellI];
            const scalar mC = m_.internalField()[cellI];
            const scalar alphaC = alpha_.internalField()[cellI];
            const scalar HC = H(p,nC,mC,alphaC);
            scalar S0 = S_0.internalField()[cellI];
            scalar Sr = S_r.internalField()[cellI];

            return  
            (   
                alphaC * mC * (S0 - Sr) / (1.0 - mC) * pow(HC, 1.0 / mC) 
                * pow
                  (
                    1.0 - pow(HC, 1.0 / mC),
                    mC
                  )
            );
        }

        scalar vanGenuchten::S(const scalar p, const label cellI)
        {
            scalar S0 = S_0.internalField()[cellI];
            scalar Sr = S_r.internalField()[cellI];
            const scalar nS = n_.internalField()[cellI];
            const scalar mS = m_.internalField()[cellI];
            const scalar alphaS = alpha_.internalField()[cellI];
            const scalar HS = H(p,nS,mS,alphaS);
            return SFunc(HS,S0,Sr);
        }

        scalar vanGenuchten::S
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            scalar S0 = S_0.boundaryField()[patchI][faceI];
            scalar Sr = S_r.boundaryField()[patchI][faceI];
            const scalar nS = n_.boundaryField()[patchI][faceI];
            const scalar mS = m_.boundaryField()[patchI][faceI];
            const scalar alphaS = alpha_.boundaryField()[patchI][faceI];
            const scalar HS = H(p,nS,mS,alphaS);
            return SFunc(HS,S0,Sr);
        }

        scalar vanGenuchten::kr(const scalar p, const label cellI)
        {
            const scalar nk = n_.internalField()[cellI];
            const scalar mk = m_.internalField()[cellI];
            const scalar alphak = alpha_.internalField()[cellI];
            const scalar Hk = H(p,nk, mk, alphak);

            return krFunc(p,Hk,mk);
        }

        scalar vanGenuchten::kr
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            const scalar nk = n_.boundaryField()[patchI][faceI];
            const scalar mk = m_.boundaryField()[patchI][faceI];
            const scalar alphak = alpha_.boundaryField()[patchI][faceI];
            const scalar Hk = H(p,nk, mk, alphak);

            return krFunc(p,Hk,mk);
        }

        /*scalar vanGenuchten::dkbydp(const scalar p, const label cellI)
        {
            const volScalarField ap = alpha_*mag(p);
            const volScalarField apn = pow(ap,n_);
            const volScalarField partone = neg(p)*(n_-1)*alpha_*pow(ap,n_-1)*pow(1+apn,-5*m_/2)*(pow(1+apn,m_)-pow(apn,m_));
            const volScalarField parttwo = (apn*(5./2.*apn/(1+apn)+2)-0.5*pow(1+apn,m_+1));

           dkbydp().ref() = - partone * parttwo;

           return dkbydp();
        }*/




    } // End of namespace saturationLaws
} // End of namespace Foam

//*********************************************************** //
