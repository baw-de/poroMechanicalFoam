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

#include "montenegroMaier.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace saturationLaws
    {
        defineTypeNameAndDebug(montenegroMaier, 0);

        addToRunTimeSelectionTable(
            saturationLaw,
            montenegroMaier,
            dictionary);

        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        scalar montenegroMaier::Sfgr(const scalar p,const scalar Sfgr0) const
        {
            return Sfgr0*(p_At_.value()/max(p_At_.value()+p,0.95*p_At_.value()));
        }

        scalar montenegroMaier::Sfw(const scalar SfgrScalar) const
        {
            return 1.0 - SfgrScalar;
        }

        scalar montenegroMaier::Sf(const scalar &HSf, const scalar S0, const scalar Sr) const
        {
            return Sr + (S0 - Sr)*HSf;
        }

        scalar montenegroMaier::SFunc(const scalar SfwScalar, const scalar SfScalar) const
        {
            return SfwScalar * SfScalar;
        } 

        scalar montenegroMaier::Ss_gr(const scalar p, const label cellI) const
        {
            scalar Sfgr0 = S_fgr0_.internalField()[cellI];
            return Sfgr0*(p_At_.value()/pow(max(p_At_.value()+p,0.95*p_At_.value()),2));
        }

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        montenegroMaier::montenegroMaier(
            const word &name,
            dictionary &poroHydraulicProperties,
            const volScalarField &pField)
            : vanGenuchten(name, poroHydraulicProperties, pField),
              montenegroMaierCoeffs_(poroHydraulicProperties.subDict(typeName + "Coeffs")),
            S_fgr0_(
              IOobject(
                  "S_fgr0",
                  mesh().time().timeName(),
                  db(),
                  IOobject::READ_IF_PRESENT,
                  IOobject::NO_WRITE),
              mesh(),
              montenegroMaierCoeffs_.get<dimensionedScalar>("S_fgr0")),
            p_At_(montenegroMaierCoeffs_.get<dimensionedScalar>("p_At")),
            Sf_(
              IOobject(
                  "Sf",
                  mesh().time().timeName(),
                  db(),
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE),
              mesh(),
              dimensionedScalar("",dimless,1.0),
                  "zeroGradient"),
            Sfw_(
              IOobject(
                  "Sfw",
                  mesh().time().timeName(),
                  db(),
                  IOobject::NO_READ,
                  IOobject::NO_WRITE),
              mesh(),
              dimensionedScalar("",dimless,1.0),
                  "zeroGradient"),
            Ss_gr_(
              IOobject(
                  "Ss_gr",
                  mesh().time().timeName(),
                  db(),
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE),
              mesh(),
              dimensionedScalar("",dimless/p_At_.dimensions(),1.0),
                  "zeroGradient")
        {
           if (p_At_.dimensions()!=pField.dimensions() && dimless/alpha_.dimensions()!=pField.dimensions())
            {
              FatalErrorIn("montenegroMaier::montenegroMaier") << "Pressure dimensions are: " << pField.dimensions() << nl
                                                    << "p_At has (wrong?) dimensions " << p_At_.dimensions() << nl
                                                    << "alpha has (wrong?) dimensions " << alpha_.dimensions() << endl;

            } 
            if (debug)
            {
                S_fgr0_.write();
            }
        }

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


        scalar montenegroMaier::C(const scalar p, const label cellI)
        {
            scalar C_vg = vanGenuchten::C(p,cellI);
            Ss_gr_.internalFieldRef()[cellI] = Ss_gr(p,cellI);

            scalar Sfgr0 = S_fgr0_.internalField()[cellI];
            scalar SfgrI = Sfgr(p,Sfgr0);
            scalar SfwI = Sfw(SfgrI);

            const scalar nS = n_.internalField()[cellI];
            const scalar mS = m_.internalField()[cellI];
            const scalar alphaS = alpha_.internalField()[cellI];
            const scalar HS = vanGenuchten::H(p,nS,mS,alphaS);

            scalar Sr = S_r.internalField()[cellI];
            scalar S0 = S_0.internalField()[cellI];
            scalar SfI = Sf(HS,S0,Sr);

            return SfwI *(C_vg + (SfI * Ss_gr_.internalFieldRef()[cellI]));
        };
 
        scalar montenegroMaier::S(const scalar p, const label cellI)
        {
            scalar Sfgr0 = S_fgr0_.internalField()[cellI];
            scalar SfgrI = Sfgr(p,Sfgr0);
            scalar SfwI = Sfw(SfgrI);
            // For output 
            Sfw_.internalFieldRef()[cellI] = SfwI;

            const scalar nS = n_.internalField()[cellI];
            const scalar mS = m_.internalField()[cellI];
            const scalar alphaS = alpha_.internalField()[cellI];
            const scalar HS = vanGenuchten::H(p,nS,mS,alphaS);

            scalar Sr = S_r.internalField()[cellI];
            scalar S0 = S_0.internalField()[cellI];
            scalar SfI = Sf(HS,S0,Sr);
            // For output 
            Sf_.internalFieldRef()[cellI] = SfI;

            return SFunc(SfwI,SfI);
        }  

        scalar montenegroMaier::S
        (
            const scalar p,
            const label patchI,
            const label faceI
        )
        {
            scalar Sfgr0 = S_fgr0_.boundaryField()[patchI][faceI];
            scalar SfgrI = Sfgr(p,Sfgr0);
            scalar SfwI = Sfw(SfgrI);
            // For output 
            Sfw_.boundaryFieldRef()[patchI][faceI] = SfwI;

            const scalar nS = n_.boundaryField()[patchI][faceI];
            const scalar mS = m_.boundaryField()[patchI][faceI];
            const scalar alphaS = alpha_.boundaryField()[patchI][faceI];
            const scalar HS = vanGenuchten::H(p,nS,mS,alphaS);
            
            scalar Sr = S_r.boundaryField()[patchI][faceI];
            scalar S0 = S_0.boundaryField()[patchI][faceI];
            scalar SfI = Sf(HS,S0,Sr);
            // For output 
            Sf_.boundaryFieldRef()[patchI][faceI] = SfI;

            return SFunc(SfwI,SfI);
        }  
        


    } // End of namespace saturationLaws
} // End of namespace Foam

//*********************************************************** //
