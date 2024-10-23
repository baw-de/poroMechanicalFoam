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

#include "Casulli.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    namespace richardsLinearizations
    {
        defineTypeNameAndDebug(Casulli, 0);

        addToRunTimeSelectionTable(
            richardsLinearization,
            Casulli,
            dictionary);

        // * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        Casulli::Casulli(
            const word &name,
            varSatPoroHydraulicModel &poroHydraulic,
            dictionary &poroFluidProperties,
            volScalarField& S)
            : richardsLinearization(name, poroHydraulic, poroFluidProperties, S),
              casulliDict_(poroFluidProperties.subDictOrAdd("Casulli")),
              r_max_(casulliDict_.lookupOrAddDefault<scalar>("CasulliResidual", 1e-6)),
              nCasInt_(casulliDict_.lookupOrAddDefault<scalar>("CasulliInnerCorrectors", 1000)),
              nCasExt_(casulliDict_.lookupOrAddDefault<scalar>("CasulliOuterCorrectors", 1000)),
                C_pStar_(
                    IOobject(
                        "C_pStar",
                        mesh_.time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE),
                    mesh_,
                    dimensionedScalar("0",C().dimensions(), 0.0),
                "zeroGradient"),
                S_pStar_(
                    IOobject(
                        "S_pStar",
                        mesh_.time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE),
                    mesh_,
                    dimensionedScalar("0", dimless, 0.0),
                "zeroGradient"),
                pStar_(
                    IOobject(
                        "pStar",
                        mesh_.time().timeName(),
                        mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE),
                    mesh_,
                    dimensionedScalar("0", dimless/C_pStar_.dimensions(), 0.0)),
              P_(),
              Q_(),
              S1_(),
              S2_(),
              f1_(),
              f2_(),
              pPrevInt_(),
              pPrevExt_(),
              extCorr_(0),
              intCorr_(0),
              totCorr_(0),
              storedTimeIndex_(0)
        {
            update_pStar(poroHydraulic.pStar());
        }

        // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

        void Casulli::initalize(volScalarField &totalP,volScalarField &pField)
        {
            ///////////////// PORE-FLUID CONTINUITY ////////////////////////////////////////////////////
            C() = poroHydraulic().C(totalP);
            
            P(totalP,C());

            Q(totalP,P_());

            S() = poroHydraulic().S(totalP);

            S1(totalP,S());

	          S1_.ref().storePrevIter();

            S2(totalP,S1_());

	          S2_.ref().storePrevIter();

            f1_.reset(new volScalarField("f1",S1_() - P_() * pField));

            f2_.reset(new volScalarField("f2",S2_() - Q_() * pField));

            pPrevInt_.reset(new volScalarField("pCasPrevIntIter",pField));
            pPrevExt_.reset(new volScalarField("pCasPrevExtIter",pField));

        // Requirement for Initial guess at the beginning of each timestep
            //if(mesh().time().timeIndex()!=storedTimeIndex_)
            //{
                totalP = min(totalP, pStar_);
                //storedTimeIndex_ = mesh().time().timeIndex();
            //}

        }

        bool Casulli::checkConvergedAndUpdate(volScalarField &totalP, volScalarField &pField)
        {
            volScalarField S1_Prev("S1Prev",S1_());
            S() = poroHydraulic().S(totalP);
            S1(totalP,S());
            scalar rint = L2Norm((P_() * (pPrevInt_() - pField) - (S1_Prev - S1_()))().primitiveField()); // Residual as in Casulli
            if (rint < r_max_ || intCorr_ > nCasInt_)
                {
                    intCorr_ = 0;
                    volScalarField S2_Prev("S2Prev",S2_());
                    S2(totalP,S1_());
                    scalar rext = L2Norm(-(Q_() * (pPrevExt_() - pField) - (S2_Prev - S2_()))().primitiveField());
                        Info << "Casulli:" << tab 
                             << "Iterations (total) =  " << totCorr_+1 << tab
                             << "S_1 Residual = " << rint << tab  
                             << "S_2 Residual = " << rext << endl;
                    if (rext < r_max_ || extCorr_ > nCasExt_)
                    {
                    	P_.clear();
                    	Q_.clear();
                    	S1_.clear();
                    	S2_.clear();
                    	f1_.clear();
                    	f2_.clear();
                    	pPrevInt_.clear();
                    	pPrevExt_.clear();
                        extCorr_ = 0;
                        totCorr_ = 0;
                        return true;
                    }
                    else
                    {
                        C() = poroHydraulic().C(totalP);
                        P(totalP,C()); // First parameter from Jordan-Decomposition of C(h)
                        f1_.ref() = S1_() - P_() * pField;                                 // RHS of Jordan-Decomposition (f^kl)
                        pPrevInt_.ref() = pField; // store PRefrgh^(k-1,l-1)
                        intCorr_++;
                        totCorr_++;
                        Q(totalP,P_()); // Second parameter from Jordan-Decomposition of C(h)
                        f2_.ref() = S2_() - Q_() * pField;                               // RHS of Jordan-Decomposition (d^k)
                        pPrevExt_.ref() = pField; // store PRefrgh^(k-1)
                        extCorr_++;
                        return false;
                    }
                }
            C() = poroHydraulic().C(totalP);
            P(totalP,C()); // First parameter from Jordan-Decomposition of C(h)
            f1_.ref() = S1_() - P_() * pField;                                 // RHS of Jordan-Decomposition (f^kl)
            pPrevInt_.ref() = pField; // store PRefrgh^(k-1,l-1)
            intCorr_++;
            totCorr_++;
            return false;
        }
                
        tmp<fvScalarMatrix> Casulli::ddtS(const volScalarField &S, volScalarField &pField)
        {
            tmp<fvScalarMatrix> ddtSTMP( new fvScalarMatrix(
                pField,
                P_().dimensions() * pField.dimensions() * dimVol / dimTime));

            scalar rDeltaT = 1.0 / mesh().time().deltaT().value();

            fvScalarMatrix &ddtS = ddtSTMP.ref();

            ddtS.diag() = rDeltaT * (P_().internalField() - Q_().internalField()) * mesh().V();
            ddtS.source() = rDeltaT * (S.oldTime() - f1_() + f2_()) * mesh().V();
            return ddtSTMP;
        }

    void Casulli::update_pStar(const volScalarField &new_pStar)
    {
        pStar_.primitiveFieldRef() = new_pStar;
        pStar_.correctBoundaryConditions();
        C_pStar_.primitiveFieldRef() = poroHydraulic().C(pStar_);
        C_pStar_.correctBoundaryConditions();
        S_pStar_.primitiveFieldRef() = poroHydraulic().S(pStar_);
        S_pStar_.correctBoundaryConditions();
    };

  void Casulli::P(const volScalarField &pField, const volScalarField &C)
  {
  	P_.clear();
  	P_.reset(
  	new volScalarField(
                IOobject(
                    "PCas",
                    mesh().time().timeName(),
                    pField.db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                mesh(),
                dimensionedScalar("",C_pStar_.dimensions(),0.0),
                "zeroGradient")
                );
       volScalarField& PRef = P_.ref();

    /*forAll(PRef.boundaryField(), patchI)
    {
      forAll(PRef.boundaryField()[patchI], faceI)
      {

        if (pField.boundaryField()[patchI][faceI] <= pStar_.boundaryField()[patchI][faceI])
        {
          PRef.boundaryFieldRef()[patchI][faceI] = C.boundaryField()[patchI][faceI];
        }
        else
        {
          PRef.boundaryFieldRef()[patchI][faceI] = C_pStar_.boundaryField()[patchI][faceI];
        }
      }
    }*/

    forAll(PRef, cellI)
    {
      if (pField[cellI] <= pStar_[cellI])
      {
        PRef[cellI] = C[cellI];
      }
      else
      {
        PRef[cellI] = C_pStar_[cellI];
      }
    }
    
    PRef.correctBoundaryConditions();

  };

  void Casulli::Q(const volScalarField &pField, const volScalarField &PRef)
  {
  	Q_.clear();
  	Q_.reset(
  	new volScalarField(
                IOobject(
                    "QCas",
                    mesh().time().timeName(),
                    pField.db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                mesh(),
                dimensionedScalar("",C_pStar_.dimensions(),0.0),
                "zeroGradient")
                );
       volScalarField& QRef = Q_.ref();
           
    /*forAll(QRef.boundaryField(), patchI)
    {
      forAll(QRef.boundaryField()[patchI], faceI)
      {
        if (pField.boundaryField()[patchI][faceI] > pStar_.boundaryField()[patchI][faceI])
        {
          QRef.boundaryFieldRef()[patchI][faceI] = PRef.boundaryField()[patchI][faceI] - CRef.boundaryField()[patchI][faceI];
        }
      }
    }*/

    forAll(QRef, cellI)
    {
      if (pField[cellI] > pStar_[cellI])
      {
        QRef[cellI] = PRef[cellI] - C()[cellI];
      }
    }
    
    QRef.correctBoundaryConditions();

  };

  void Casulli::S1(const volScalarField &pField, const volScalarField &S)
  {
    	S1_.clear();
  	S1_.reset(
  	new volScalarField(
                IOobject(
                    "S1",
                    mesh().time().timeName(),
                    pField.db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                mesh(),
                dimensionedScalar("",S_pStar_.dimensions(),0.0),
                "zeroGradient")
                );
       volScalarField& S1Ref = S1_.ref();

    /*forAll(S1Ref.boundaryField(), patchI)
    {
      forAll(S1Ref.boundaryField()[patchI], faceI)
      {
        if (pField.boundaryField()[patchI][faceI] <= pStar_.boundaryField()[patchI][faceI])
        {
          S1Ref.boundaryFieldRef()[patchI][faceI] = S.boundaryField()[patchI][faceI];
        }
        else
        {
          S1Ref.boundaryFieldRef()[patchI][faceI] = S_pStar_.boundaryField()[patchI][faceI] + C_pStar_.boundaryField()[patchI][faceI] * (pField.boundaryField()[patchI][faceI] - pStar_.boundaryField()[patchI][faceI]);
        }
      }
    }*/


    forAll(S1Ref, cellI)
    {
      if (pField[cellI] <= pStar_[cellI])
      {
        S1Ref[cellI] = S[cellI];
      }
      else
      {
        S1Ref[cellI] = S_pStar_[cellI] + C_pStar_[cellI] * (pField[cellI] - pStar_[cellI]);
      }
    }
    
    S1Ref.correctBoundaryConditions();
  };

  void Casulli::S2(const volScalarField &pField, const volScalarField &S1Ref)
  {
  
      	S2_.clear();
  	S2_.reset(
  	new volScalarField(
                IOobject(
                    "S2",
                    mesh().time().timeName(),
                    pField.db(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                mesh(),
                dimensionedScalar("",S_pStar_.dimensions(),0.0),
                "zeroGradient")
                );
       volScalarField& S2Ref = S2_.ref();
       

    /*forAll(S2Ref.boundaryField(), patchI)
    {
      forAll(S2Ref.boundaryField()[patchI], faceI)
      {
        if (pField.boundaryField()[patchI][faceI] > pStar_.boundaryField()[patchI][faceI])
        {
          S2Ref.boundaryFieldRef()[patchI][faceI] = S1Ref.boundaryField()[patchI][faceI] - SRef.boundaryField()[patchI][faceI];
        }
      }
    }*/

    forAll(S2Ref, cellI)
    {
      if (pField[cellI] > pStar_[cellI])
      {
        S2Ref[cellI] = S1Ref[cellI] - S()[cellI];
      }
    }

      S2Ref.correctBoundaryConditions();
  };


    } // End of namespace richardsLinearizations
} // End of namespace Foam

//*********************************************************** //
