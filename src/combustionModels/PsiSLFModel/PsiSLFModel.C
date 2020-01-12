/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

Contributors/Copyright
    2014 Hagen Müller <hagen.mueller@unibw.de> Universität der Bundeswehr München
    2014 Likun Ma <L.Ma@tudelft.nl> TU Delft
    2019 Lim Wei Xian <weixian001@e.ntu.edu.sg> NTUsg

\*---------------------------------------------------------------------------*/

#include "PsiSLFModel.H"
#include "reactingMixture.H"
#include "volFields.H"
#include "hashedWordList.H"
//#include "physicoChemicalConstants.H"

namespace Foam
{
namespace combustionModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType>
PsiSLFModel<CombThermoType>::PsiSLFModel
(
    const word& modelType, const fvMesh& mesh, const word& phaseName
)
:
    CombThermoType(modelType, mesh, phaseName),
    solver_(tableSolver(mesh, tables(), Lambda_table())),	//
    Y_(this->thermo().composition().Y()),
    dq_(this->thermo().dq()),
    Srr_(this->thermo().Srr()),			//added             
    Z_(this->thermo().Z()),
    Pv_(this->thermo().Pv()),				//added
    varZ_(this->thermo().varZ()),
    //Chi_(this->thermo().Chi()),
    ubIF_(mesh.cells().size()),
    ubP_(),
    posIF_(mesh.cells().size()),
    posP_(),
    //useScalarDissipation_(this->coeffs().lookup("useScalarDissipation")),
    useProgressVariable_(this->coeffs().lookup("useProgressVariable")),		//added
    useMixtureFractionVariance_(this->coeffs().lookup("useMixtureFractionVariance"))
{
	const polyBoundaryMesh& patches = mesh.boundaryMesh();
	int patchSize = 0;
    forAll(patches, patchI)
    {
    	const polyPatch& pp = patches[patchI];
    	if (pp.size() > patchSize) patchSize = pp.size();
    }

    ubP_.setSize(patchSize);
    posP_.setSize(patchSize);
}

// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType>
PsiSLFModel<CombThermoType>::~PsiSLFModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CombThermoType>
hashedWordList PsiSLFModel<CombThermoType>::tables()
{
	hashedWordList tableNames = this->thermo().composition().species();
	//tableNames.append("chi");
	tableNames.append("dq");
	tableNames.append("Srr");		//added
	
	return tableNames;
}

template<class CombThermoType>
word PsiSLFModel<CombThermoType>::Lambda_table()
{
        word LambdaNames = "Pv_lambda_table";
        return LambdaNames;
}


template<class CombThermoType>
void PsiSLFModel<CombThermoType>::correct()
{
    
    // bound the progress variable
    scalar LambdaMax = solver_.maxPv();

    const scalarField& ZCells = Z_.internalField();
    const scalarField& varZCells = varZ_.internalField();
    const scalarField& PvCells = Pv_.internalField();

    //scalarField& chiCells = Chi_.internalField();
    scalarField& dqCells = dq_.internalField();
    scalarField& SrrCells = Srr_.internalField();

    //- Update the species and enthalpy field and reaction rate
//    if(this->active())
    {
       scalarList x(3, 0.0);
       scalarList y(3, 0.0);
       double Zeta;
       double Lambda;

       // Interpolate for internal Field 
       forAll(Y_, i)
       {
    	  scalarField& YCells = Y_[i].internalField();

          forAll(ZCells, cellI)
          {
        	 if (i == 0)
        	 {
        	 Zeta = sqrt(varZCells[cellI]/max(ZCells[cellI]*(1 - ZCells[cellI]), SMALL));
                 y[0] = min(Zeta, 0.99);
                 y[1] = ZCells[cellI];
                 y[2] = max(0.0, PvCells[cellI]);

                 List<int> ubCL_ = solver_.upperBounds_Lambda(y);
                 scalarList posCL_ = solver_.position_Lambda(ubCL_, y);

                 scalarList Lambda_ = solver_.interpolateS(ubCL_, posCL_);
                 Lambda = solver_.upperBoundsLambda_(y[2],Lambda_);

                 //if (useProgressVariable_)   x[0] = min(chiMax, max(0,chiCells[cellI]));
                 if (useProgressVariable_)   x[0] = min(Lambda, LambdaMax);
                 if (useMixtureFractionVariance_) x[1] = min(Zeta, 0.99);
                 x[2] = ZCells[cellI];	

                 ubIF_[cellI] = solver_.upperBounds(x);
                 posIF_[cellI] = solver_.position(ubIF_[cellI], x);
                 
            	 dqCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 2));
            	 SrrCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], (solver_.sizeTableNames() - 1));		//added
        	 }

        	 YCells[cellI] = solver_.interpolate(ubIF_[cellI], posIF_[cellI], i);
          }
       }

       // Interpolate for patches
       forAll(dq_.boundaryField(), patchi)   
       {
          const fvPatchScalarField& pPv = Pv_.boundaryField()[patchi];
          const fvPatchScalarField& pvarZ = varZ_.boundaryField()[patchi];
          const fvPatchScalarField& pZ = Z_.boundaryField()[patchi];

          fvPatchScalarField& pdq = dq_.boundaryField()[patchi];
          fvPatchScalarField& pSrr = Srr_.boundaryField()[patchi];	//added

          forAll(Y_, i)
          {
        	  fvPatchScalarField& pY = Y_[i].boundaryField()[patchi];

              forAll(pY , facei)
              {
             	 if (i == 0)
             	 {
                     Zeta = sqrt(pvarZ[facei]/max(pZ[facei]*(1 - pZ[facei]), SMALL));
                     y[0] = min(Zeta, 0.99);
                     y[1] = pZ[facei];
                     y[2] = max(0.0, pPv[facei]);

                     List<int> ubCL_ = solver_.upperBounds_Lambda(y);
                     scalarList posCL_ = solver_.position_Lambda(ubCL_, y);

                     scalarList Lambda_ = solver_.interpolateS(ubCL_, posCL_);
                     Lambda = solver_.upperBoundsLambda_(y[2],Lambda_);

                     //if (useProgressVariable_)   x[0] = min(chiMax, max(0, pChi[facei]));
                     if (useProgressVariable_)   x[0] = min(Lambda, LambdaMax);
                     if (useMixtureFractionVariance_) x[1] = min(Zeta, 0.99);
                     x[2] = pZ[facei];

                     ubP_[facei] = solver_.upperBounds(x);
                     posP_[facei] = solver_.position(ubP_[facei], x);

                     //pChi[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 3));
                     pdq[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 2));
                     pSrr[facei] = solver_.interpolate(ubP_[facei], posP_[facei], (solver_.sizeTableNames() - 1));			//added
             	 }

            	 pY[facei] = solver_.interpolate(ubP_[facei], posP_[facei], i);
             }
          }
       }
       //this->thermo().correct();
    Info << "----------> Srr min/max =" << min(Srr_).value() << ", "
        << max(Srr_).value() << endl;
    Info << "----------> dq min/max =" << min(dq_).value() << ", "
        << max(dq_).value() << endl;   
    }
}

template<class CombThermoType>
Switch PsiSLFModel<CombThermoType>::correctDensity()
{
	return true;
}

template<class CombThermoType>
Foam::tmp<Foam::fvScalarMatrix>
PsiSLFModel<CombThermoType>::R
(
    volScalarField& Y              
) const
{
    tmp<fvScalarMatrix> tSu(new fvScalarMatrix(Y, dimMass/dimTime));
    return tSu;
}

template<class CombThermoType>
Foam::tmp<Foam::volScalarField>
PsiSLFModel< CombThermoType>::Sh() const
{
    tmp<volScalarField> tSh
    (
        new volScalarField
        (
            IOobject
            (
                "Sh",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tSh;
}

template<class CombThermoType>
Foam::tmp<Foam::volScalarField>
PsiSLFModel< CombThermoType>::dQ() const
{
    tmp<volScalarField> tdQ
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("dQ", dimEnergy/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tdQ;
}

template<class CombThermoType>
bool PsiSLFModel<CombThermoType>::read()
{
    if (CombThermoType::read())
    {
        //this->coeffs().lookup("useScalarDissipation") >> useScalarDissipation_;
        this->coeffs().lookup("useProgressVariable") >> useProgressVariable_;
        this->coeffs().lookup("useMixtureFractionVariance") >> useMixtureFractionVariance_;
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
