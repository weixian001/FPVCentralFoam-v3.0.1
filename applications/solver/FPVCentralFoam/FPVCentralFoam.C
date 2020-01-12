/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "bound.H"	//added
#include "fvCFD.H"
#include "psiCombustionModel.H"		//added
//#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

	surfaceScalarField rhoZmix_pos(interpolate(rhoZmix, pos, Zmix.name()));//added
        surfaceScalarField rhoZmix_neg(interpolate(rhoZmix, neg, Zmix.name()));//added

        surfaceScalarField rhoPv_pos(interpolate(rhoPv, pos, Pv.name()));//added
        surfaceScalarField rhoPv_neg(interpolate(rhoPv, neg, Pv.name()));//added

        surfaceScalarField rhovarZ_pos(interpolate(rhovarZ, pos, varZ.name()));//added
        surfaceScalarField rhovarZ_neg(interpolate(rhovarZ, neg, varZ.name()));//added

        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());

        volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"
        #include "readTimeControls.H"

/*        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {*/
            #include "setDeltaT.H"
        //}

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );

        surfaceScalarField phiPv      //added
        (
            "phiPv",
            aphiv_pos*rhoPv_pos + aphiv_neg*rhoPv_neg
        );

        surfaceScalarField phiZmix	//added
        (
            "phiZmix",
            aphiv_pos*rhoZmix_pos + aphiv_neg*rhoZmix_neg
        );

        surfaceScalarField phivarZ      //added
        (
            "phivarZ",
            aphiv_pos*rhovarZ_pos + aphiv_neg*rhovarZ_neg
        );

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() == rho.boundaryField()*U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

	// --- Solve Pv or Chi (Progress Variable)
	solve(fvm::ddt(rhoPv) + fvc::div(phiPv)); // add source here or in diff-cor part?

        Pv.internalField() =
            rhoPv.dimensionedInternalField()
           /rho.dimensionedInternalField();
        Pv.correctBoundaryConditions();
        rhoPv.boundaryField() == rho.boundaryField()*Pv.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, Pv) - fvc::ddt(rho, Pv)
              - fvm::laplacian(turbulence->DZEff(), Pv)
	      - thermo.Srr() // add here or inviscid part?
            );
	    Info<< "----------> Pv min/max   = " << min(Pv).value() << ", "
            << max(Pv).value() << endl;
            rhoPv = rho*Pv;
        }


	// --- Solve Zmix (Mixture Fraction)
        solve(fvm::ddt(rhoZmix) + fvc::div(phiZmix));

        //Zmix.dimensionedInternalField() =
        Zmix.internalField() =
            rhoZmix.dimensionedInternalField()
           /rho.dimensionedInternalField();
        Zmix.correctBoundaryConditions();
        rhoZmix.boundaryField() == rho.boundaryField()*Zmix.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, Zmix) - fvc::ddt(rho, Zmix)
              - fvm::laplacian(turbulence->DZEff(), Zmix)
            );
            //bound(Zmix, 0.0);
	    Info<< "----------> Zmix min/max   = " << min(Zmix).value() << ", "
            << max(Zmix).value() << endl;
            rhoZmix = rho*Zmix;
        }

	// --- Solve varZ (Mixture Fraction Variance)

        solve(fvm::ddt(rhovarZ) + fvc::div(phivarZ));

        varZ.internalField() =
            rhovarZ.dimensionedInternalField()
           /rho.dimensionedInternalField();
        varZ.correctBoundaryConditions();
        rhovarZ.boundaryField() == rho.boundaryField()*varZ.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, varZ) - fvc::ddt(rho, varZ)
              - fvm::laplacian(turbulence->DZEff(), varZ)
	      - 2.0*turbulence->Production()*magSqr(fvc::grad(Zmix))
	       + turbulence->Destruction()*varZ
            );
	    bound(varZ, 0.0);
	    Info<< "----------> varZ min/max   = " << min(varZ).value() << ", "
            << max(varZ).value() << endl;
            rhovarZ = rho*varZ;
        }

	// --- Update the Srr() and dq() 
        combustion->correct();

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + (mesh.Sf() & fvc::interpolate(tauMC))
            )
            & (a_pos*U_pos + a_neg*U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );

        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryField() ==
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence->alphaEff(), e)
	      - thermo.dq()
            );
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
        }
        
	p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() == psi.boundaryField()*p.boundaryField();

        //combustion->correct();
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
