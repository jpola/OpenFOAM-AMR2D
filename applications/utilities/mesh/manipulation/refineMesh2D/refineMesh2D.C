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
    refineMesh2D

Description
    Adaptive mesh refinement for 2D mesh

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "Time.H"
#include "dynamicFvMesh.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "fvCFD.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    //I put the from #includes explicitly here for better code management.

    //#include "setRootCase.H"
    Foam::argList args(argc, argv);
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    //#include "createTime.H"
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);

    runTime.functionObjects().off();

    //#include "createDynamicFvMesh.H"
    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    autoPtr<dynamicFvMesh> meshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                dynamicFvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
    Info << "\nLook here in the tutorial case RubiksCube/constant I defined dynamicMeshDict\n"
         << "There in line #9. dynamicFvMesh dynamicRefineFvMesh2D caused that\n"
         << "our mes is of type dynamicRefineFvMesh2D. Briliant! :)\n" << endl;

    dynamicFvMesh& mesh = meshPtr();

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Create phi
    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    );

    phi = linearInterpolate(U) & mesh.Sf();

    //Create alpha field
    Info<< "Reading transportProperties\n" << endl;
    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2
    );
    rho.oldTime();

    runTime++;

    mesh.update();

    runTime.write();

	Info << "Refine Mesh 2D End\n" << endl;
	return 0;

}
