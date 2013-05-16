#include <iostream>
#include <cmath>
#include <openbabel/obconversion.h>
#include <openbabel/conformersearch.h>
#include <openbabel/mol.h>
#include <openbabel/data.h>

#ifdef __APPLE__
    #include <Accelerate/Accelerate.h>
    #define FORTRANINT int
#else
    #include <f2c.h>
    #include <blaswrap.h>
    #include <cblas.h>
    #include <clapack.h>
    #define FORTRANINT long int
#endif

using namespace std;
using namespace OpenBabel;

#include "Utils.cpp"

double volumeOverlap(const vector<double> &coordsMoleculeA, const vector<double> &coordsMoleculeB, const vector<int> &atomicNumbersA, const vector<int> &atomicNumbersB) {
    if (coordsMoleculeA.size() != atomicNumbersA.size() * 3 or coordsMoleculeB.size() != atomicNumbersB.size() * 3) { cerr << "ERROR: INCORRECT MATCHING OF NUMBER OF COORDINATES AND ATOMIC NUMBERS; EXITING" << endl; abort(); }

    double totalVolumeOverlap = 0;
    const double constP = 2.0 * M_SQRT2;
    const double A = 4.0 * M_PI * constP / 3.0;
    const double B = -M_PI * pow(0.75 * constP * M_1_PI, 2.0/3.0);

    for (unsigned int i=0; i < coordsMoleculeA.size(); i+=3) {
        double vdwRA = etab.GetVdwRad(atomicNumbersA[i / 3]); // atomic number vectir is 3x shorter than 3d coordinates vector

        for (unsigned int j=0; j < coordsMoleculeB.size(); j+=3) {
            double vdwRB = etab.GetVdwRad(atomicNumbersB[j / 3]);
            
            double sqvA = vdwRA * vdwRA;
            double sqvB = vdwRB * vdwRB;
            double C = sqvA + sqvB;

            double distanceSquared = pow(coordsMoleculeB[j]-coordsMoleculeA[i], 2) + pow(coordsMoleculeB[j+1]-coordsMoleculeA[i+1], 2) + pow(coordsMoleculeB[j+2]-coordsMoleculeA[i+2], 2);
            totalVolumeOverlap += A * pow(sqvA * sqvB  / C, 1.5) * exp(B * distanceSquared / C);
        }
    }
    return totalVolumeOverlap;
}

double volumeOverlap (OBMol &moleculeA, OBMol &moleculeB) {
    vector<double> coordsA, coordsB;
    vector<int> atomNumsA, atomNumsB;
    generateCoordsMatrixFromMolecule(coordsA, moleculeA);
    generateCoordsMatrixFromMolecule(coordsB, moleculeB);
    generateAtomicNumbersListFromMolecule(atomNumsA, moleculeA);
    generateAtomicNumbersListFromMolecule(atomNumsB, moleculeB);
    return volumeOverlap(coordsA, coordsB, atomNumsA, atomNumsB); 
}

void generateEigenMatrix(vector<double> &eigenvectors, vector<double> &eigenvalues, const vector<double> &matrix) {
    // http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen.html#gaeed8a131adf56eaa2a9e5b1e0cce5718  |  http://www.ualberta.ca/~kbeach/lapack.html
    FORTRANINT matrixOrder = 3;
    eigenvalues.clear(); eigenvalues.resize(3, 0);
    eigenvectors.clear(); eigenvectors = matrix;
    FORTRANINT lwork = 102, info; vector<double> work(lwork); // LAPACK parameters; work[0] will show the optimum value for lwork when call completes
    int lapackErrCode = dsyev_("V", "L", &matrixOrder, &eigenvectors[0], &matrixOrder, &eigenvalues[0], &work[0], &lwork, &info);
    if (lapackErrCode) { cerr << "WARNING: LAPACK DSYEV() RETURNED NON-ZERO ERROR CODE " << lapackErrCode << endl; }
}

void generateCovarMatrixFromMolecule(vector<double> &matrix, OBMol &molecule) {
    double uX = 0, uY = 0, uZ = 0, *moleculeCoords = molecule.GetCoordinates();
    for (unsigned int i=0; i < 3 * molecule.NumAtoms(); i+=3) { uX += moleculeCoords[i]; uY += moleculeCoords[i+1]; uZ += moleculeCoords[i+2]; }
    uX /= molecule.NumAtoms(); uY /= molecule.NumAtoms(); uZ /= molecule.NumAtoms();

    double cXX = 0, cYY = 0, cZZ = 0, cXY = 0, cXZ = 0, cYZ = 0;
    for (unsigned int i=0; i < 3 * molecule.NumAtoms(); i+=3) {
        double atomWeight = molecule.GetAtom(i/3 + 1)->GetAtomicMass(); // getAtoms is 1-based instead of 0-based
        cXX += pow(moleculeCoords[i] - uX, 2) * atomWeight;
        cYY += pow(moleculeCoords[i+1] - uY, 2) * atomWeight;
        cZZ += pow(moleculeCoords[i+2] - uZ, 2) * atomWeight;
        cXY += (moleculeCoords[i] - uX) * (moleculeCoords[i+1] - uY) * atomWeight;
        cXZ += (moleculeCoords[i] - uX) * (moleculeCoords[i+2] - uZ) * atomWeight;
        cYZ += (moleculeCoords[i+1] - uY) * (moleculeCoords[i+2] - uZ) * atomWeight;
    }

    matrix.clear(); matrix.resize(9);
    matrix[0] = cXX; 
    matrix[1] = matrix[3] = cXY; 
    matrix[2] = matrix[6] = cXZ; 
    matrix[4] = cYY;
    matrix[5] = matrix[7] = cYZ;  
    matrix[8] = cZZ;

    for (unsigned int i=0; i < matrix.size(); i++) matrix[i] /= molecule.GetMolWt();
}

void generatePCARotationMatrix(vector<double> &rotMatrix, unsigned int optCode, const vector<double> &U, const vector<double> &V) {
    if (U.size() != 9 or V.size() != 9 or optCode > 3) { cerr << "ERROR: SIZE OF U OR V IS INCORRECT, OR OPTCODE IS WRONG; EXITING" << endl; abort(); }
    rotMatrix.clear(); rotMatrix.resize(9);

    if (optCode == 0) { 
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, 3, 1, &V[0], 3, &U[0], 3, 0, &rotMatrix[0], 3);
    } else {
        vector<double> tmpU(U);
        if (optCode == 3) {
            for (int i=0; i<6; i++) tmpU[i] = -tmpU[i];
        } else if (optCode == 2) {
            for (int i=0; i<3; i++) tmpU[i] = -tmpU[i];
            for (int i=6; i<tmpU.size(); i++) tmpU[i] = -tmpU[i];
        } else if (optCode == 1) {
            for (int i=3; i<tmpU.size(); i++) tmpU[i] = -tmpU[i];
        }  
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, 3, 1, &V[0], 3, &tmpU[0], 3, 0, &rotMatrix[0], 3);
    }
}

void findBestPCAOrientation(OBMol &moleculeA, OBMol &moleculeB) {
    cout << endl << "BEGIN INITIAL ORIENTATION SEARCH" << endl
        << "Searching for the best initial orientation matrix..." << endl;
    vector<double> coordA, coordB, comA, comB, covA, covB, eVectA, eVectB, eValA, eValB, tempR, bestR, bestA;
    vector<int> atomicNumsA, atomicNumsB;
    map<int, string> RTable; RTable[0] = "R0"; RTable[1] = "Rx"; RTable[2] = "Ry"; RTable[3] = "Rz";

    generateCoordsMatrixFromMolecule(coordA, moleculeA);
    generateCoordsMatrixFromMolecule(coordB, moleculeB);
    cout << "\nORIGINAL COORDS OF MOLECULE A:" << endl; printMoleculeCoords(moleculeA);
    cout << "\nORIGINAL COORDS OF MOLECULE B:" << endl; printMoleculeCoords(moleculeB);

    generateAtomicNumbersListFromMolecule(atomicNumsA, moleculeA);
    generateAtomicNumbersListFromMolecule(atomicNumsB, moleculeB);

    getMoleculeCenterCoords(comA, moleculeA);
    getMoleculeCenterCoords(comB, moleculeB);

    generateCovarMatrixFromMolecule(covA, moleculeA);
    generateCovarMatrixFromMolecule(covB, moleculeB);

    generateEigenMatrix(eVectA, eValA, covA);
    generateEigenMatrix(eVectB, eValB, covB);

    double bestVolumeOverlap=0; int bestRcode=0;
    for (int i=0; i<4; i++) {
        generatePCARotationMatrix(tempR, i, eVectA, eVectB);
        vector<double> currentPCACoordA = coordA;

        translate3DMatrixCoordinates(currentPCACoordA, -comA[0], -comA[1], -comA[2]);
        rotate3DMatrixCoordinates(currentPCACoordA, tempR);
        translate3DMatrixCoordinates(currentPCACoordA, comB[0], comB[1], comB[2]);

        double curVolOverlap = volumeOverlap(currentPCACoordA, coordB, atomicNumsA, atomicNumsB);
        if (curVolOverlap > bestVolumeOverlap) {
            bestRcode = i;
            bestVolumeOverlap = curVolOverlap;
            bestR = tempR;
            bestA = currentPCACoordA;
        }
    }

    cout << "\nThe best initial orientation matrix is: " << RTable[bestRcode] << ", which produces a volume overlap of " << bestVolumeOverlap << endl;
    printMatrix(bestR, 3, 3);
    cout << "\nRESULTING A:" << endl; printMatrix(bestA, bestA.size() / 3, 3, false);
    cout << "\nEND INITIAL ORIENTATION SEARCH.  SAVING COORDINATES TO MOLECULE A...\n\n";
    saveCoordsMatrixToMolecule(moleculeA, bestA);    

}






#include "SteepestDescent.cpp"



void PCAPlusSteepestDescent(OBMol &moleculeA, OBMol &moleculeB, double alpha, double betaRadians) {
    cout << endl << "BEGIN ORIENTATION SEARCH" << endl;
    vector<double> coordA, coordB, comA, comB, covA, covB, eVectA, eVectB, eValA, eValB, tempR, bestR, bestA;
    vector<int> atomicNumsA, atomicNumsB;
    map<int, string> RTable; RTable[0] = "R0"; RTable[1] = "Rx"; RTable[2] = "Ry"; RTable[3] = "Rz";


    generateCoordsMatrixFromMolecule(coordA, moleculeA); cout << "\nORIGINAL COORDS OF MOLECULE A:" << endl; printMoleculeCoords(moleculeA);
    generateCoordsMatrixFromMolecule(coordB, moleculeB); cout << "\nORIGINAL COORDS OF MOLECULE B:" << endl; printMoleculeCoords(moleculeB);
    generateAtomicNumbersListFromMolecule(atomicNumsA, moleculeA);
    generateAtomicNumbersListFromMolecule(atomicNumsB, moleculeB);
    getMoleculeCenterCoords(comA, moleculeA);
    getMoleculeCenterCoords(comB, moleculeB);
    
    generateCovarMatrixFromMolecule(covA, moleculeA);
    generateCovarMatrixFromMolecule(covB, moleculeB);
    generateEigenMatrix(eVectA, eValA, covA);
    generateEigenMatrix(eVectB, eValB, covB);

    double bestVolumeOverlap=0; int bestRcode=0;
    for (int i=0; i<4; i++) {
        generatePCARotationMatrix(tempR, i, eVectA, eVectB);
        vector<double> currentPCACoordA = coordA, currentSDCoordA;

        translate3DMatrixCoordinates(currentPCACoordA, -comA[0], -comA[1], -comA[2]);
        rotate3DMatrixCoordinates(currentPCACoordA, tempR);
        translate3DMatrixCoordinates(currentPCACoordA, comB[0], comB[1], comB[2]);

        cout << endl << "BEGIN STEEPEST DESCENT SEARCH..." << endl;
        double curVolOverlap = steepestDescentEngine(currentSDCoordA, currentPCACoordA, coordB, atomicNumsA, atomicNumsB, comA, alpha, betaRadians);
        
        cout << "END STEEOEST DESCENT SEARCH.\nThe convergent solution starting with position PCA0 produces a volume overlap of " << curVolOverlap << endl;
        if (curVolOverlap > bestVolumeOverlap) {
            bestRcode = i;
            bestVolumeOverlap = curVolOverlap;
            bestR = tempR;
            bestA = currentSDCoordA;
        }
    }

    cout << "\nThe best initial orientation matrix is: " << RTable[bestRcode] << ", which, after running Steepest Descent, produces a volume overlap of " << bestVolumeOverlap << endl;
    cout << "PCA matrix " << RTable[bestRcode] << ":" << endl; printMatrix(bestR, 3, 3);
    cout << "\nResulting A after Steepest Descent:" << endl; printMatrix(bestA, bestA.size() / 3, 3, false);
    cout << "\nEND ORIENTATION SEARCH.  SAVING COORDINATES TO MOLECULE A...\n\n";
    saveCoordsMatrixToMolecule(moleculeA, bestA);
}


void generateConformers(OBMol &molecule, int numConformers=50, int numChildren=10, int mutability=5, int convergence=30) {
    OBConformerSearch cs;
    cs.Setup(molecule, numConformers, numChildren, mutability, convergence); // numConformers 30 // numChildren 5 // mutability 5 // convergence 25
    cs.SetScore(new OBEnergyConformerScore);
    cs.Search(); cs.GetConformers(molecule);
}


int main (int argc, char **argv) {
    if(argc < 3) {
        cout << "Usage: ProgrameName InputFileName OutputFileName\n";
        return 1;
    }

    vector<OBMol> molecules;
    importMoleculesFromFile(molecules, argv[1]);
    generateConformers(molecules[0]);
    writeMoleculeConformersToFile(argv[2], molecules[0], true);


/*    cout << "VEC SIZE: " << molecules.size() << endl;
    cout << "NUM CONFORMERS: " << molecules[0].NumConformers() << endl;
    cout << "NUM ATOMS: " << molecules[0].NumAtoms() << endl;
*/
    //generateConformers(molecules[1]);
    //molecules[1].SetConformer(10);


    //findBestPCAOrientation(molecules[0], molecules[1]);
    //runSteepestDescent(molecules[0], molecules[1], 0.5, 10.0 * M_PI / 180.0);
    
    //PCAPlusSteepestDescent(molecules[0], molecules[1], 0.5, 10.0 * M_PI / 180.0);
    //writeMoleculeToFile(argv[3], molecules[0], true);
    //writeMoleculeToFile(argv[3], molecules[1]);

    //sampleTest(mol);
    //testRot();
    //testEigen();

    return 0;
}
