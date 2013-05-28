#include <iostream>
#include <cmath>
#include <openbabel/obconversion.h>
#include <openbabel/conformersearch.h>
#include <openbabel/forcefield.h>
#include <openbabel/op.h>
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

double volumeOverlap(const vector<double> &coordsMoleculeA, const vector<double> &coordsMoleculeB, const vector<double> &VDWsA, const vector<double> &VDWsB) {
    if (coordsMoleculeA.size() != VDWsA.size() * 3 or coordsMoleculeB.size() != VDWsB.size() * 3) { cerr << "ERROR: INCORRECT MATCHING OF NUMBER OF COORDINATES AND ATOMIC NUMBERS; EXITING" << endl; abort(); }

    double totalVolumeOverlap = 0;
    static const double constP = 2.0 * M_SQRT2;
    static const double A = 4.0 * M_PI * constP / 3.0;
    static const double B = -M_PI * pow(0.75 * constP * M_1_PI, 2.0/3.0);

    for (unsigned int i=0; i < coordsMoleculeA.size(); i+=3) {
        double vdwRA = VDWsA[i / 3]; // atomic number vectir is 3x shorter than 3d coordinates vector
        double sqvA = vdwRA * vdwRA;

        for (unsigned int j=0; j < coordsMoleculeB.size(); j+=3) {
            double vdwRB = VDWsB[j / 3];            
            double sqvB = vdwRB * vdwRB;
            double C = sqvA + sqvB;

            double distanceSquared = pow(coordsMoleculeB[j]-coordsMoleculeA[i], 2) + pow(coordsMoleculeB[j+1]-coordsMoleculeA[i+1], 2) + pow(coordsMoleculeB[j+2]-coordsMoleculeA[i+2], 2);
            totalVolumeOverlap += A * pow(sqvA * sqvB  / C, 1.5) * exp(B * distanceSquared / C);
        }
    }
    return totalVolumeOverlap;
}

double similarityIndex(const vector<double> &coordsMoleculeA, const vector<double> &coordsMoleculeB, const vector<double> &VDWsA, const vector<double> &VDWsB) {
    double overlapAA = volumeOverlap(coordsMoleculeA, coordsMoleculeA, VDWsA, VDWsA);
    double overlapBB = volumeOverlap(coordsMoleculeB, coordsMoleculeB, VDWsB, VDWsB);
    double overlapAB = volumeOverlap(coordsMoleculeA, coordsMoleculeB, VDWsA, VDWsB);
    return 2.0 * overlapAB / (overlapAA + overlapBB);
}

double volumeOverlap(OBMol &moleculeA, OBMol &moleculeB) {
    vector<double> coordsA, coordsB, VDWsA, VDWsB;
    generateCoordsMatrixFromMolecule(coordsA, moleculeA);
    generateCoordsMatrixFromMolecule(coordsB, moleculeB);
    generateVDWRadiusListFromMolecule(VDWsA, moleculeA);
    generateVDWRadiusListFromMolecule(VDWsB, moleculeB);
    return volumeOverlap(coordsA, coordsB, VDWsA, VDWsB); 
}

double similarityIndex(OBMol &moleculeA, OBMol &moleculeB) {
    vector<double> coordsA, coordsB, VDWsA, VDWsB;
    generateCoordsMatrixFromMolecule(coordsA, moleculeA);
    generateCoordsMatrixFromMolecule(coordsB, moleculeB);
    generateVDWRadiusListFromMolecule(VDWsA, moleculeA);
    generateVDWRadiusListFromMolecule(VDWsB, moleculeB);
    return similarityIndex(coordsA, coordsB, VDWsA, VDWsB); 
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
    vector<double> coordA, coordB, comA, comB, covA, covB, eVectA, eVectB, eValA, eValB, tempR, bestR, bestA, VDWsA, VDWsB;
    map<int, string> RTable; RTable[0] = "R0"; RTable[1] = "Rx"; RTable[2] = "Ry"; RTable[3] = "Rz";

    generateCoordsMatrixFromMolecule(coordA, moleculeA);
    generateCoordsMatrixFromMolecule(coordB, moleculeB);
    cout << "\nORIGINAL COORDS OF MOLECULE A:" << endl; printMoleculeCoords(moleculeA);
    cout << "\nORIGINAL COORDS OF MOLECULE B:" << endl; printMoleculeCoords(moleculeB);

    generateVDWRadiusListFromMolecule(VDWsA, moleculeA);
    generateVDWRadiusListFromMolecule(VDWsB, moleculeB);

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

        double curVolOverlap = volumeOverlap(currentPCACoordA, coordB, VDWsA, VDWsB);
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



double PCAPlusSteepestDescent(OBMol &moleculeA, OBMol &moleculeB, double alpha, double betaRadians, bool verbose=false) {
    if (verbose) cout << endl << "BEGIN ORIENTATION SEARCH" << endl;
    vector<double> coordA, coordB, comA, comB, covA, covB, eVectA, eVectB, eValA, eValB, tempR, bestR, bestA, VDWsA, VDWsB;
    map<int, string> RTable; RTable[0] = "R0"; RTable[1] = "Rx"; RTable[2] = "Ry"; RTable[3] = "Rz";


    generateCoordsMatrixFromMolecule(coordA, moleculeA); if (verbose) { cout << "\nORIGINAL COORDS OF MOLECULE A:" << endl; printMoleculeCoords(moleculeA); }
    generateCoordsMatrixFromMolecule(coordB, moleculeB); if (verbose) { cout << "\nORIGINAL COORDS OF MOLECULE B:" << endl; printMoleculeCoords(moleculeB); }
    generateVDWRadiusListFromMolecule(VDWsA, moleculeA);
    generateVDWRadiusListFromMolecule(VDWsB, moleculeB);
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

        if (verbose) cout << endl << "BEGIN STEEPEST DESCENT SEARCH..." << endl;
        double curVolOverlap = steepestDescentEngine(currentSDCoordA, currentPCACoordA, coordB, VDWsA, VDWsB, comA, alpha, betaRadians, verbose);
        
        if (verbose) cout << "END STEEOEST DESCENT SEARCH.\nThe convergent solution starting with position PCA0 produces a volume overlap of " << curVolOverlap << endl;
        if (curVolOverlap > bestVolumeOverlap) {
            bestRcode = i;
            bestVolumeOverlap = curVolOverlap;
            bestR = tempR;
            bestA = currentSDCoordA;
        }
    }

    if (verbose) {
        cout << "\nThe best initial orientation matrix is: " << RTable[bestRcode] << ", which, after running Steepest Descent, produces a volume overlap of " << bestVolumeOverlap << endl;
        cout << "PCA matrix " << RTable[bestRcode] << ":" << endl; printMatrix(bestR, 3, 3);
        cout << "\nResulting A after Steepest Descent:" << endl; printMatrix(bestA, bestA.size() / 3, 3, false);
        cout << "\nEND ORIENTATION SEARCH.  SAVING COORDINATES TO MOLECULE A...\n\n";
    }
    saveCoordsMatrixToMolecule(moleculeA, bestA);
    return bestVolumeOverlap;
}

void PCAEngine(vector<double> &finalCoordsA, vector<double> &coordsMoleculeA, vector<double> &coordsMoleculeB, vector<double> &eVectA, vector<double> &eVectB, vector<double> &comA, vector<double> &comB, vector<double> &VDWsA, vector<double> &VDWsB) {
    vector<double> tempR; double bestVolumeOverlap = 0;

    for (unsigned int k=0; k < 4; k++) {
        generatePCARotationMatrix(tempR, k, eVectA, eVectB);

        vector<double> currentPCACoordA = coordsMoleculeA;
        translate3DMatrixCoordinates(currentPCACoordA, -comA[0], -comA[1], -comA[2]);
        rotate3DMatrixCoordinates(currentPCACoordA, tempR);
        translate3DMatrixCoordinates(currentPCACoordA, comB[0], comB[1], comB[2]);

        double currentVolumeOverlap = volumeOverlap(currentPCACoordA, coordsMoleculeB, VDWsA, VDWsB);
        if (currentVolumeOverlap > bestVolumeOverlap) {
            finalCoordsA = currentPCACoordA;
            bestVolumeOverlap = currentVolumeOverlap;
        }
    }
}


void PCACovarianceMatrix(vector<double> &matrix, vector<double> &moleculeCoords, vector<double> &atomicMasses, double molecularWeight) {
    unsigned int numAtoms = moleculeCoords.size() / 3;
    double uX = 0, uY = 0, uZ = 0;
    for (unsigned int i=0; i < 3 * numAtoms; i+=3) { uX += moleculeCoords[i]; uY += moleculeCoords[i+1]; uZ += moleculeCoords[i+2]; }
    uX /= numAtoms; uY /= numAtoms; uZ /= numAtoms;

    double cXX = 0, cYY = 0, cZZ = 0, cXY = 0, cXZ = 0, cYZ = 0;
    for (unsigned int i=0; i < 3 * numAtoms; i+=3) {
        double atomWeight = atomicMasses[i / 3]; // getAtoms is 1-based instead of 0-based
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

    for (unsigned int i=0; i < matrix.size(); i++) matrix[i] /= molecularWeight;
}

void getPCAEigenMatrix(vector<double> &eVects, vector<double> &moleculeCoords, vector<double> &atomicMasses, double molecularWeight) {
    vector<double> covA, eValA;
    PCACovarianceMatrix(covA, moleculeCoords, atomicMasses, molecularWeight);
    generateEigenMatrix(eVects, eValA, covA);
}

void runComformerComparisons(OBMol &moleculeA, OBMol &moleculeB, bool verbose=false) { // A is the target and B is the reference
    vector< vector<double> > coordAs, coordBs, comAs, comBs, eVectAs, eVectBs;
    vector<double> VDWsA, VDWsB, massesA, massesB, currentPCACoordA, currentSDCoordA, bestCoordsA;

    double molecularWeightA = moleculeA.GetMolWt(), molecularWeightB = moleculeB.GetMolWt();
    double bestVolumeOverlap = -1;
    int bestJ = -1, bestI = -1, stepCount = 0;
    
    generateVDWRadiusListFromMolecule(VDWsA, moleculeA);
    generateVDWRadiusListFromMolecule(VDWsB, moleculeB);
    generateAtomicMassesListFromMolecule(massesA, moleculeA);
    generateAtomicMassesListFromMolecule(massesB, moleculeB);
    generateCoordsMatrixFromMoleculeConformers(coordAs, moleculeA);
    generateCoordsMatrixFromMoleculeConformers(coordBs, moleculeB);
    getMoleculeConformerCenterCoords(comAs, moleculeA);
    getMoleculeConformerCenterCoords(comBs, moleculeB);
    eVectAs.resize( moleculeA.NumConformers() ), eVectBs.resize( moleculeB.NumConformers() );
    for (unsigned int k=0; k < moleculeA.NumConformers(); k++) getPCAEigenMatrix(eVectAs[k], coordAs[k], massesA, molecularWeightA);
    for (unsigned int k=0; k < moleculeB.NumConformers(); k++) getPCAEigenMatrix(eVectBs[k], coordBs[k], massesB, molecularWeightB);
    cout << "Finished setting up data; running search...\n";

    for (unsigned int j=0; j < moleculeB.NumConformers(); j++) {
        for (unsigned int i=0; i < moleculeA.NumConformers(); i++) {
            PCAEngine(currentPCACoordA, coordAs[i], coordBs[j], eVectAs[i], eVectBs[j], comAs[i], comBs[j], VDWsA, VDWsB);
            double currentVolumeOverlap = steepestDescentEngine(currentSDCoordA, currentPCACoordA, coordBs[j], VDWsA, VDWsB, comAs[i], 1.0, 10.0 * M_PI / 180.0);
            if (currentVolumeOverlap > bestVolumeOverlap) {
                bestVolumeOverlap = currentVolumeOverlap;
                bestCoordsA = currentSDCoordA;
                bestI = i; bestJ = j;
            }
            if (verbose) cout << "finished round " << ++stepCount << " (A#" << i << " and B#" << j << ")" << endl;
        }
    }    
    cout << "\nThe best overlap is between conformer A#" << bestI << " and B#" << bestJ << ", which, after PCA followed by Steepest Descent, produces a volume overlap of " << bestVolumeOverlap << endl;
    
    //saveCoordsMatrixToMolecule(moleculeA, bestCoordsA);

    addConformerToMolecule(moleculeA, bestCoordsA); moleculeA.SetConformer(moleculeA.NumConformers() - 1);
    moleculeB.SetConformer(bestJ);
}


void runComparisons(int argc, char **argv) {
    if(argc < 4) { cout << "Usage: " << argv[0] << " ConformerSet1 ConformerSet2 OutputFileName\n"; exit(-1); }
    vector<OBMol> molecules;
    importMoleculeConformersFromFile(molecules, argv[1]);
    importMoleculeConformersFromFile(molecules, argv[2]);
    cout << "Finished importing molecules\n";

    runComformerComparisons(molecules[0], molecules[1]);

    cout << "\nSaving those best conformers to file..." << endl;
    writeMoleculeToFile(argv[3], molecules[0], true);
    writeMoleculeToFile(argv[3], molecules[1]);
}

void runRMSDTest(int argc, char **argv) {
    if(argc != 4) { cout << "Usage: " << argv[0] << " BeginningPosition Reference XRayDeterminedFinalPosition outputPDBFile\n"; exit(-1); }
    vector<OBMol> moleculesList;
    importMoleculesFromFile(moleculesList, argv[1]);
    importMoleculesFromFile(moleculesList, argv[2]);
    importMoleculesFromFile(moleculesList, argv[3]);
 
    PCAPlusSteepestDescent(moleculesList[0], moleculesList[1], 1.0, 10.0 * M_PI / 180.0);

    cout << "RMSD BETWEEN CRYSTAL STRUCtURE POSITION AND CALCULATED POSITION IS "<< calculateRMSD(moleculesList[0], moleculesList[2]) << endl;
    cout << "CALCULATED VOLUME OVERLAP IS " << volumeOverlap(moleculesList[0], moleculesList[1]) << " FOR THIS PROGRAM AND " << volumeOverlap(moleculesList[2], moleculesList[1]) << " FOR CHIMERA" << endl;
}

void runRMSDTest2(int argc, char **argv) {
    if(argc != 5) { cout << "Usage: " << argv[0] << " TargetBeginningPosition TargetBeginningPositionConformers Reference TargetXRayMatch\n"; exit(-1); }
    vector<OBMol> moleculesList;
    importMoleculesFromFile(moleculesList, argv[1]);
    importMoleculeConformersFromFile(moleculesList, argv[2]);
    importMoleculesFromFile(moleculesList, argv[3]);
    importMoleculesFromFile(moleculesList, argv[4]);

    PCAPlusSteepestDescent(moleculesList[0], moleculesList[2], 1.0, 10.0 * M_PI / 180.0);

    runComformerComparisons(moleculesList[1], moleculesList[2]);
    double bestConformerRMSD = calculateRMSD(moleculesList[1], moleculesList[3]);

    //cout << "conformer rmsd's:\n";
    cout << "conformer overlaps's:\n";
    for (unsigned int i=0; i < moleculesList[1].NumConformers() - 1; i++) {
        moleculesList[1].SetConformer(i);
        PCAPlusSteepestDescent(moleculesList[1], moleculesList[2], 1.0, 10.0 * M_PI / 180.0);
        //double conformerRMSD = calculateRMSD(moleculesList[1], moleculesList[3]);
        //cout << conformerRMSD << endl;
        cout << volumeOverlap(moleculesList[1], moleculesList[3]) << endl;

    }

    double RMSD = calculateRMSD(moleculesList[0], moleculesList[3]);

    cout << "RMSD BETWEEN CRYSTAL STRUCtURE FINAL POSITION AND CALCULATED POSITION OF TARGET IS "<< RMSD << endl 
         << "RMSD BETWEEN CRYSTAL STRUCtURE FINAL POSITION AND THE BEST CALCULATED CONFORMER POSITION OF TARGET IS "<< bestConformerRMSD << endl
         << "AN RMSD DIFFERENCE OF " << abs(bestConformerRMSD - RMSD) << " HAS BEEN OBSERVED.\n";
}


void runRMSDTest3(int argc, char **argv) {
    if(argc != 7) { cout << "Usage: " << argv[0] << " TargetBeginningPosition TargetBeginningPositionConformers Reference TargetXRayMatch originalTargetOuput conformerTargetOutput\n"; exit(-1); }
    vector<OBMol> moleculesList;
    importMoleculesFromFile(moleculesList, argv[1]);
    importMoleculeConformersFromFile(moleculesList, argv[2]);
    importMoleculesFromFile(moleculesList, argv[3]);
    importMoleculesFromFile(moleculesList, argv[4]);

    PCAPlusSteepestDescent(moleculesList[0], moleculesList[2], 1.0, 10.0 * M_PI / 180.0);
    runComformerComparisons(moleculesList[1], moleculesList[2]);
    double RMSD = calculateRMSD(moleculesList[0], moleculesList[3]);
    double bestConformerRMSD = calculateRMSD(moleculesList[1], moleculesList[3]);

    cout << "conformer RMSDs and overlaps:\n";
    for (unsigned int i=0; i < moleculesList[1].NumConformers() - 1; i++) {
        moleculesList[1].SetConformer(i); PCAPlusSteepestDescent(moleculesList[1], moleculesList[2], 1.0, 10.0 * M_PI / 180.0);
        cout << calculateRMSD(moleculesList[1], moleculesList[3]) << "\t" << volumeOverlap(moleculesList[1], moleculesList[2]) << endl;
    }

    cout << RMSD << "\t" << volumeOverlap(moleculesList[0], moleculesList[2]) << endl;
    cout << "0\t" << volumeOverlap(moleculesList[3], moleculesList[2]) << endl;

    cout << "RMSD BETWEEN CRYSTAL STRUCtURE FINAL POSITION AND THE BEST CALCULATED CONFORMER POSITION OF TARGET IS "<< bestConformerRMSD << endl;

/*    cout << "RMSD BETWEEN CRYSTAL STRUCtURE FINAL POSITION AND CALCULATED POSITION OF TARGET IS "<< RMSD << endl 
         << "RMSD BETWEEN CRYSTAL STRUCtURE FINAL POSITION AND THE BEST CALCULATED CONFORMER POSITION OF TARGET IS "<< bestConformerRMSD << endl
         << "VOLUME OVERLAP BETWEEN REFERENCE AND CALCULATED POSITION OF TARGET IS " << volumeOverlap(moleculesList[0], moleculesList[2]) << endl
         << "VOLUME OVERLAP BETWEEN THE REFERENCE AND THE BEST CALCULATED CONFORMER POSITION OF TARGET IS " << volumeOverlap(moleculesList[1], moleculesList[2]) << endl
         << "AN RMSD DIFFERENCE OF " << abs(bestConformerRMSD - RMSD) << " HAS BEEN OBSERVED.\n";
*/
    writeMoleculeToFile(argv[5], moleculesList[0], true);
    writeMoleculeToFile(argv[6], moleculesList[1], true);
}


void cleanUpAndGenerateConformers(int argc, char **argv) {
    if(argc != 4) { cout << "Usage: " << argv[0] << " MoleculeInputFile MoleculeOutputFile ConformerOutputFile\n"; exit(-1); }
    vector<OBMol> moleculesList;
    importMoleculesFromFile(moleculesList, argv[1]);

    removeNonBondedAtomsInMolecule(moleculesList[0]);
    moleculesList[0].AddHydrogens(false, true);
    writeMoleculeToFile(argv[2], moleculesList[0], true);

    //generateConformers(moleculesList[0]);
    //writeMoleculeConformersToFile(argv[3], moleculesList[0], true);
}




int main (int argc, char **argv) {
    //runComparisons(argc, argv);
    runRMSDTest3(argc, argv);


    //cleanUpAndGenerateConformers(argc, argv);



/*    vector<OBMol> moleculesList;
    importMoleculesFromFile(moleculesList, argv[1]);
    generateConformers(moleculesList[0]);
    writeMoleculeConformersToFile(argv[2], moleculesList[0], true);
*/    //cout << volumeOverlap (molecules[0], molecules[1]) << endl;


    //runComparisons(argc, argv);
    //if(argc < 3) { cout << "Usage: ProgrameName InputFileName OutputFileName\n"; return 1; }

    //vector<OBMol> molecules;
    //importMoleculesFromFile(molecules, argv[1]);
    //removeNonBondedAtomsInMolecule(molecules[0]);
    //generateConformers(molecules[0]);
    //writeMoleculeConformersToFile(argv[2], molecules[0], true);

/*    cout << "VEC SIZE: " << molecules.size() << endl;
    cout << "NUM CONFORMERS: " << molecules[0].NumConformers() << endl;
    cout << "NUM ATOMS: " << molecules[0].NumAtoms() << endl;
*/
    //molecules[1].SetConformer(10);


    //findBestPCAOrientation(molecules[0], molecules[1]);
    //runSteepestDescent(molecules[0], molecules[1], 0.5, 10.0 * M_PI / 180.0);
    //PCAPlusSteepestDescent(molecules[0], molecules[1], 0.5, 10.0 * M_PI / 180.0);

    return 0;
}
