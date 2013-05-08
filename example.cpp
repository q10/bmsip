#include <iostream>
#include <limits>
#include <cmath>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/data.h>

#ifdef __APPLE__
    #include <Accelerate/Accelerate.h>
#else
    #include <cblas.h>
    #include <clapack.h>
    #include <f2c.h>
    #include <blaswrap.h>
#endif

using namespace std;
using namespace OpenBabel;

void printMatrix(vector<double> &matrix, unsigned int rows, unsigned int columns, bool columnMajorOrder = true) {
    if (matrix.size() != rows * columns) { cerr << "ERROR: INCORRECT MATCHING OF ROWS AND COLUMNS WITH ACTUAL VECTOR SIZE; EXITING" << endl; abort(); }
    if (columnMajorOrder) {
        for (unsigned int i=0; i < rows; i++) {
            for (unsigned int j=i; j < matrix.size(); j+=rows)
                cout << matrix[j] << " ";                
            cout << endl;
        }
    } else {
        for (unsigned int i=0; i < matrix.size(); i+=columns) {
            for (unsigned int j=0; j < columns; j++)
                cout << matrix[i + j] << " ";
            cout << endl;
        }
    }
}

void generateCoordsMatrixFromMolecule(vector<double> &matrix, OBMol &molecule) {     // generates column-order matrix of coordinates
    matrix.clear(); matrix.insert(matrix.end(), molecule.GetCoordinates(), &molecule.GetCoordinates()[3*molecule.NumAtoms()]);
}

void saveCoordsMatrixToMolecule(OBMol &molecule, vector<double> &matrix) {
    if (matrix.size() != molecule.NumAtoms() * 3) { cerr << "ERROR: INCORRECT MATCHING OF NUMBER OF COORDINATES; EXITING" << endl; abort(); }
    molecule.SetCoordinates(&matrix[0]);
}

void generateAtomicNumbersListFromMolecule(vector<int> &numList, OBMol &molecule) {
    numList.clear();
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++)
        numList.push_back((*iter)->GetAtomicNum());
}


double volumeOverlap(vector<double> &coordsMoleculeA, vector<double> &coordsMoleculeB, vector<int> &atomicNumbersA, vector<int> &atomicNumbersB) {
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
            totalVolumeOverlap += A * pow(sqvA * sqvB  / C, 1.5) * exp(B * distanceSquared / C );
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
    int matrixOrder = 3;
    eigenvalues.clear(); eigenvalues.resize(3, 0);
    eigenvectors.clear(); eigenvectors = matrix;
    int lwork = 102, info; vector<double> work(lwork); // LAPACK parameters; work[0] will show the optimum value for lwork when call completes
    if (dsyev_("V", "L", &matrixOrder, &eigenvectors[0], &matrixOrder, &eigenvalues[0], &work[0], &lwork, &info)) { cerr << "ERROR: NO MATRIX INITIALIZED IN POINTER; EXITING" << endl; abort(); }
}

void getMoleculeCenterCoords(vector<double> &centerCoords, OBMol &molecule) {
    centerCoords.clear(); centerCoords.resize(3, 0);
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++) {
        double *tmpCoords = (*iter)->GetCoordinate(); double tmpMass = (*iter)->GetAtomicMass();  
        centerCoords[0] += tmpMass*tmpCoords[0]; centerCoords[1] += tmpMass*tmpCoords[1]; centerCoords[2] += tmpMass*tmpCoords[2];
    }
    for (unsigned int i=0; i < centerCoords.size(); i++) centerCoords[i] /= molecule.GetMolWt();
}

inline void translate3DMatrixCoordinates(vector<double> &matrix, double x, double y, double z) {
    for (unsigned int i=0; i < matrix.size(); i+=3) { matrix[i] += x; matrix[i+1] += y; matrix[i+2] += z; }
}

void rotate3DMatrixCoordinates(vector<double> &matrix, vector<double> &rotationMatrix) {
    // both matrices must be of column-order
    vector<double> resultMatrix(matrix.size(), 0);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, matrix.size()/3, 3, 1, &rotationMatrix[0], 3, &matrix[0], 3, 0, &resultMatrix[0], 3);
    matrix = resultMatrix;
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

void generateOptimalRotationMatrix(vector<double> &rotMatrix, unsigned int optCode, const vector<double> &U, const vector<double> &V) {
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



void findBestInitialOrientation(OBMol &moleculeA, OBMol &moleculeB) {
    cout << endl << "BEGIN INITIAL ORIENTATION SEARCH" << endl
        << "Searching for the best initial orientation matrix..." << endl;
    vector<double> coordA, coordB, comA, comB, covA, covB, eVectA, eVectB, eValA, eValB, tempR, bestR;
    vector<int> atomicNumsA, atomicNumsB;
    map<int, string> RTable; RTable[0] = "R0"; RTable[1] = "Rx"; RTable[2] = "Ry"; RTable[3] = "Rz";

    generateCoordsMatrixFromMolecule(coordA, moleculeA);
    generateCoordsMatrixFromMolecule(coordB, moleculeB);

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
        generateOptimalRotationMatrix(tempR, i, eVectA, eVectB);
        vector<double> tempA = coordA;

        translate3DMatrixCoordinates(tempA, -comA[0], -comA[1], -comA[2]);
        rotate3DMatrixCoordinates(tempA, tempR);
        translate3DMatrixCoordinates(tempA, comB[0], comB[1], comB[2]);

        double curVolOverlap = volumeOverlap(tempA, coordB, atomicNumsA, atomicNumsB);
        if (curVolOverlap > bestVolumeOverlap) {
            bestRcode = i;
            bestVolumeOverlap = curVolOverlap;
            bestR = tempR;
        }
    }

    cout << "The best initial orientation matrix is: " << RTable[bestRcode] << ", which produced a volume overlap of " << bestVolumeOverlap << endl;
    printMatrix(bestR, 3, 3);
    cout << "END INITIAL ORIENTATION SEARCH" << endl << endl;
/*
    generateOptimalRotationMatrix(Rx, 1, eVectA, eVectB);
    generateOptimalRotationMatrix(Ry, 2, eVectA, eVectB);
    generateOptimalRotationMatrix(Rz, 3, eVectA, eVectB);


    // the transformation itself, R(p - p0) + q0
    translate3DMatrixCoordinates(coordA, -comA[0], -comA[1], -comA[2]);
    rotate3DMatrixCoordinates(coordA, Rx);
    translate3DMatrixCoordinates(coordA, comB[0], comB[1], comB[2]);

    saveCoordsMatrixToMolecule(moleculeA, coordA);
    OBConversion obconversion;
    obconversion.SetOutFormat("sdf");


    obconversion.WriteFile (&moleculeA, "newA.sdf");
    cout << "OOOOKKKKKKK" << endl;
*/
}







void sampleTest(OBMol &molecule) {
    cout << "BESGIN TEST" << endl;
    vector<double> abc;
    generateCovarMatrixFromMolecule(abc, molecule);
    for (int i=0; i<abc.size(); i+=3) { cout << "{" << abc[i] << "," << abc[i+1] << "," << abc[i+2] << "},"; } cout << endl;
    vector<double> eigenvectors, eigenvalues;
    generateEigenMatrix(eigenvectors, eigenvalues, abc);

    for (int i=0; i<eigenvalues.size(); i++) cout << eigenvalues[i] << " ";
        cout << endl;
    for (int i=0; i<eigenvectors.size(); i++) cout << eigenvectors[i] << " ";
        cout << endl;
    cout << "END TEST" << endl << endl;
}


void testEigen() {
    double xyz[] = {7, 9, 2, 9, 1, 6, 2, 6, 10}; 
    vector<double> sample(xyz, &xyz[9]);
    for (int i=0; i<sample.size(); i++) cout << sample[i] << " ";
        cout << endl;

    vector<double> eigenvectors, eigenvalues;
    generateEigenMatrix(eigenvectors, eigenvalues, sample);
    for (int i=0; i<eigenvalues.size(); i++) cout << eigenvalues[i] << " ";
        cout << endl;
    for (int i=0; i<eigenvectors.size(); i++) cout << eigenvectors[i] << " ";
        cout << endl << endl << endl << endl;
}

void testGenRot() {
    double a[] = {1, 2, 3, 2, 4, 5, 3, 5, 6}; 
    vector<double> vA(a, &a[9]);
    double b[] = {7, 9, 2, 2, 1, 0, 2, 6, 10}; 
    vector<double> vB(b, &b[9]);
    vector<double> c;
    generateOptimalRotationMatrix(c, 3, vA, vB);
    cout << endl << endl;
    for (int i=0; i<vA.size(); i++) cout << vA[i] << " ";
        cout << endl;
    for (int i=0; i<vB.size(); i++) cout << vB[i] << " ";
        cout << endl;
    for (int i=0; i<c.size(); i++) cout << c[i] << " ";
        cout << endl;
}

void testRot() {
    cout << "BEGIN ROTATE TEST" << endl;
    double a[] = {1, 2, 3, 2, 4, 5, 3, 5, 6}; 
    vector<double> vA(a, &a[9]);
    double b[] = {7, 9, 2, 2, 1, 0, 2, 6, 10}; 
    vector<double> vB(b, &b[9]);
    vector<double> c;
    for (int i=0; i<vA.size(); i++) cout << vA[i] << " ";
        cout << endl;
    for (int i=0; i<vB.size(); i++) cout << vB[i] << " ";
        cout << endl;
    rotate3DMatrixCoordinates(vA, vB);
    for (int i=0; i<vA.size(); i++) cout << vA[i] << " ";
        cout << endl;
    for (int i=0; i<vB.size(); i++) cout << vB[i] << " ";
        cout << endl;
    vB.push_back(199);
    cout << endl;
    printMatrix(vB, 5, 2);
    cout << "END ROTATE TEST" << endl;
}




int main (int argc, char **argv) {
    if(argc < 3) {
        cout << "Usage: ProgrameName InputFileName InputFileName2\n";
        return 1;
    }



    OBConversion obconversion;
    obconversion.SetInFormat("sdf");

    OBMol mol2;
    obconversion.ReadFile(&mol2, argv[2]);

    OBMol mol;
    bool notatend = obconversion.ReadFile(&mol, argv[1]);

    cout << "VOL OVERLAP = " << volumeOverlap (mol, mol2) << endl;


    sampleTest(mol);
    findBestInitialOrientation(mol, mol2);


    while (notatend) {
        std::cout << "Molecular Weight: " << mol.GetMolWt() << std::endl;
        for (OBAtomIterator iter = mol.BeginAtoms(); iter != mol.EndAtoms(); iter++) {
            cout << (*iter)->GetVector() << ", " << (*iter)->GetCoordinate()[0] << endl;
        }
        mol.Clear();
        notatend = obconversion.Read(&mol);
    }
    
    testRot();
    //testEigen();

    return 0;
}

/*

// older implementation of volumeOverlap

double volumeOverlap (OBMol &moleculeA, OBMol &moleculeB) {
    double totalVolumeOverlap = 0;

    const double constP = 2.0 * M_SQRT2;
    const double A = 4.0 * M_PI * constP / 3.0;
    const double B = -M_PI * pow(0.75 * constP * M_1_PI, 2.0/3.0);


    for (OBAtomIterator iterA = moleculeA.BeginAtoms(); iterA != moleculeA.EndAtoms(); iterA++) {
        double *coordsOfAtomI = (*iterA)->GetCoordinate();
        double vdwRA = etab.GetVdwRad((*iterA)->GetAtomicNum());

        for (OBAtomIterator iterB = moleculeB.BeginAtoms(); iterB != moleculeB.EndAtoms(); iterB++) {
            double *coordsOfAtomJ = (*iterB)->GetCoordinate();
            double vdwRB = etab.GetVdwRad((*iterB)->GetAtomicNum());
            
            double sqvA = vdwRA * vdwRA;
            double sqvB = vdwRB * vdwRB;
            double C = sqvA + sqvB;

            double distanceSquared = pow(coordsOfAtomJ[0]-coordsOfAtomI[0], 2) + pow(coordsOfAtomJ[1]-coordsOfAtomI[1], 2) + pow(coordsOfAtomJ[2]-coordsOfAtomI[2], 2);


            totalVolumeOverlap += A * pow(sqvA * sqvB  / C, 1.5) * exp(B * distanceSquared / C );
        }
    }
    return totalVolumeOverlap;
}
*/



/*
   ifstream ifs(argv[1]);
    if(!ifs) {
        cout << "Cannot open input file\n";
        return 1;
    }

    ofstream ofs(argv[2]);
    if(!ofs) {
        cout << "Cannot open output file\n";
        return 1;
    }
    OpenBabel::OBConversion conv(&ifs, &ofs);
    if(!conv.SetInAndOutFormats("CML","MOL")) {
        cout << "Formats not available\n";
        return 1;
    }
    int n = conv.Convert();
    cout << n << " molecules converted\n";

*/
