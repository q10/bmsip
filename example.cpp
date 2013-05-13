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

void printMoleculeCoords(OBMol &molecule) {
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++) cout << (*iter)->GetVector() << endl;
}

unsigned int importMoleculesFromFile(vector<OBMol> &moleculesList, const string &fileName, const string &format) { // does not clear list, only appends new molecules to it
    OBConversion obconversion;
    obconversion.SetInFormat(format.c_str());

    unsigned int numMoleculesInFile = 0;
    OBMol mol; bool notAtEnd = obconversion.ReadFile(&mol, fileName.c_str());
    while (notAtEnd) {
        numMoleculesInFile++;
        moleculesList.push_back(mol);
        mol.Clear(); notAtEnd = obconversion.Read(&mol);
    }
}

void writeMoleculeToFile(const string &fileName, const string &format, OBMol &molecule, bool rewriteFile=false) {
    ios_base::openmode fileMode = rewriteFile ? ios::out : (ios::out|ios::app);
    std::ofstream ofs(fileName.c_str(), fileMode);
    OBConversion obconversion; obconversion.SetOutFormat(format.c_str());
    obconversion.Write(&molecule, &ofs);  // obconversion.WriteFile(&molecule, fileName.c_str());
}

void writeAllMoleculeConformersToFile(const string &fileName, const string &format, OBMol &molecule) {
    for (int i = molecule.NumConformers()-1; i >= 0; i--) { // DO NOT USE UNSIGNED INT i!!!
        molecule.SetConformer(i);
        writeMoleculeToFile(fileName, format, molecule);
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
            bestA = tempA;
        }
    }

    cout << "\nThe best initial orientation matrix is: " << RTable[bestRcode] << ", which produces a volume overlap of " << bestVolumeOverlap << endl;
    printMatrix(bestR, 3, 3);
    cout << "\nRESULTING A:" << endl; printMatrix(bestA, bestA.size() / 3, 3, false); 
    cout << "\nEND INITIAL ORIENTATION SEARCH.  SAVING COORDINATES TO MOLECULE A...\n\n";
    saveCoordsMatrixToMolecule(moleculeA, bestA);


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


    obconversion.WriteFile(&moleculeA, "newA.sdf");
    cout << "OOOOKKKKKKK" << endl;
*/
}


void generateConformers(OBMol &molecule) {
    OBConformerSearch cs;
    cs.Setup(molecule); // numConformers 30 // numChildren 5 // mutability 5 // convergence 25
    //cs.SetScore(new OBEnergyConformerScore);
    cs.Search();

/*
    cout << "Setting up conformer searching..." << endl
            << "   conformers:  30" << endl
            << "   children:    5" << endl
            << "   mutability:  5" << endl
            << "   convergence: 25" << endl;

    // Print the rotor keys
    RotorKeys keys = cs.GetRotorKeys();
    for (RotorKeys::iterator key = keys.begin(); key != keys.end(); ++key) {
        for (unsigned int i = 1; i < key->size(); ++i) cout << key->at(i) << " ";
        cout << endl;
    }
*/
    // Get the conformers
    cout << "NUM CONFORMERS: " << molecule.NumConformers() << endl;
    cs.GetConformers(molecule);
    cout << "NUM CONFORMERS: " << molecule.NumConformers() << endl << endl;
    printMoleculeCoords(molecule); cout << endl << endl;
    molecule.SetConformer(0);
    printMoleculeCoords(molecule); cout << endl << endl;
    molecule.SetConformer(11);
    printMoleculeCoords(molecule);
       
}



#include "SteepestDescent.cpp"




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
    if(argc < 4) {
        cout << "Usage: ProgrameName InputFileName InputFileName2 OutputFileName\n";
        return 1;
    }

    vector<OBMol> molecules;
    importMoleculesFromFile(molecules, argv[1], "sdf");
    importMoleculesFromFile(molecules, argv[2], "sdf");

    cout << "VOL OVERLAP = " << volumeOverlap(molecules[0], molecules[1]) << endl;

    //findBestInitialOrientation(molecules[0], molecules[1]);
    //generateConformers(molecules[1]);
    //molecules[1].SetConformer(10);


    runSteepestDescent(molecules[0], molecules[1], 0.5, 10.0 * M_PI / 180.0);
    //cout << "WRITING MOLECULE A TO FILE '" << argv[3] << "' (WILL OVERWRITE EXISTING FILE IF ANY)...\n";
    //writeMoleculeToFile(argv[3], "sdf", molecules[0]);

    //generateConformers(molecules[1]);
    //writeAllMoleculeConformersToFile(argv[3], "sdf", molecules[1]);




    //sampleTest(mol);
    //testRot();
    //testEigen();

    return 0;
}
