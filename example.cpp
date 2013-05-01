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
#endif

using namespace std;
using namespace OpenBabel;


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



void generateCoordsMatrixFromMolecule(vector<double> &matrix, OBMol &molecule) {     // generates column-order matrix of coordinates
    matrix.clear(); matrix.insert(matrix.end(), molecule.GetCoordinates(), &molecule.GetCoordinates()[3*molecule.NumAtoms()]);
}

void saveCoordsMatrixToMolecule(OBMol &molecule, vector<double> &matrix) {
    //if (matrix.size()/3 != )
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

void translate3DMatrixCoordinates(vector<double> &matrix, double x, double y, double z) {
    for (unsigned int i=0; i < matrix.size(); i+=3) { matrix[i] += x; matrix[i+1] += y; matrix[i+2] += z; }
}

void rotate3DMatrixCoordinates(vector<double> *matrix, vector<double> &rotationMatrix) {
    // both matrices must be of column-order
    if (not matrix) { cerr << "ERROR: NO MATRIX INITIALIZED IN POINTER; EXITING" << endl; abort(); }
    vector<double> *resultMatrix = new vector<double>(matrix->size());
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, matrix->size()/3, 3, 1, &rotationMatrix[0], 3, &(*matrix)[0], 3, 0, &(*resultMatrix)[0], 3);
    delete matrix; matrix = resultMatrix;
}

void generateCovarMatrixFromMolecule(vector<double> &matrix, OBMol &molecule) {
    double uX = 0, uY = 0, uZ = 0, *moleculeCoords = molecule.GetCoordinates();
    for (unsigned int i=0; i < 3 * molecule.NumAtoms(); i+=3) { uX += moleculeCoords[i]; uY += moleculeCoords[i+1]; uZ += moleculeCoords[i+2]; }
    uX /= molecule.NumAtoms(); uY /= molecule.NumAtoms(); uZ /= molecule.NumAtoms();

    double cXX = 0, cYY = 0, cZZ = 0, cXY = 0, cXZ = 0, cYZ = 0;
    for (unsigned int i=0; i < 3 * molecule.NumAtoms(); i+=3) {
        double atomWeight = molecule.GetAtom(i/3)->GetAtomicMass();
        cXX += pow(moleculeCoords[i] - uX, 2) * atomWeight;
        cYY += pow(moleculeCoords[i+1] - uY, 2) * atomWeight;
        cZZ += pow(moleculeCoords[i+2] - uZ, 2) * atomWeight;
        cXY += (moleculeCoords[i] - uX) * (moleculeCoords[i+1] - uY) * atomWeight;
        cXZ += (moleculeCoords[i] - uX) * (moleculeCoords[i+2] - uZ) * atomWeight;
        cYZ += (moleculeCoords[i+1] - uY) * (moleculeCoords[i+2] - uZ) * atomWeight;
    }
    
    matrix.resize(9);
    matrix[0] = cXX; 
    matrix[1] = matrix[3] = cXY; 
    matrix[2] = matrix[6] = cXZ; 
    matrix[4] = cYY;
    matrix[5] = matrix[7] = cYZ;  
    matrix[8] = cZZ;

    for (unsigned int i=0; i < matrix.size(); i++) matrix[i] /= molecule.GetMolWt();
}




int main (int argc, char **argv) {
    if(argc < 3) {
        cout << "Usage: ProgrameName InputFileName InputFileName2\n";
        return 1;
    }


    double xyz[] = {1, 2, 3, 2, 4, 5, 3, 5, 6}; 
    vector<double> sample(xyz, &xyz[9]);
    for (int i=0; i<sample.size(); i++) cout << sample[i] << " ";
        cout << endl << endl << endl << endl;

    vector<double> eigenvectors, eigenvalues;
    generateEigenMatrix(eigenvectors, eigenvalues, sample);
    for (int i=0; i<eigenvalues.size(); i++) cout << eigenvalues[i] << " ";
        cout << endl << endl << endl << endl;
    for (int i=0; i<eigenvectors.size(); i++) cout << eigenvectors[i] << " ";
        cout << endl << endl << endl << endl;


    OBConversion obconversion;
    obconversion.SetInFormat("sdf");

    OBMol mol2;
    obconversion.ReadFile(&mol2, argv[2]);

    OBMol mol;
    bool notatend = obconversion.ReadFile(&mol, argv[1]);

    cout << "VOL OVERLAP = " << volumeOverlap (mol, mol2) << endl;

    while (notatend) {
        std::cout << "Molecular Weight: " << mol.GetMolWt() << std::endl;
        for (OBAtomIterator iter = mol.BeginAtoms(); iter != mol.EndAtoms(); iter++) {
            cout << (*iter)->GetVector() << ", " << (*iter)->GetCoordinate()[0] << endl;
        }
        mol.Clear();
        notatend = obconversion.Read(&mol);
    }
    
    return 0;
}

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
