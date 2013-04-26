#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

using namespace std;
using namespace OpenBabel;

/*
double volumeOverlap (OBMol &moleculeA, OBMol &moleculeB) {
    for (OBAtomIterator iterA = moleculeA.BeginAtoms(); iterA != moleculeA.EndAtoms(); iterA++) {
        double *coordsOfAtomI = (*iterA)->GetCoordinate();

        for (OBAtomIterator iterB = moleculeB.BeginAtoms(); iterB != moleculeB.EndAtoms(); iterB++) {
            double *coordsOfAtomJ = (*iter)->GetCoordinate();

            4/3 pi P ( R1^2 R2^2  / (R1^2+R2^2) )^(3/2)    exp(    -pi (3P/4pi)^(2/3)   (ri - rj)^2 /(R1^2+R2^2)    )


        }
    }
}
*/



int main (int argc, char **argv) {
    if(argc<2) {
        cout << "Usage: ProgrameName InputFileName\n";
        return 1;
    }



    OBConversion obconversion;
    obconversion.SetInFormat("sdf");

    OBMol mol;

    bool notatend = obconversion.ReadFile(&mol, argv[1]);
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
