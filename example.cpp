#include <iostream>
#include <cmath>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/data.h>

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
