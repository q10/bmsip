void steepestDescentEngine(OBMol &moleculeA, OBMol &moleculeB, double alpha, double betaRadians) {
    cout << endl << "BEGIN STEEPEST DESCENT SEARCH" << endl
        << "alpha = " << alpha << endl << "beta = " << betaRadians << endl << "Initializing data..." << endl;

    vector<double> coordsMoleculeA, coordsMoleculeB, centerOfMassA, force, torque, delT, delR, tempA, bestAInThisStep;
    vector<int> atomicNumsA, atomicNumsB;

    generateCoordsMatrixFromMolecule(coordsMoleculeA, moleculeA);
    generateCoordsMatrixFromMolecule(coordsMoleculeB, moleculeB);
    generateAtomicNumbersListFromMolecule(atomicNumsA, moleculeA);
    generateAtomicNumbersListFromMolecule(atomicNumsB, moleculeB);
    getMoleculeCenterCoords(centerOfMassA, moleculeA);

    double currentStepH = 1, bestVolumeOverlapSoFar = volumeOverlap(coordsMoleculeA, coordsMoleculeB, atomicNumsA, atomicNumsB);

    while (currentStepH > 0) {
        cout << "stepping... ";
        double bestHInThisStep = 0;
        bestAInThisStep = coordsMoleculeA;

        // calculate force and torque, and normalize them
        calculateForceAndTorqueVectors(force, torque, coordsMoleculeA, coordsMoleculeB, atomicNumsA, atomicNumsB, centerOfMassA);
        normalizeVector(force); normalizeVector(torque);

        for (double h=0; h < 1.0; h += 0.01) {
            generateSteepestDescentTranslationalVector(delT, force, h, alpha); // generate deviation translational vector
            generateSteepestDescentRotationMatrix(delR, torque, h, betaRadians); // generate deviation rotational matrix

            // rotate and translate
            tempA = coordsMoleculeA;
            translate3DMatrixCoordinates(tempA, -centerOfMassA[0], -centerOfMassA[1], -centerOfMassA[2]);
            rotate3DMatrixCoordinates(tempA, delR);
            translate3DMatrixCoordinates(tempA, centerOfMassA[0] + delT[0], centerOfMassA[1] + delT[1], centerOfMassA[2] + delT[2]); // move coordinates back FROM center of rotation and add the translational displacement

            double currentVolOverlap = volumeOverlap(tempA, coordsMoleculeB, atomicNumsA, atomicNumsB);
            if (currentVolOverlap > bestVolumeOverlapSoFar) { // check overlap goodness
                bestHInThisStep = h;
                bestVolumeOverlapSoFar = currentVolOverlap;
                bestAInThisStep = tempA;
            }
        }
        coordsMoleculeA = bestAInThisStep; currentStepH = bestHInThisStep;
        cout << "Best h in this step is " << bestHInThisStep << " and volume overlap is at " << bestVolumeOverlapSoFar << endl;
    }

    cout << "\nThe convergent solution coordinates of A produces a volume overlap of " << bestVolumeOverlapSoFar << endl;
    cout << "\nRESULTING A:" << endl; printMatrix(coordsMoleculeA, coordsMoleculeA.size() / 3, 3, false);

    cout << "\nEND STEEOEST DESCENT SEARCH.  SAVING COORDINATES TO MOLECULE A...\n\n";
    saveCoordsMatrixToMolecule(moleculeA, coordsMoleculeA);
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
    generatePCARotationMatrix(c, 3, vA, vB);
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
    /*rotate3DMatrixCoordinates(vA, vB);
    for (int i=0; i<vA.size(); i++) cout << vA[i] << " ";
        cout << endl;
    for (int i=0; i<vB.size(); i++) cout << vB[i] << " ";
        cout << endl;
    vB.push_back(199);
    cout << endl;
    printMatrix(vB, 5, 2);
    */
    vector<double> C; vA.resize(3); vB.resize(3);
    generateMatrixWFromNormalizedVectorW(C, vB);
    for (int i=0; i<C.size(); i++) cout << C[i] << " ";
        cout << endl;
    cout << "END ROTATE TEST" << endl;
}


/*
    generatePCARotationMatrix(Rx, 1, eVectA, eVectB);
    generatePCARotationMatrix(Ry, 2, eVectA, eVectB);
    generatePCARotationMatrix(Rz, 3, eVectA, eVectB);


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
    OBConversion obconversion;
    obconversion.SetInFormat("sdf");

    OBMol mol2;
    obconversion.ReadFile(&mol2, argv[2]);

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


/*
    // Print the rotor keys
    RotorKeys keys = cs.GetRotorKeys();
    for (RotorKeys::iterator key = keys.begin(); key != keys.end(); ++key) {
        for (unsigned int i = 1; i < key->size(); ++i) cout << key->at(i) << " ";
        cout << endl;
    }
*/

    //sampleTest(mol);
    //testRot();
    //testEigen();

/*
    vector< vector<double> > coordAs, coordBs, comAs, comBs, eVectAs, eVectBs;
    vector<double> VDWsA, VDWsB, massesA, massesB, currentPCACoordA, currentSDCoordA, bestCoordsA;

    double molecularWeightA = molecules[0].GetMolWt(), molecularWeightB = molecules[1].GetMolWt();
    double bestVolumeOverlap = -1;
    int bestJ = -1, bestI = -1, stepCount = 0;
    
    generateVDWRadiusListFromMolecule(VDWsA, molecules[0]);
    generateVDWRadiusListFromMolecule(VDWsB, molecules[1]);
    generateAtomicMassesListFromMolecule(massesA, molecules[0]);
    generateAtomicMassesListFromMolecule(massesB, molecules[1]);
    generateCoordsMatrixFromMoleculeConformers(coordAs, molecules[0]);
    generateCoordsMatrixFromMoleculeConformers(coordBs, molecules[1]);
    getMoleculeConformerCenterCoords(comAs, molecules[0]);
    getMoleculeConformerCenterCoords(comBs, molecules[1]);    
    eVectAs.resize( molecules[0].NumConformers() ), eVectBs.resize( molecules[1].NumConformers() );
    for (unsigned int k=0; k < molecules[0].NumConformers(); k++) getPCAEigenMatrix(eVectAs[k], coordAs[k], massesA, molecularWeightA);
    for (unsigned int k=0; k < molecules[1].NumConformers(); k++) getPCAEigenMatrix(eVectBs[k], coordBs[k], massesB, molecularWeightB);
    cout << "Finished setting up data; running search...\n";

    for (unsigned int j=0; j < molecules[1].NumConformers(); j++) {
        for (unsigned int i=0; i < molecules[0].NumConformers(); i++) {
            PCAEngine(currentPCACoordA, coordAs[i], coordBs[j], eVectAs[i], eVectBs[j], comAs[i], comBs[j], VDWsA, VDWsB);
            double currentVolumeOverlap = steepestDescentEngine(currentSDCoordA, currentPCACoordA, coordBs[j], VDWsA, VDWsB, comAs[i], 1.0, 10.0 * M_PI / 180.0);
            if (currentVolumeOverlap > bestVolumeOverlap) {
                bestVolumeOverlap = currentVolumeOverlap;
                bestCoordsA = currentSDCoordA;
                bestI = i; bestJ = j;
            }
            cout << "finished round " << ++stepCount << endl;
        }
    }

    cout << "\nThe best overlap is between conformer A#" << bestI << " and B#" << bestJ << ", which, after PCA followed by Steepest Descent, produces a volume overlap of " << bestVolumeOverlap << endl;
    cout << "\nSaving those best conformers to file..." << endl;

    saveCoordsMatrixToMolecule(molecules[0], bestCoordsA);
    molecules[1].SetConformer(bestJ);
*/
