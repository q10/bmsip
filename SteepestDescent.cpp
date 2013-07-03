inline void normalizeVector(vector<double> &vec) {
	double normalizationFactor = 1.0 / cblas_dnrm2(vec.size(), &vec[0], 1);
	for (unsigned int i=0; i < vec.size(); i++) vec[i] *= normalizationFactor;
}

inline void crossProduct(vector<double> &product, const vector<double> &U, const vector<double> &V) {
	if (U.size() != 3 or V.size() != 3) { cerr << "ERROR: size of U or V is not 3!"; abort(); }
	product.clear(); product.resize(3);
	product[0] = U[1] * V[2] - U[2] * V[1];
	product[1] = U[2] * V[0] - U[0] * V[2];
	product[2] = U[0] * V[1] - U[1] * V[0];
}

void generateSteepestDescentTranslationalVector(vector<double> &delT, const vector<double> &normalizedF, double h, double alphaTranslation) {
	double factor = h * alphaTranslation;
	delT.resize(normalizedF.size());
	for (unsigned int i=0; i < normalizedF.size(); i++) delT[i] = factor * normalizedF[i];
}

void generateMatrixWFromNormalizedVectorW(vector<double> &matW, const vector<double> &vecW) {
	if (vecW.size() != 3) { cerr << "ERROR: vecW size is not 3!"; abort(); }
	matW.resize(9); matW[0] = matW[4] = matW[8] = 0;
	matW[1] = vecW[2]; matW[2] = -vecW[1];
	matW[3] = -vecW[2]; matW[5] = vecW[0];
	matW[6] = vecW[1]; matW[7] = -vecW[0];
}

void generateSteepestDescentRotationMatrix(vector<double> &matR, const vector<double> &normalizedVectorT, double h, double beta) {
	if (normalizedVectorT.size() != 3) { cerr << "ERROR: vecW size is not 3!"; abort(); }

	double theta = h * beta;
	double tsin = sin(theta), tcos = 1 - cos(theta);
	matR.clear(); matR.resize(9, 0); matR[0] = matR[4] = matR[8] = 1;
	vector<double> matW;

	generateMatrixWFromNormalizedVectorW(matW, normalizedVectorT);

	for (unsigned int i=0; i < matR.size(); i++) matR[i] += tsin * matW[i];
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, tcos, &matW[0], 3, &matW[0], 3, 1, &matR[0], 3);

}

void calculateForceAndTorqueVectors(vector<double> &force, vector<double> &torque, vector<double> &coordsMoleculeA, vector<double> &coordsMoleculeB, const vector<double> &VDWsA, const vector<double> &VDWsB, const vector<double> &centerOfMassA, const vector< vector<double> > &atomMatchScoringTable) {
	if (coordsMoleculeA.size() != VDWsA.size() * 3 or coordsMoleculeB.size() != VDWsB.size() * 3) { cerr << "ERROR: INCORRECT MATCHING OF NUMBER OF COORDINATES AND ATOMIC NUMBERS; EXITING" << endl; abort(); }
	force.clear(); force.resize(3, 0); torque.clear(); torque.resize(3, 0);
	vector<double> differenceVector(3), cProduct(3), PmP(3);

    static const double constP = 2.0 * M_SQRT2;
    static const double A = 4.0 * M_PI * constP / 3.0;
    static const double B = -M_PI * pow(0.75 * constP * M_1_PI, 2.0/3.0);


	for (unsigned int i=0; i < coordsMoleculeA.size(); i+=3) {
		vector<double> delSdelPI(3, 0); // clear the partial derivative sum
        double vdwRA = VDWsA[i / 3]; // atomic number vectir is 3x shorter than 3d coordinates vector
        double sqvA = vdwRA * vdwRA;

		for (unsigned int j=0; j < coordsMoleculeB.size(); j+=3) {
            double vdwRB = VDWsB[j / 3];
            double sqvB = vdwRB * vdwRB;
            double C = sqvA + sqvB;
            double D = B / C;

			for (unsigned int k=0; k<3; k++) differenceVector[k] = coordsMoleculeA[i+k] - coordsMoleculeB[j+k];
            double distanceSquared = pow( cblas_dnrm2(differenceVector.size(), &differenceVector[0], 1), 2 );
			double expPmQ2 = 2 * D * A * pow(sqvA * sqvB  / C, 1.5) * exp( D * distanceSquared ); // the partial derivative evaluated
			for (unsigned int k=0; k<3; k++) delSdelPI[k] +=  atomMatchScoringTable[i/3][j/3] * expPmQ2 * differenceVector[k];
		}
		
		for (unsigned int k=0; k<3; k++) force[k] += delSdelPI[k]; // add to force
		
		for (unsigned int k=0; k<3; k++) PmP[k] = coordsMoleculeA[i+k] - centerOfMassA[k];
		crossProduct(cProduct, PmP, delSdelPI);
		for (unsigned int k=0; k<3; k++) torque[k] += cProduct[k]; // add to torque
	}
}

double steepestDescentEngine(vector<double> &finalCoordsA, vector<double> &coordsMoleculeA, vector<double> &coordsMoleculeB, vector<double> &VDWsA, vector<double> &VDWsB, vector<double> &centerOfMassA, const vector< vector<double> > &atomMatchScoringTable, double alphaTranslation, double betaRadians, bool verbose=false) {
	finalCoordsA = coordsMoleculeA;
	vector<double> bestAInThisStep, force, torque, delT, delR, tempA;
	double currentStepH = 1, bestVolumeOverlapSoFar = volumeOverlap(finalCoordsA, coordsMoleculeB, VDWsA, VDWsB, atomMatchScoringTable);

    while (currentStepH > 0) {
    	if (verbose) cout << "stepping... ";
		double bestHInThisStep = 0;
		bestAInThisStep = finalCoordsA;

		// calculate force and torque, and normalize them
		calculateForceAndTorqueVectors(force, torque, finalCoordsA, coordsMoleculeB, VDWsA, VDWsB, centerOfMassA, atomMatchScoringTable);
		normalizeVector(force); normalizeVector(torque);

		for (double h=0; h < 1.0; h += 0.01) {
			generateSteepestDescentTranslationalVector(delT, force, h, alphaTranslation); // generate deviation translational vector
			generateSteepestDescentRotationMatrix(delR, torque, h, betaRadians); // generate deviation rotational matrix

			// rotate and translate
			tempA = finalCoordsA;
			translate3DMatrixCoordinates(tempA, -centerOfMassA[0], -centerOfMassA[1], -centerOfMassA[2]);
			rotate3DMatrixCoordinates(tempA, delR);
			translate3DMatrixCoordinates(tempA, centerOfMassA[0] + delT[0], centerOfMassA[1] + delT[1], centerOfMassA[2] + delT[2]); // move coordinates back FROM center of rotation and add the translational displacement

		    double currentVolOverlap = volumeOverlap(tempA, coordsMoleculeB, VDWsA, VDWsB, atomMatchScoringTable);
		    if (currentVolOverlap > bestVolumeOverlapSoFar) { // check overlap goodness
		    	bestHInThisStep = h;
		    	bestVolumeOverlapSoFar = currentVolOverlap;
		    	bestAInThisStep = tempA;
		    }
		}
		finalCoordsA = bestAInThisStep; currentStepH = bestHInThisStep;
		if (verbose) cout << "Best h in this step is " << bestHInThisStep << " and volume overlap is at " << bestVolumeOverlapSoFar << endl;
	}
	return bestVolumeOverlapSoFar;
}

void runSteepestDescent(OBMol &moleculeA, OBMol &moleculeB, double alphaTranslation, double betaRadians, bool verbose) {
    if (verbose) cout << endl << "BEGIN STEEPEST DESCENT SEARCH" << endl << "alphaTranslation = " << alphaTranslation << endl << "beta = " << betaRadians << endl << "Initializing data..." << endl;

	vector<double> coordsMoleculeA, coordsMoleculeB, centerOfMassA, finalCoordsA, VDWsA, VDWsB;
    vector< vector<double> > atomMatchScoringTable;

    generateCoordsMatrixFromMolecule(coordsMoleculeA, moleculeA);
    generateCoordsMatrixFromMolecule(coordsMoleculeB, moleculeB);
    generateVDWRadiusListFromMolecule(VDWsA, moleculeA);
    generateVDWRadiusListFromMolecule(VDWsB, moleculeB);
    getMoleculeCenterCoords(centerOfMassA, moleculeA);
    generateAtomMatchScoringTableFromTwoMolecules(atomMatchScoringTable, moleculeA, moleculeB);


    double bestVolumeOverlap = steepestDescentEngine(finalCoordsA, coordsMoleculeA, coordsMoleculeB, VDWsA, VDWsB, centerOfMassA, atomMatchScoringTable, alphaTranslation, betaRadians, verbose);

    if (verbose) {
	    cout << "\nThe convergent solution coordinates of A produces a volume overlap of " << bestVolumeOverlap << endl;
	    cout << "\nRESULTING A:" << endl; printMatrix(finalCoordsA, finalCoordsA.size() / 3, 3, false);
	    cout << "\nEND STEEOEST DESCENT SEARCH.  SAVING COORDINATES TO MOLECULE A...\n\n";
	}
	saveCoordsMatrixToMolecule(moleculeA, finalCoordsA);
}
