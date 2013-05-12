inline void normalizeVector(vector<double> &vec) {
	double normalizationFactor = 1.0 / cblas_dnrm2(vec.size(), &vec[0], 1);
	for (unsigned int i=0; i < vec.size(); i++) vec[i] *= normalizationFactor;
}

inline crossProduct(vector<double> &product, const vector<double> &U, const vector<double> &V) {
	if (U.size() != 3 or V.size() != 3) { cerr << "ERROR: size of U or V is not 3!"; abort(); }
	product.clear(); product.resize(3);
	product[0] = U[1] * V[2] - U[2] * V[1];
	product[1] = U[2] * V[0] - U[0] * V[2];
	product[0] = U[0] * V[1] - U[1] * V[0];
}

void generateSteepestDescentTranslationalVector(vector<double> &delT, vector<double> &F, double h, double alpha) {
	double factor = h * alpha / cblas_dnrm2(F.size(), &F[0], 1);
	delT.resize(F.size());
	for (unsigned int i=0; i < F.size(); i++) delT[i] = factor * F[i];
}

void generateMatrixWFromNormalizedVectorW(vector<double> &matW, const vector<double> &vecW) {
	if (vecW.size() != 3) { cerr << "ERROR: vecW size is not 3!"; abort(); }
	matW.resize(9); matW[0] = matW[4] = matW[8] = 0;
	matW[1] = vecW[2]; matW[2] = -vecW[1];
	matW[3] = -vecW[2]; matW[5] = vecW[0];
	matW[6] = vecW[1]; matW[7] = -vecW[0];
}

void generateSteepestDescentRotationMatrix(vector<double> &matR, const vector<double> &normalizedVectorT, double h, double beta) {
	if (vecW.size() != 3) { cerr << "ERROR: vecW size is not 3!"; abort(); }

	double theta = h * beta;
	double tsin = sin(theta), tcos = 1 - cos(theta);
	matR.clear(); matR.resize(9, 0); matR[0] = matR4] = matR[8] = 1;
	vector<double> matW;
	generateMatrixWFromNormalizedVectorW(matW, normalizedVectorT);

	for (unsigned int i=0; i < matR.size(); i++) matR[i] += tsin * matW[i];
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, tcos, &matW[0], 3, &matW[0], 3, 1, &matR[0], 3);
}

void calculateForceAndTorqueVectors(vector<double> &force, vector<double> &torque, vector<double> &coordA, vector<double> &coordB, vector<double> &comA) {
	force.clear(); force.resize(3, 0);
	torque.clear(); torque.resize(3, 0);
	vector<double> delSdelP, differenceVector(3), cProduct(3), PmP(3);

	for (unsigned int i=0; i < coordA.size(); i+=3) {
		delSdelP.clear(); delSdelP.resize(3, 0);	
		for (unsigned int j=0; j < coordB.size(); j+=3) {
			// not correct yet
			double expPmQ2 = -2 * exp( pow(cblas_dnrm2(differenceVector.size(), &differenceVector[0], 1), 2) );
			for (unsigned int k=0; k<3; k++) differenceVector[k] = coordA[i+k] - coordB[j+k];
			for (unsigned int k=0; k<3; k++) delSdelP[k] += expPmQ2 * differenceVector[k];
		}
		
		for (unsigned int k=0; k<3; k++) force[k] += delSdelP[k]; // add to force
		
		for (unsigned int k=0; k<3; k++) PmP[k] = coordA[i+k] - comA[k];
		crossProduct(cProduct, PmP, delSdelP);
		for (unsigned int k=0; k<3; k++) force[k] += cProduct[k]; // add to torque
	}
}
