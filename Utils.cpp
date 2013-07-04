void extractFileExtension(string &format, const string &fileName) {
    size_t npos = fileName.find_last_of(".");
    if(npos != std::string::npos) format = fileName.substr(npos+1);    
}

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

unsigned int importMoleculesFromFile(vector<OBMol> &moleculesList, const string &fileName) { // does not clear list, only appends new molecules to it
    string format; extractFileExtension(format, fileName);
    OBConversion obconversion;
    if (not obconversion.SetInFormat(format.c_str())) { cerr << "ERROR: OpenBabel does not recognize the following format: '" << format << "'; exiting" << endl; abort(); }

    unsigned int numMoleculesInFile = 0;
    OBMol mol; bool notAtEnd = obconversion.ReadFile(&mol, fileName.c_str());
    while (notAtEnd) {
        numMoleculesInFile++;
        moleculesList.push_back(mol);
        mol.Clear(); notAtEnd = obconversion.Read(&mol);
    }
    return numMoleculesInFile;
}

void writeMoleculeToFile(const string &fileName, OBMol &molecule, bool rewriteFile=false) {
    string format; extractFileExtension(format, fileName);

    ios_base::openmode fileMode = rewriteFile ? ios::out : (ios::out|ios::app);
    std::ofstream ofs(fileName.c_str(), fileMode);
    
    OBConversion obconversion;
    if (not obconversion.SetOutFormat(format.c_str())) { 
        cerr << "WARNING: OpenBabel does not recognize the following format: '" << format << "'; will write to SDF format" << endl;
        obconversion.SetOutFormat("sdf");
    }
    cout << "WRITING MOLECULE TO FILE '" << fileName << "'";
    if (rewriteFile) cout << " (WILL OVERWRITE EXISTING FILE IF ANY)...\n";
    else cout << "...\n";
    obconversion.Write(&molecule, &ofs);  // obconversion.WriteFile(&molecule, fileName.c_str());
}

// currently takes the first molecule in file, and if sequential molecules have different formulas, simply disposes of it
void importMoleculeConformersFromFile(vector<OBMol> &moleculesList, const string &fileName) {
    vector<OBMol> tempList;
    importMoleculesFromFile(tempList, fileName);
    OBMol *tempMoleculeBuild = new OBMol(tempList[0]);
    for (unsigned int i=1; i < tempList.size(); i++) {
        if (tempMoleculeBuild->GetFormula().compare( tempList[i].GetFormula() ) == 0) tempMoleculeBuild->AddConformer(tempList[i].GetCoordinates());
    }
    moleculesList.push_back(*tempMoleculeBuild);
}

void addConformerToMolecule(OBMol &molecule, vector<double> &coordinates) {
    double *copyInMemory = new double [coordinates.size()];
    memcpy(copyInMemory, &coordinates[0], sizeof(double) * coordinates.size());
    molecule.AddConformer(copyInMemory);
}

void tempImportMoleculesFromFile(vector<OBMol> &moleculesList, const string &fileName) {
    string format; extractFileExtension(format, fileName);
    OBConversion obconversion;
    if (not obconversion.SetInFormat(format.c_str())) { cerr << "ERROR: OpenBabel does not recognize the following format: '" << format << "'; exiting" << endl; abort(); }

    OBMol mol; bool notAtEnd = obconversion.ReadFile(&mol, fileName.c_str());
    vector<OBBase*> tempList;
    while (notAtEnd) {
        tempList.push_back(new OBMol(mol));
        mol.Clear(); notAtEnd = obconversion.Read(&mol);
    }

    //for (int i=0; i<tempPtr.size(); i++) cout << dynamic_cast<OBMol*>(tempPtr[i])->NumConformers() << endl;

    cout << tempList.size() << endl;

    OBOp* pOp = OBOp::FindType("readconformer");
    if(pOp)
        pOp->ProcessVec(tempList);


    cout << tempList.size() << endl;
    for (int i=0; i<tempList.size(); i++) cout << dynamic_cast<OBMol*>(tempList[i])->GetFormula() << endl;
}

void writeMoleculeConformersToFile(const string &fileName, OBMol &molecule, bool rewriteFile=false) {
    for (int i = molecule.NumConformers()-1; i >= 0; i--) { // DO NOT USE UNSIGNED INT i!!!
        molecule.SetConformer(i);
        writeMoleculeToFile(fileName, molecule, rewriteFile);
        rewriteFile = false; // if rewriteFile flag is on, keep it on only for the first conformer's write
    }
}

void printMoleculeCoords(OBMol &molecule) {
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++) cout << (*iter)->GetVector() << endl;
}

void generateCoordsMatrixFromMolecule(vector<double> &matrix, OBMol &molecule) {     // generates column-order matrix of coordinates
    matrix.clear(); matrix.insert(matrix.end(), molecule.GetCoordinates(), &molecule.GetCoordinates()[3*molecule.NumAtoms()]);
}

void generateCoordsMatrixFromMoleculeConformers(vector< vector<double> > &matrixes, OBMol &molecule) {     // generates column-order matrix of coordinates
    matrixes.resize( molecule.NumConformers() );
    for (unsigned int i=0; i < molecule.NumConformers(); i++) {
        molecule.SetConformer(i);
        generateCoordsMatrixFromMolecule(matrixes[i], molecule);
    }
    molecule.SetConformer(0);
}

void saveCoordsMatrixToMolecule(OBMol &molecule, vector<double> &matrix) {
    if (matrix.size() != molecule.NumAtoms() * 3) { cerr << "ERROR: INCORRECT MATCHING OF NUMBER OF COORDINATES; EXITING" << endl; abort(); }
    molecule.SetCoordinates(&matrix[0]);
}

void writeTemporaryMoleculeCoordsToFile(const string &fileName, OBMol &molecule, vector<double> &tempCoords, bool rewriteFile=false) {
    vector<double> oldCoords;
    generateCoordsMatrixFromMolecule(oldCoords, molecule);
    saveCoordsMatrixToMolecule(molecule, tempCoords);
    writeMoleculeToFile(fileName, molecule, rewriteFile);
    saveCoordsMatrixToMolecule(molecule, oldCoords);
}

void generateAtomicNumbersListFromMolecule(vector<int> &numList, OBMol &molecule) {
    numList.clear();
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++)
        numList.push_back((*iter)->GetAtomicNum());
}

void generateVDWRadiusListFromMolecule(vector<double> &VDWList, OBMol &molecule) {
    VDWList.clear();
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++)
        VDWList.push_back( etab.GetVdwRad((*iter)->GetAtomicNum()) );
}

void generateAtomicMassesListFromMolecule(vector<double> &massList, OBMol &molecule) {
    massList.clear();
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++)
        massList.push_back( (*iter)->GetAtomicMass() );
}

void generateAtomMatchScoringTableFromTwoMolecules(vector< vector<double> > &atomMatchScoringTable, OBMol &moleculeA, OBMol &moleculeB, double alpha=1.0, double beta=-1.0) {
    atomMatchScoringTable.resize(moleculeA.NumAtoms()); unsigned int i=0;
    for(unsigned int k=0; k < atomMatchScoringTable.size(); k++) atomMatchScoringTable[k].clear();

    for (OBAtomIterator iterA = moleculeA.BeginAtoms(); iterA != moleculeA.EndAtoms(); iterA++, i++) {
        OBAtom *atomA = *iterA; int a = 0;
        if (atomA->IsHbondDonor() and atomA->IsHbondAcceptor()) a = 6; // class B
        else if (atomA->IsHbondDonor()) a = 5; // class D
        else if (atomA->IsHbondAcceptor()) a = 4; // class A
        else if (atomA->IsAromatic()) a = 3; // class R
        else if (atomA->IsCarbon()) a = 2; // class L
        else if (not atomA->IsHydrogen()) a = 1; // class X

        for (OBAtomIterator iterB = moleculeB.BeginAtoms(); iterB != moleculeB.EndAtoms(); iterB++) {
            OBAtom *atomB = *iterB; int b = 0;
            if (atomB->IsHbondDonor() and atomB->IsHbondAcceptor()) b = 6; // class B
            else if (atomB->IsHbondDonor()) b = 5; // class D
            else if (atomB->IsHbondAcceptor()) b = 4; // class A
            else if (atomB->IsAromatic()) b = 3; // class R
            else if (atomB->IsCarbon()) b = 2; // class L
            else if (not atomB->IsHydrogen()) b = 1; // class X

            if (a == 0 or b == 0) atomMatchScoringTable[i].push_back(0); // ignore the hydrogens if they have not been removed yet
            else if (a == b) atomMatchScoringTable[i].push_back(alpha);
            else if ((a == 6 and (b == 5 or b == 4))  or  (b == 6 and (a == 5 or a == 4))) atomMatchScoringTable[i].push_back(alpha); // Scoring between HBond donors/acceptors 
            else { atomMatchScoringTable[i].push_back(beta); }
        }
    }
}

void getMoleculeCenterCoords(vector<double> &centerCoords, OBMol &molecule) {
    centerCoords.clear(); centerCoords.resize(3, 0);
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++) {
        double *tmpCoords = (*iter)->GetCoordinate(); double tmpMass = (*iter)->GetAtomicMass();  
        centerCoords[0] += tmpMass*tmpCoords[0]; centerCoords[1] += tmpMass*tmpCoords[1]; centerCoords[2] += tmpMass*tmpCoords[2];
    }
    for (unsigned int i=0; i < centerCoords.size(); i++) centerCoords[i] /= molecule.GetMolWt();
}

void getMoleculeConformerCenterCoords(vector< vector<double> > &centerCoords, OBMol &molecule) {
    centerCoords.resize( molecule.NumConformers() );
    for (unsigned int i=0; i < molecule.NumConformers(); i++) {
        molecule.SetConformer(i);
        getMoleculeCenterCoords(centerCoords[i], molecule);
    }
    molecule.SetConformer(0);
}

void printMoleculeCenterCoords(OBMol &molecule) {
    vector<double> centerCoords; getMoleculeCenterCoords(centerCoords, molecule);
    cout << "CENTER COORDS: ";
    for (unsigned int i=0; i < centerCoords.size(); i++) { cout << centerCoords[i] << " "; }
    cout << endl;
}

inline void translate3DMatrixCoordinates(vector<double> &matrix, double x, double y, double z) {
    for (unsigned int i=0; i < matrix.size(); i+=3) { matrix[i] += x; matrix[i+1] += y; matrix[i+2] += z; }
}

void rotate3DMatrixCoordinates(vector<double> &matrix, vector<double> &rotationMatrix) {
    vector<double> resultMatrix(matrix.size(), 0); // both matrices must be of column-order
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, matrix.size()/3, 3, 1, &rotationMatrix[0], 3, &matrix[0], 3, 0, &resultMatrix[0], 3);
    matrix = resultMatrix;
}

void removeNonBondedAtomsInMolecule(OBMol &molecule) {
    cout << "Deleting unbonded atoms... "; int numDeleted = 0;
    for (OBAtomIterator iter = molecule.BeginAtoms(); iter != molecule.EndAtoms(); iter++) {
        OBAtom *currentAtom = *iter;
        if ((not currentAtom->HasSingleBond()) and (not currentAtom->HasNonSingleBond()) and molecule.DeleteAtom(currentAtom)) {
            numDeleted++; iter--;
        }
    }
    cout << "deleted " << numDeleted << " unbonded atoms." << endl;
}

void generateConformers(OBMol &molecule, int numConformers=50, int numChildren=10, int mutability=5, int convergence=30) {
    OBConformerSearch cs;
    cs.Setup(molecule, numConformers, numChildren, mutability, convergence); // numConformers 30 // numChildren 5 // mutability 5 // convergence 25
    cs.SetScore(new OBEnergyConformerScore);
    cs.Search(); cs.GetConformers(molecule);

    //pmol->AddHydrogens(false, true); // Add some hydrogens before running MMFF

    OBForceField* pFF = OBForceField::FindForceField("MMFF94");
    if (not pFF) {
        pFF = OBForceField::FindForceField("UFF");
        if (not pFF) { cerr << "ERROR: CANNOT USE EITHER MMFF94 OR UFF FORCE FIELDS; EXITING\n"; exit(-1); }
        cerr << "WARNING: CANNOT USE MMFF94, USING UFF FORCE FIELD FOR CALCULATIONS\n";
    }
    
    for (int i = molecule.NumConformers()-1; i >= 0; i--) { // DO NOT USE UNSIGNED INT i!!!
        molecule.SetConformer(i);
        pFF->Setup(molecule);
        if (not pFF->Setup(molecule)) { cerr << "ERROR: CANNOT LOAD MOLECULE TO FORCEFIELDSOLVER; EXITING\n"; exit(-1); }
        //pFF->SteepestDescent(500, 1.0e-4);
        //pFF->WeightedRotorSearch(200, 25);
        pFF->ConjugateGradients(1000, 1.0e-6);
        pFF->GetCoordinates(molecule);
    }
}

double calculateRMSD(vector<double> &list1, vector<double> &list2, unsigned int dimensions=3) {
    if (list1.size() != list2.size() or list1.size() % dimensions != 0) { 
        cerr << "ERROR: Cannot run RMSD - list sizes don't match or are incorrect; list1 has size "
             << list1.size() << " while list2 has size " << list2.size() << "; exiting." << endl;
        abort(); 
    }
    double total = 0;
    for (unsigned int i=0; i < list1.size(); i++) {
        double difference = list2[i] - list1[i];
        total += difference * difference;
    }
    return sqrt(total * dimensions / list1.size());
}

double calculateRMSD(OBMol &moleculeA, OBMol &moleculeB) {
    if (moleculeA.NumAtoms() != moleculeB.NumAtoms()) { 
        cerr << "ERROR: Cannot run RMSD - number of atoms in molecule don't match or are incorrect; moleculeA has size "
             << moleculeA.NumAtoms() << " while moleculeB has size " << moleculeB.NumAtoms() << "; exiting." << endl;
        abort(); 
    }
    vector<double> matrix1, matrix2;
    generateCoordsMatrixFromMolecule(matrix1, moleculeA);
    generateCoordsMatrixFromMolecule(matrix2, moleculeB);
    return calculateRMSD(matrix1, matrix2);
}

void readPDBFile(vector< vector<string> > &PDBFile, const string &filename) {
    ifstream filestream(filename.c_str());
    if (not filestream.is_open()) { cerr << "Could not open input PDB file " << filename << "; exiting" << endl; abort(); }
    PDBFile.clear(); string line;
    while (getline(filestream, line)) {
        vector<string> tokenedPDBLine; string temp; istringstream iss(line);
        while (iss >> temp) { tokenedPDBLine.push_back(temp); }
        PDBFile.push_back(tokenedPDBLine);
    }
}

void writeMDPDBFile(const string &filename, vector< vector<string> > &PDBFile, vector< vector<double> > *coordSets) {
    ofstream outfile(filename.c_str());
    if (not outfile.is_open()) { cerr << "Could not open output PDB file " << filename << " for writing; exiting" << endl; abort(); }
    outfile << fixed << setprecision(3);

    int coordCounter = 0, setCounter = 0;
    for (unsigned int i=0; i < PDBFile.size(); i++) {
        vector<string> &line = PDBFile[i];
        if (line[0].compare("MODEL") == 0) {
            coordCounter = 0;
            outfile << line[0] << setw(16) << line[1] << "\n";
        } else if (line[0].compare("ENDMDL") == 0) {
            setCounter++;
            outfile << line[0] << "\n";
        } else {
            vector<double> &coords = (*coordSets)[setCounter];
            outfile << line[0] 
                       << right << setw(7) << line[1]
                       << " " << left << setw(5) << line[2]
                       << left << setw(4) << line[3]
                       << line[4]
                       << right << setw(4) << line[5]
                       << right << setw(12) << coords[coordCounter]
                       << right << setw(8) << coords[coordCounter+1]
                       << right << setw(8) << coords[coordCounter+2] << "\n";
                       coordCounter+=3;
        }
    }
    outfile.close();
}


void extractAllCoordsFromMDPDB(vector< vector<double> > **coordSetsHandle, const vector< vector<string> > &PDBFile) {
    *coordSetsHandle = new vector< vector<double> >();
    vector< vector<double> > &coordSets = **coordSetsHandle;
    for (unsigned int i=0; i < PDBFile.size(); i++) {
        if (PDBFile[i][0].compare("MODEL") == 0) coordSets.push_back( vector<double>() );
        else if ((PDBFile[i][0].compare("ATOM") == 0)) {
            coordSets.back().push_back( atof(&PDBFile[i][6][0]) );
            coordSets.back().push_back( atof(&PDBFile[i][7][0]) );
            coordSets.back().push_back( atof(&PDBFile[i][8][0]) );
        }
    }
}

void extractBackboneDataFromMDPDB(vector< vector<double> > **coordSetsHandle, vector< vector<double> > &comSets, vector<double> &VDWList, 
                                vector<double> &atomicMasses, const vector< vector<string> > &PDBFile) {
    /*
    ATOM      1  N   CYS A   1       1.943   5.540   6.216
    ATOM      2  H1  CYS A   1       1.256   5.445   6.951
    ATOM      3  H2  CYS A   1       2.854   5.693   6.624
    */
    *coordSetsHandle = new vector< vector<double> >();
    vector< vector<double> > &coordSets = **coordSetsHandle;
    VDWList.clear(); atomicMasses.clear(); bool flag = true;
    for (unsigned int i=0; i < PDBFile.size(); i++) {
        if (PDBFile[i][0].compare("MODEL") == 0) coordSets.push_back( vector<double>() );
        else if ((PDBFile[i][0].compare("ATOM") == 0) and (PDBFile[i][2].compare("N") == 0 or PDBFile[i][2].compare("CA") == 0 or 
                    PDBFile[i][2].compare("C") == 0 or PDBFile[i][2].compare("O") == 0)) {
            coordSets.back().push_back( atof(&PDBFile[i][6][0]) );
            coordSets.back().push_back( atof(&PDBFile[i][7][0]) );
            coordSets.back().push_back( atof(&PDBFile[i][8][0]) );


            if (flag) {
                if (PDBFile[i][2] == "N") {
                    VDWList.push_back(1.55); atomicMasses.push_back(14.0067); 
                } else if (PDBFile[i][2].compare("CA") == 0 or PDBFile[i][2].compare("C") == 0) {
                    VDWList.push_back(1.70); atomicMasses.push_back(12.0107);
                } else if (PDBFile[i][2] == "O") {
                    VDWList.push_back(1.52); atomicMasses.push_back(15.9994);
                }
            }
        }
        else if (PDBFile[i][0].compare("ENDMDL") == 0) flag = false;
    }

    comSets.clear();
    for (unsigned int i=0; i < coordSets.size(); i++) {
        double totalMass = 0; vector<double> centerCoords(3, 0);
        for (unsigned int j=0; j < coordSets[i].size();) {
            double tmpMass = atomicMasses[j]; totalMass += tmpMass;
            centerCoords[0] += tmpMass*coordSets[i][j++]; centerCoords[1] += tmpMass*coordSets[i][j++]; centerCoords[2] += tmpMass*coordSets[i][j++];
        }
        for (unsigned int i=0; i < centerCoords.size(); i++) centerCoords[i] /= totalMass;
        comSets.push_back(centerCoords);
    }
}

void generateCovarMatrixFromTables(vector<double> &matrix, const vector<double> &coords, const vector<double> &masses) {
    double uX = 0, uY = 0, uZ = 0, totalMass = 0; unsigned int j=0;
    for (unsigned int i=0; i < coords.size(); i+=3) {
        totalMass += masses[j++];
        uX += coords[i]; uY += coords[i+1]; uZ += coords[i+2];
    }
    uX /= masses.size(); uY /= masses.size(); uZ /= masses.size();

    double cXX = 0, cYY = 0, cZZ = 0, cXY = 0, cXZ = 0, cYZ = 0; j=0;
    for (unsigned int i=0; i < coords.size(); i+=3) {
        double atomWeight = masses[j++];
        cXX += pow(coords[i] - uX, 2) * atomWeight;
        cYY += pow(coords[i+1] - uY, 2) * atomWeight;
        cZZ += pow(coords[i+2] - uZ, 2) * atomWeight;
        cXY += (coords[i] - uX) * (coords[i+1] - uY) * atomWeight;
        cXZ += (coords[i] - uX) * (coords[i+2] - uZ) * atomWeight;
        cYZ += (coords[i+1] - uY) * (coords[i+2] - uZ) * atomWeight;
    }

    matrix.clear(); matrix.resize(9);
    matrix[0] = cXX; 
    matrix[1] = matrix[3] = cXY; 
    matrix[2] = matrix[6] = cXZ; 
    matrix[4] = cYY;
    matrix[5] = matrix[7] = cYZ;  
    matrix[8] = cZZ;

    for (unsigned int i=0; i < matrix.size(); i++) matrix[i] /= totalMass;
}
