system "make clean all"

SHAEP = "../Downloads/shaep"
files = ["1TMN", "3TMN", "4TMN", "1THL", "1TLP"]
files.product(files).each do |x,y|
	reference = "RMSD_TEST/CHIMERA_SUPERPOSITION_1TMN_AS_REFERENCE/" + x + "_ligand.mol2"
	outputfile = "RMSD_TEST/SHAEP_SUPERPOSITION/out" + x + "ref_" + y + "_query.sdf"
	logfile = "RMSD_TEST/log" + x + "_" + y + ".txt"
	targetconformers = "RMSD_TEST/ORIGINAL_PDBS/conformer" + y + ".mol2"
	system SHAEP + " --onlyshape -q " + reference + " " + targetconformers + " -s " + outputfile + " " + logfile

	xraytarget = "RMSD_TEST/CHIMERA_SUPERPOSITION_1TMN_AS_REFERENCE/"+ y + "_ligand.sdf"
	system "./example " + xraytarget + " " + outputfile + " >>SHAEP_NUMBERS"
end