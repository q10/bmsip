system "make clean all"

SHAEP = "../Downloads/shaep"
files = ["1TMN", "3TMN", "4TMN", "1THL", "1TLP"]
files.product(files).each do |x,y|
	reference = "RMSD_TEST/1TMN_AS_REFERENCE_IN_CHIMERA_SUPERPOSITION/" + x + "_ligand.mol2"
	outputfile = "RMSD_TEST/SHAEP_SUPERPOSITION/out" + x + "ref_" + y + ".sdf"
	logfile = "RMSD_TEST/log" + x + "_" + y + ".txt"
	targetconformers = "RMSD_TEST/ORIGINAL_PDBS/conformer" + y + ".mol2"
	system SHAEP + " -q " + reference + " " + targetconformers + " -s " + outputfile + " " + logfile

	xraytarget = "RMSD_TEST/1TMN_AS_REFERENCE_IN_CHIMERA_SUPERPOSITION/"+ y + "_ligand.sdf"
	system "./example " + xraytarget + " " + outputfile + " >>SHAEP_NUMBERS"
end