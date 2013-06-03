system "make clean all"

[["3TMN", "1TMN"], 
["4TMN", "1TMN"], 
["4TMN", "3TMN"], 
["1THL", "1TMN"], 
["1THL", "3TMN"], 
["1THL", "4TMN"], 
["1TLP", "1TMN"], 
["1TLP", "3TMN"], 
["1TLP", "4TMN"], 
["1TLP", "1THL"]].each do |p|
	x, y = p[0..1]
	targetorig = "RMSD_TEST/ORIGINAL_PDBS/ligand" + y + ".sdf"
	targetconformer = "RMSD_TEST/ORIGINAL_PDBS/conformer" + y + ".sdf"

	reference = "RMSD_TEST/1TMN_AS_REFERENCE_IN_CHIMERA_SUPERPOSITION/" + x + "_ligand.sdf"
	chimeratarget = "RMSD_TEST/1TMN_AS_REFERENCE_IN_CHIMERA_SUPERPOSITION/" + y + "_ligand.sdf"

	outputxray = "RMSD_TEST/ROKS_SUPERPOSITION/" + y + "_" + x + "reference.sdf"
	outputconformer = "RMSD_TEST/ROKS_SUPERPOSITION/conformer" + y + "_" + x + "reference.sdf"

	system ["./example", targetorig, targetconformer, reference, chimeratarget, outputxray, outputconformer, ">>RMSD_OUTPUT3 2>&1"].join " "
end

system "mv RMSD_OUTPUT3 RMSD_TEST/RMSD_OUTPUT3 && git add RMSD_TEST/RMSD_OUTPUT3 && git add RMSD_TEST/ROKS_SUPERPOSITION/"
system "git commit -a -m \"completing the rmsd table\" && git push"