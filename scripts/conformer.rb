require 'Utils'

=begin
Dir.glob("PRELIM2DTO3D/*") do |f|
	outfile = "CONFORMERS/" + f.split("/").last.split(".").first + ".mol2"
	system "./example " + f + " " + outfile
end

%w(Ambrisentan.sdf Atrasentan.sdf Avosentan.sdf BMS193884.sdf Bosentan.sdf BQ123.sdf BQ788.sdf Clazosentan.mol Darusentan.sdf Edonentan.sdf Enrasentan.sdf FR139317.mol2 J104132.sdf Macitentan.sdf Nebentan.sdf SB209670.sdf TAK044.sdf TBC3711.sdf TezosentanDisodium.sdf Zibotentan.sdf).each do |f|
	outfile = "CONFORMERS/" + f.split(".").first + ".mol2"
	system "./example ../Downloads/ANALOGS/FORMAL/" + f + " " + outfile
end


jjobs = Dir.glob("../CONFORMERS/*").delete_if { |x| not File.file?(x) or x =~ /BQ123|234551|TAK044/ }.collect do |fl|
	ligand = File.basename( fl, ".*" ) # get base filename without the file extension
	reference = "../CONFORMERS/BQ123.mol2"
	final_struct = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".mol2"
	logfile = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".log"

	["../example", reference, fl, final_struct, "&>", logfile].join " "
end

runJobs(jjobs)

#system "git add ../ALL_PAIRS_BQ123_AS_REFERENCE/ && git commit -a -m \™results of all-pairs\™ && git push"
=end


# Grab Tanimoto Indexes and conformers
=begin
conformerhis = []
Dir.globfiles("../ALL_PAIRS_BQ123_AS_REFERENCE/*.log").each do |fl|
	tanimoto = open(fl).grep(/Tanimoto/)[0].gsub(/\n/,"").split(" ")[-1]
	conformers = open(fl).grep(/conformer\ A#[0-9]+\ and\ B#[0-9]+/)[0].split(" ").select { |w| w =~ /#/}.collect { |x| x.gsub(/[^a-zA-Z0-9]/, "").gsub(/[a-zA-Z]/, "") }.join "\t"
	puts [fl.basename, tanimoto, conformers].join "\t"
	conformerhis.push conformers.split("\t")[-1]
end
puts conformerhis.inspect
((0...50).collect{ |x| x.to_s }.to_a + conformerhis).histogram.sort { |x, y| x[0].to_i <=> y[0].to_i }.each { |x| puts x[0] + "\t" + (x[1]-1).to_s }


overlapContributions = []
ligands = []
Dir.globfiles("../ALL_PAIRS_BQ123_AS_REFERENCE/*.mol2").each do |fl|
	ligands.push fl.basename.split("_")[-1]
	overlapContributions.push `../example #{fl}`.split("\n").collect { |x| x.to_f }.normalize.collect { |x| x.to_s }
end
overlapContributions.transpose.unshift(ligands).each { |x| puts x.join "\t" }

a = Dir.globfiles("../../Downloads/ANALOGS/2D/*").delete_if { |x| x=~/mat.dat/ }.collect #do |fl|
	#puts fl.basename #['obabel', fl, '-O', fl.basename+'.png'].join " "
#end

puts a.each_slice(6).each { |x| puts x.collect {|y| y.basename}.join "\t" }
=end


system "cd .. && make clean all"

files = %w(1TMN 3TMN 4TMN 1THL 1TLP)
jjobs = files.product(files).collect do |x, y|
	targetBeginningPosition = "../RMSD_TEST/ORIGINAL_PDBS/ligand" + y + ".mol2"
	targetBeginningPositionConformers = "../RMSD_TEST/ORIGINAL_PDBS/conformer" + y + ".mol2"
	reference = "../RMSD_TEST/CHIMERA_SUPERPOSITION_1TMN_AS_REFERENCE/" + x + "_ligand.mol2"
	targetXRayMatch = "../RMSD_TEST/CHIMERA_SUPERPOSITION_1TMN_AS_REFERENCE/"+ y + "_ligand.mol2"
	originalTargetOuput = "../RMSD_TEST/ROKS2/" + x + "-" + y + "_original.sdf"
	conformerTargetOutput = "../RMSD_TEST/ROKS2/" + x + "-" + y + "_conformers.sdf"
	logfile = "../RMSD_TEST/ROKS2/" + x + "-" + y + ".log"
	["./example", targetBeginningPosition, targetBeginningPositionConformers, 
		reference, targetXRayMatch, originalTargetOuput, conformerTargetOutput, ">", logfile, "2>&1"].join " "
end
puts jjobs
#runJobs(jjobs)