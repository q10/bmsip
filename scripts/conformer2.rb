require 'Utils'

#system "cd .. && make clean all"
=begin
jjobs = Dir.glob("../CONFORMERS/*").delete_if { |x| not File.file?(x) or x =~ /BQ123/ }.collect do |fl|
	ligand = File.basename( fl, ".*" ) # get base filename without the file extension
	reference = "../CONFORMERS/BQ123.mol2"
	final_struct = "../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/BQ123-" + ligand + ".mol2"
	logfile = "../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/BQ123-" + ligand + ".log"

	["../example", reference, fl, final_struct, "&>", logfile].join " "
end

runJobs(jjobs)


overlapContributions = []
ligands = []
Dir.globfiles("../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/*.mol2").each do |fl|
	ligands.push fl.basename.split("_")[-1]
	overlapContributions.push `../example #{fl}`.split("\n").collect { |x| x.to_f }.normalize.collect { |x| x.to_s }#.first(44)
end
overlapContributions.transpose.unshift(ligands).each { |x| puts x.join "\t" }

=end

files = []
File.open("../../MD/bensan/11_movie/movie2.pdb").each do |ln|
	if ln =~ /MODEL/
		files.push []
	elsif (not ln =~ /ENDMDL/)
		files.last.push ln
	end
end

files.each_with_index { |fl, i| File.open("MOVIES2/movie2_#{i}.pdb", 'w') { |n| n.write(fl.join) }  }