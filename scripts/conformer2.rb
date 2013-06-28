require 'Utils'

system "cd .. && make clean all"

jjobs = Dir.glob("../CONFORMERS/*").delete_if { |x| not File.file?(x) or x =~ /BQ123/ }.collect do |fl|
	ligand = File.basename( fl, ".*" ) # get base filename without the file extension
	reference = "../CONFORMERS/BQ123.mol2"
	final_struct = "../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/BQ123-" + ligand + ".mol2"
	logfile = "../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/BQ123-" + ligand + ".log"

	["../example", reference, fl, final_struct, "&>", logfile].join " "
end

runJobs(jjobs)
