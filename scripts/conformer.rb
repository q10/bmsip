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
=end


jjobs = Dir.glob("../CONFORMERS/*").delete_if { |x| x =~ /BQ123|234551|TAK044/ }.collect do |fl|
	ligand = File.basename( fl, ".*" ) # get base filename without the file extension
	reference = "../CONFORMERS/BQ123.mol2"
	final_struct = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".mol2"
	logfile = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".log"

	["../example", reference, fl, final_struct, "&>", logfile].join " "
end

runJobs(jjobs)

#system "git add ../ALL_PAIRS_BQ123_AS_REFERENCE/ && git commit -a -m \™results of all-pairs\™ && git push"
