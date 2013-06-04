
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

THREADS = 10
jjobs = []
pids = []

Dir.glob("../CONFORMERS/*").delete_if { |x| x =~ /BQ123|234551|TAK044/ }.each do |fl|
	ligand = File.basename( fl, ".*" ) # get base filename without the file extension
	reference = "../CONFORMERS/BQ123.mol2"
	final_struct = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".mol2"
	logfile = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".log"

	jjobs.push ["../example", reference, fl, final_struct, "&>", logfile].join " "
end

THREADS.times do
	pids.push fork { exec jjobs.pop }
	system "renice -20 " + pids[-1].to_s
end

while jjobs.size > 0 do 
	pids.delete_if { |x| Process.waitpid(x, Process::WNOHANG) }

	remaining = (THREADS - pids.size) > jjobs.size ? jjobs.size : (THREADS - pids.size) # in case number of open threads is larger than number of un-run trials
	remaining.times do
		pids.push fork { exec jjobs.pop }
		system "renice -20 " + pids[-1].to_s
	end
	sleep 60
end

system "git add ../ALL_PAIRS_BQ123_AS_REFERENCE/ && git commit -a -m \™results of all-pairs\™ && git push"