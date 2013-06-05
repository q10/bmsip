
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

def processorCount
	case RbConfig::CONFIG['host_os']
	when /darwin9/
		`hwprefs cpu_count`.to_i
    when /darwin/
		((`which hwprefs` != '') ? `hwprefs thread_count` : `sysctl -n hw.ncpu`).to_i
	when /linux/
		`cat /proc/cpuinfo | grep processor | wc -l`.to_i
	when /freebsd/
		`sysctl -n hw.ncpu`.to_i
	when /mswin|mingw/
		require 'win32ole'
		wmi = WIN32OLE.connect("winmgmts://")
		cpu = wmi.ExecQuery("select NumberOfCores from Win32_Processor") # TODO count hyper-threaded in this
		cpu.to_enum.first.NumberOfCores
	end
end



def runJobs(jobList, numThreads, priority=-20, sleepTime=60)
	raise ArgumentError unless (jobList.kind_of? Array or jobList.kind_of? String)
	raise ArgumentError unless (numThreads.kind_of? Integer and numThreads > 0)
	raise ArgumentError unless (numThreads.kind_of? Integer and numThreads >= -20 and numThreads < 20)
	raise ArgumentError unless ((sleepTime.kind_of? Integer or sleepTime.kind_of? Float) and  sleepTime > 0)
	jobList = [jobList] if jobList.kind_of? String

	pids = []
	while jobList.size > 0 do
		pids.delete_if { |x| Process.waitpid(x, Process::WNOHANG) }

		remaining = (numThreads - pids.size) > jobList.size ? jobList.size : (numThreads - pids.size) # in case number of open threads is larger than number of un-run trials
		remaining.times do
			text = jobList.pop
			pids.push fork { exec text.to_s }
			system "renice " + priority.to_s + " " + pids[-1].to_s
		end
		sleep sleepTime
	end
end


jjobs = Dir.glob("../CONFORMERS/*").delete_if { |x| x =~ /BQ123|234551|TAK044/ }.collect do |fl|
	ligand = File.basename( fl, ".*" ) # get base filename without the file extension
	reference = "../CONFORMERS/BQ123.mol2"
	final_struct = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".mol2"
	logfile = "../ALL_PAIRS_BQ123_AS_REFERENCE/BQ123_" + ligand + ".log"

	["../example", reference, fl, final_struct, "&>", logfile].join " "
end

runJobs(jjobs, 10)

system "git add ../ALL_PAIRS_BQ123_AS_REFERENCE/ && git commit -a -m \™results of all-pairs\™ && git push"
