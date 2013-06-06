require "RbConfig"

def self.processorCount
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

def runJobs(jobList, numThreads=nil, priority=-20, sleepTime=60)
	raise ArgumentError unless (jobList.kind_of? Array or jobList.kind_of? String)
	raise ArgumentError unless (numThreads.nil? or (numThreads.kind_of? Integer and numThreads > 0))
	raise ArgumentError unless (priority.kind_of? Integer and (-20...20) === priority)
	raise ArgumentError unless ((sleepTime.kind_of? Integer or sleepTime.kind_of? Float) and  sleepTime > 0)
	jobList = [jobList] if jobList.kind_of? String
	if numThreads.nil?
		ideal = processorCount() - 2
		numThreads = ideal if ideal > 0 else 1 
	end

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
