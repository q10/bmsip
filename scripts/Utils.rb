#require "RbConfig"

module Enumerable
	def to_f
		self.collect { |x| x.to_f }
	end
#	def to_s
#		self.collect { |x| x.to_s }
#	end
	def sum
		self.inject(:+)
	end
	def mean
		self.sum / self.length.to_f
	end
	def normalize
		s = self.sum
		self.collect{ |x| x / s }
	end
    def variance
		m = self.mean
		sum = self.inject(0){|accum, i| accum + (i-m)**2 }
		sum/(self.length - 1).to_f
    end
    def stddev
		return Math.sqrt(self.variance)
    end
    def median
		temp, len = self.sort, self.length
		return len % 2 == 1 ? temp[len/2] : (temp[len/2 - 1] + temp[len/2]) / 2
	end
    def dataStats
		arr = self.collect { |i| i.to_f }
		return [arr.mean, arr.stddev, arr.median, arr.min, arr.max]
	end
	def histogram
		his = inject(Hash.new(0)) { |h, x| h[x] += 1; h }
		return his.keys.zip( his.values )
  	end
	def split_by
		result = [a=[]]
		each { |o| yield(o) ? (result << a=[o]) : (a << o) }
		result.pop if a.empty?
		result.shift
		result
	end
end

class String
	def basename
		File.basename(self, ".*")
	end
end

class Dir
	def self.globfiles(path)
		self.glob(path).delete_if { |x| not File.file? x }
	end
end
=begin
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
=end
def self.runJobs(jobList, numThreads=11, priority=-20, sleepTime=60)
	raise ArgumentError unless (jobList.kind_of? Array or jobList.kind_of? String)
	raise ArgumentError unless (numThreads.kind_of? Integer and (numThreads+2) > 0)
	raise ArgumentError unless (priority.kind_of? Integer and (-20...20) === priority)
	raise ArgumentError unless ((sleepTime.kind_of? Integer or sleepTime.kind_of? Float) and  sleepTime > 0)

	jobList = [jobList] if jobList.kind_of? String
	numThreads = 1 if numThreads < 1
	pids = []
	startTime = Time.now
	while jobList.size > 0 do
		pids.delete_if { |x| Process.waitpid(x, Process::WNOHANG) }

		remaining = (numThreads - pids.size) > jobList.size ? jobList.size : (numThreads - pids.size) # in case number of open threads is larger than number of un-run trials
		remaining.times do
			text = jobList.pop
			pids.push fork { exec text.to_s }
			puts text
			system "renice " + priority.to_s + " " + pids[-1].to_s
		end
		sleep sleepTime
	end

	while pids.size > 0 do # finished submitting jobs in earlier loop; now waiting for jobs to complete
		pids.delete_if { |x| Process.waitpid(x, Process::WNOHANG) }
		sleep sleepTime
	end

	puts "Job batch finished in approximately #{Time.now - startTime} seconds."
end

=begin
egrep 'out of|SCF Don|Converged|NO|YES' Bosentan.log

# retains only protein heavy atoms, no waters, ions, or remarks
sed '/WAT/d' ../../7_min_PME/all_min3.pdb | sed '/CIP/d' | sed '/CIM/d' | sed '/REMARK/d' | awk '$3!~/^H/' > temp.pdb

# backbone atoms only
sed '/WAT/d' ../../7_min_PME/all_min3.pdb | sed '/CIP/d' | sed '/CIM/d' | sed '/REMARK/d' | awk '$3~/(^CA$)|(^C$)|(^N$)|(^O$)/' > temp.pdb

=end
