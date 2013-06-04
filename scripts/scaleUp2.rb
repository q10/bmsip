

=begin
(50..500).step(50) do |num|
	arrayList = []
	["1TMN", "3TMN", "4TMN", "1THL", "1TLP"].each do |name|
		input = `python rmsd/MolsRMSD.py RMSD_TEST/CHIMERA_SUPERPOSITION_1TMN_AS_REFERENCE/#{name}_ligand.sdf RMSD_TEST/MORE_CONFORMERS/#{name}conformer#{num}.sdf`
		processed = input.split("\t").find_all {|x| x.include? "\n"}.find_all {|x| x.gsub(/\n0/,"").match(/\A[+-]?\d+?(\.\d+)?\Z/)}.collect { |x| x.gsub(/\n0/,"").gsub(/\n/, "") }
		arrayList.push processed
	end
	output = File.open( "RMSD_TEST/combined_#{num}.rmsd", "w" )
	arrayList.transpose.each { |x| output << x.join("\t") + "\n" }
	output.close
end
=end


module Enumerable
    def sum
      self.inject(:+)
    end

    def mean
      self.sum / self.length.to_f
    end

    def sample_variance
      m = self.mean
      sum = self.inject(0){|accum, i| accum + (i-m)**2 }
      sum/(self.length - 1).to_f
    end

    def stddev
      return Math.sqrt(self.sample_variance)
    end

    def median
    	temp = self.sort
    	len = self.length
		return len % 2 == 1 ? temp[len/2] : (temp[len/2 - 1] + temp[len/2]) / 2
    end

end


def stats(temparr)
	arr = temparr.collect { |i| i.to_f }
	return [arr.mean, arr.stddev, arr.median, arr.min, arr.max]
end


["1TMN", "3TMN", "4TMN", "1THL", "1TLP"].each do |name|
	current_array = []

	(50..500).step(50) do |num|
		input = `python rmsd/MolsRMSD.py RMSD_TEST/CHIMERA_SUPERPOSITION_1TMN_AS_REFERENCE/#{name}_ligand.sdf RMSD_TEST/MORE_CONFORMERS/#{name}conformer#{num}.sdf`
		rmsd_values = input.split("\t").find_all {|x| x.include? "\n"}.find_all {|x| x.gsub(/\n0/,"").match(/\A[+-]?\d+?(\.\d+)?\Z/)}.collect { |x| x.gsub(/\n0/,"").gsub(/\n/, "") }
		current_array.push [num] + stats(rmsd_values)
	end

	output = File.open( "RMSD_TEST/combined_#{name}.rmsd_stats", "w" )
	current_array.each { |x| output << x.join("\t") + "\n" }
	output.close
end