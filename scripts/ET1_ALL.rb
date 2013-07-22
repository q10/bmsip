require 'Utils'


system "cd ../src && make clean all"


=begin
jjobs = Dir.globfiles("../CONFORMERS/*.mol2").collect do |ligand|
	backbone = "../ET1/ET1_17-21_50CONFORMERS.pdb"
	peptide = "../ET1/ET1_17-21_50CONFORMERS_NOBACKBONE.pdb"
	filenamePrefix = "../ET1_SUPERIMPOSITIONS/A1BN1_ALL_NOBACKBONE/ET1__"+ligand.basename
	["../src/example", peptide, ligand, filenamePrefix, backbone, "&>", filenamePrefix+".log"].join " "
end

#puts jjobs
runJobs(jjobs, 11, -20, 60)
system "cd .. && git add ET1_SUPERIMPOSITIONS/ && git commit -a -m \"no backbone, all comparison\" && git push"


overlapContributions = Dir.globfiles("../ET1_SUPERIMPOSITIONS/A1BN1_ALL_NOBACKBONE/*.mol2").collect { |fl| `../src/example #{fl}`.split("\n").to_f.normalize }
overlapContributions.transpose.each {|x| puts x.mean*100 }
=end


conformerhis = []
Dir.globfiles("../ET1_SUPERIMPOSITIONS/A1BN1_ALL_NOBACKBONE/*.log").each do |fl|
	tanimoto = open(fl).grep(/Tanimoto/)[0].gsub(/\n/,"").split(" ")[-1]
	conformers = open(fl).grep(/conformer\ A#[0-9]+\ and\ B#[0-9]+/)[0].split(" ").select { |w| w =~ /#/}.collect { |x| x.gsub(/[^a-zA-Z0-9]/, "").gsub(/[a-zA-Z]/, "") }.join "\t"
	puts [fl.basename, tanimoto, conformers].join "\t"
	conformerhis.push conformers.split("\t")[-1]
end
puts conformerhis.inspect
((0...50).collect{ |x| x.to_s }.to_a + conformerhis).histogram.sort { |x, y| x[0].to_i <=> y[0].to_i }.each { |x| puts x[0] + "\t" + (x[1]-1).to_s }
