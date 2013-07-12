require 'Utils'
=begin
system "cd .. && make clean all"


jjobs = Dir.globfiles("../CONFORMERS/*").delete_if { |x| x =~ /Bosentan/ }.collect do |ligand|
	reference = "BOSENTAN_50CONFORMERS.pdb"
	filename = "../BOSENTAN_SUPERIMPOSITIONS/"+reference.basename+"__"+ligand.basename
	["../example", reference, ligand, filename+".mol2", "&>", filename+".log"].join " "
end


#puts jjobs
runJobs(jjobs, 11, -20, 60)
system "cd .. && git add BOSENTAN_SUPERIMPOSITIONS/ && git commit -a -m \"bosentan superimposition\" && git push"


Dir.globfiles("../BOSENTAN_SUPERIMPOSITIONS/*.log").each do |filename|
	ligand = filename.basename.split("__")[-1]
	tanimoto = open(filename).grep(/Tanimoto/)[0].gsub(/\n/,"").split(" ")[-1]
	puts [ligand, tanimoto].join "\t"
end

=end



overlapContributions = []
ligands = []
Dir.globfiles("../BOSENTAN_SUPERIMPOSITIONS/*_0_*.mol2").each do |fl|
	ligands.push fl.basename.split("__")[-1].split(".")[0]
	overlapContributions.push `../example #{fl}`.split("\n").collect { |x| x.to_f }.normalize.collect { |x| x.to_s }
end
overlapContributions.transpose.unshift(ligands).each { |x| puts x.join "\t" }
