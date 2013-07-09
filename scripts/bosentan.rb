require 'Utils'

system "cd .. && make clean all"


jjobs = Dir.globfiles("../CONFORMERS/*").delete_if { |x| x =~ /Bosentan/ }.collect do |ligand|
	reference = "BOSENTAN_50CONFORMERS.pdb"
	filename = "../BOSENTAN_SUPERIMPOSITIONS/"+reference.basename+"__"+ligand.basename
	["../example", reference, ligand, filename+".mol2", "&>", filename+".log"].join " "
end


#puts jjobs
runJobs(jjobs, 11, -20, 60)
system "cd .. && git add BOSENTAN_SUPERIMPOSITIONS/ && git commit -a -m \"bosentan superimposition\" && git push"
