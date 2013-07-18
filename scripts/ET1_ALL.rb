require 'Utils'


system "cd ../src && make clean all"

jjobs = Dir.globfiles("../CONFORMERS/*.mol2").delete_if {|x| x =~ /234551/ }.collect do |ligand|
	backbone = "../ET1/ET1_17-21_50CONFORMERS.pdb"
	peptide = "../ET1/ET1_17-21_50CONFORMERS_NOBACKBONE.pdb"
	filenamePrefix = "../ET1_SUPERIMPOSITIONS/A1BN1_ALL_NOBACKBONE/ET1__"+ligand.basename
	["../src/example", peptide, ligand, filenamePrefix, backbone, "&>", filenamePrefix+".log"].join " "
end

puts jjobs
#runJobs(jjobs.first(1), 11, -20, 60)
#system "cd .. && git add 1EDN_SUPERIMPOSITIONS/ && git commit -a -m \"no backbone, all comparison\" && git push"
