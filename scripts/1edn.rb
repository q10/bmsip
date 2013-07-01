require 'Utils'




jjobs = Dir.globfiles("../1EDN/*").product( Dir.globfiles("../CONFORMERS/*") ).collect do |peptide, ligand|
	["../example", peptide, ligand, "../1EDN_SUPERIMPOSITIONS/"+peptide.basename+"__"+ligand.basename+".pdb"].join " "
end

runJobs(jjobs)
