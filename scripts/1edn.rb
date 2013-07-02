require 'Utils'




jjobs = Dir.globfiles("../1EDN/*").product( Dir.globfiles("../CONFORMERS/*") ).collect do |peptide, ligand|
  filename = "../1EDN_SUPERIMPOSITIONS/"+peptide.basename+"__"+ligand.basename

	["../example", peptide, ligand, filename+".pdb", "&>", filename+".log"].join " "
end

runJobs(jjobs)
#puts jjobs
