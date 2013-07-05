require 'Utils'


=begin
jjobs = Dir.globfiles("../1EDN/*").product( Dir.globfiles("../CONFORMERS/*") ).collect do |peptide, ligand|
	filename = "../1EDN_SUPERIMPOSITIONS/"+peptide.basename+"__"+ligand.basename
	["../example", peptide, ligand, filename+".mol2", "&>", filename+".log"].join " "
end





runJobs(jjobs)
#puts jjobs


Dir.globfiles("../1EDN_SUPERIMPOSITIONS/*.pdb").each do |fl|
	system ["obabel", fl, "-O", "../1EDN_SUPERIMPOSITIONS/"+fl.basename+".sdf"].join " "
end


conformerhis = []
Dir.globfiles("../1EDN_SUPERIMPOSITIONS/PCA3_17-21*.log").each do |fl|
	tanimoto = open(fl).grep(/Tanimoto/)[0].gsub(/\n/,"").split(" ")[-1]
	conformers = open(fl).grep(/conformer\ A#[0-9]+\ and\ B#[0-9]+/)[0].split(" ").select { |w| w =~ /#/}.collect { |x| x.gsub(/[^a-zA-Z0-9]/, "").gsub(/[a-zA-Z]/, "") }.join "\t"
	puts [fl.basename, tanimoto, conformers].join "\t"
	conformerhis.push conformers.split("\t")[-1]
end
puts conformerhis.inspect
((0...50).collect{ |x| x.to_s }.to_a + conformerhis).histogram.sort { |x, y| x[0].to_i <=> y[0].to_i }.each { |x| puts x[0] + "\t" + (x[1]-1).to_s }


overlapContributions = []
ligands = []

%w(XRAY 3STEP PCA1 PCA2 PCA3).each do |typ|
	Dir.globfiles("../1EDN_SUPERIMPOSITIONS/"+typ+"*.mol2").each do |fl|
		ligands.push fl.basename.split("_")[-1]
		overlapContributions.push `../example #{fl}`.split("\n").collect { |x| x.to_f }.normalize.collect { |x| x.to_s }
	end
	overlapContributions.transpose.unshift(ligands).each { |x| puts x.join "," }
	overlapContributions = []
	ligands = []
	puts "\n"
end

=end

system "cd .. && make clean all"
jjobs = ["XRAY_17-21", "3STEPMIN_17-21", "PCA1_17-21", "PCA2_17-21", "PCA3_17-21"].product( Dir.globfiles("../CONFORMERS/*") ).collect do |peptide, ligand|
	original = "../1EDN/"+peptide+".pdb"
	cut =      "../1EDN/"+peptide+"_NOBACKBONE.pdb"
	filename = "../1EDN_SUPERIMPOSITIONS/ALPHA1_BETAN1_NOBACKBONE/"+peptide+"__"+ligand.basename
	["../example", cut, ligand, filename+".mol2", original, "&>", filename+".log"].join " "

end
runJobs(jjobs, 11, -20, 60)
#puts jjobs
