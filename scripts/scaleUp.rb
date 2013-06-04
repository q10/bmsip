

def stats(arr)
	lowest = arr.min
	highest = arr.max
	len = arr.length
	mean = arr.inject(:+) / len # to_f so we don't get an integer result
	sorted = arr.sort
	median = len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
	return [mean, median, lowest, highest]
end





#input = `python rmsd/MolsRMSD.py RMSD_TEST/CHIMERA_SUPERPOSITION_1TMN_AS_REFERENCE/1TMN_ligand.sdf RMSD_TEST/ORIGINAL_PDBS/conformer1TMN.sdf`
#input.split("\t").find_all {|x| x.include? "\n"}.find_all {|x| x.gsub(/\n0/,"").match(/\A[+-]?\d+?(\.\d+)?\Z/)}.collect { |x| x.gsub(/\n0/,"").gsub(/\n/, "").to_f }

system "make clean all"

["1TMN", "3TMN", "4TMN", "1THL", "1TLP"].each do |x|
	orig = "RMSD_TEST/ORIGINAL_PDBS/ligand" + x + ".sdf"
	conformerout = "RMSD_TEST/MORE_CONFORMERS/" + x + "conformer"
	system "./example " + orig + " " + conformerout
end
