require "Utils"

system "rm -rf random.pdb"
File.open("/Users/bm/roks/BosentanConformers/MABENSON/BEGIN_MD/cluster_v2/1_supimp/out.pdb").collect.split_by { |x| x =~ /MODEL/ }.shuffle.first(50).each do |chunk|	
	File.open("BOSENTAN_50CONFORMERS.pdb", "a+") { |f| chunk.each { |x| f.write x } }
end