require "Utils"

system "rm -rf random.pdb"
File.open("../BosentanConformers/MABENSON/BEGIN_MD/cluster_v2/1_supimp/out.pdb").collect.split_by { |x| x =~ /MODEL/ }.shuffle.first(200).each do |chunk|	
	File.open("random.pdb", "a+") { |f| chunk.each { |x| f.write x } }
end