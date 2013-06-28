require 'Utils'
conformerhis = []
all_tables = (0...50).collect.product((0...50).collect).transpose.reverse
names = []

Dir.globfiles("../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/*.log").sort {|x, y| x.basename.split(".")[-1] <=> y.basename.split(".")[-1]}.each do |fl|
	tanimoto = open(fl).grep(/Tanimoto/)[0].gsub(/\n/,"").split(" ")[-1]
	conformers = open(fl).grep(/conformer\ A#[0-9]+\ and\ B#[0-9]+/)[0].split(" ").select { |w| w =~ /#/}.collect { |x| x.gsub(/[^a-zA-Z0-9]/, "").gsub(/[a-zA-Z]/, "") }.join "\t"
	
  topMatchingBs = open(fl).grep(/ROUND/).collect { |x| ar=x.split(" "); [ar[-1], ar[4][2..-1]] }.sort {|y, x| x[0].to_f <=> y[0].to_f }.collect {|x| x[1] }.first(50)
  #puts topMatchingBs.inspect

  puts [fl.basename.split("-")[-1], tanimoto, conformers].join "\t"
#	conformerhis.push conformers.split("\t")[-1]
  conformerhis = conformerhis + topMatchingBs

#  names.push fl.basename.split("-")[-1]
#  all_tables.push open(fl).grep(/ROUND/).collect { |x| x.split(" ")[-1].to_f }
end

puts conformerhis.inspect
((0...50).collect{ |x| x.to_s }.to_a + conformerhis).histogram.sort { |x, y| x[0].to_i <=> y[0].to_i }.each { |x| puts x[0] + "\t" + (x[1]-1).to_s }

#all_tables = all_tables.transpose


#names.each_with_index do |x, i|
#  puts x + "\t" + all_tables.sort {|y, x| x[i+2] <=> y[i+2]}.first(20).collect{|x| x[1] }.join("\t")
#end

#Dir.globfiles("../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/*.log").delete_if {|x| x=~/BQ123-TAK044/ }.sort {|x, y| x.basename.split(".")[-1] <=> y.basename.split(".")[-1]}. each do
#all_tables.sort {|y, x| x[0] <=> y[0] }.each {|x| puts x[0]}

=begin
table =  open("../ALL_PAIRS_BQ123_AS_REFERENCE/ROKS2/BQ123-234551.log").grep(/ROUND/).collect do |x| 
  arr = x.split(" ")
  #a = arr[2][2..-1]
  #b = arr[4][2..-1]
  c = arr[-1]
  #[a, b, c]
c
end
=end
#table.transpose.each { |s| puts s.join "\t" }
#(0...50).collect.product((0...50).collect).each { |x, y| puts [y, x].join " " }
