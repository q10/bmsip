import os
from multiprocessing import Pool

def filesIn(mypath):
    return [os.path.join(mypath, f) for f in os.listdir(mypath)]

def fileNameWOExtension(path):
	return os.path.split(path)[-1].split('.')[0]


#os.system('make clean all')

files = filesIn("./CONFORMERS");
#for i in range( len(files) - 1 ):
#	for j in range (i+1, len(files)):
#		combined = fileNameWOExtension(files[i]) + '_' + fileNameWOExtension(files[j])
#		print "./example", files[i], files[j], combined + '.sdf', '&>', combined + '.log'


filesList = []
for i in range( len(files) - 1 ):
	for j in range (i+1, len(files)):
		filesList.append( (files[i], files[j]) )


def runProgram(files):
	combined = fileNameWOExtension(files[0]) + '_' + fileNameWOExtension(files[1])
	print "./example", files[0], files[1], combined + '.sdf', '&>', combined + '.log'

pool = Pool(processes=4)              # start 4 worker processes
pool.map(runProgram, filesList)
