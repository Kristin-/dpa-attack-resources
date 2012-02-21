from AlgoDPA import *
import sys, sqlite3

if len(sys.argv) < 5:
    print 'Usage: <keys> <num_traces> <num_bits> <ptext_file> <traces_file>'
    exit()

keys = [int(key, 16) for key in sys.argv[1].split()]
numTraces = int(sys.argv[2])
numBits = int(sys.argv[3])
ptextFile = sys.argv[4]
tracesFile = sys.argv[5]

dpa = AlgoDPA(keys[0], numTraces, numBits)

print 'reading files..'
f = open(ptextFile, 'r')

for j, line in enumerate(f):
    if j == numTraces:
        print 'done with simulations...'
	break
    if j%100 == 0 and j!=0:  print 'Generating Simulation for '+str(j)+'th plain-text...'
    ptext = int(line, 2)
    dpa.generatePowerSimulationModel1(ptext)

f1 = open(tracesFile, 'r')

print 'finding peak values...'
peaks = dpa.findPeaks(f1)
print 'attacking...'
correlations = dpa.attackModel1(peaks)
results = dpa.findKey(correlations)

for item in results:
 print '('+hex(item[0])+', '+str(item[1])+')'

f.close()
f1.close()
print 'finished.'
