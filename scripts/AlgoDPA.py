import sys, itertools, math, operator, sqlite3, string
import scipy.stats

#################################
#######   AlgoDPA Class  ########
#################################
class AlgoDPA:
    def __init__(self, keys, numTraces, numBits, outFile=None):
        self.keys = keys
        self.numTraces = numTraces
        self.numBits = numBits
        self.possibleKeys = 2**numBits
        self.hammings = []
        self.simTraces = []
        self.peaks = {}
	self.conn = None
        self.cursor = None
    def init_database(self):
	self.conn = sqlite3.connect('../db/database')
        self.cursor = self.conn.cursor()
        self.cursor.execute('DROP TABLE IF EXISTS hammings')
        self.cursor.execute('CREATE TABLE hammings(ptext1 INTEGER, ptext2 INTEGER, key1 INTEGER, key2 INTEGER, round INTEGER, hd1 INTEGER, hd2 INTEGER, combined INTEGER);')
    def close_database(self):
        self.conn.commit()
        self.conn.close()
    def binseq(self, k):
        lst = [int(''.join(x),2) for x in itertools.product('01', repeat=k)]
        return lst
    def hammingDistance(self, d1, d2):
        return bin(d1 ^ d2)[2:].count('1')
    def generatePowerSimulationModel1(self, ptext):
        keys = self.binseq(self.numBits)
        data = ptext
        mask = 0b1
        d = {}
     
        for key in keys:
            d[key] = []
            tmpKey = key
            for i in xrange(self.numBits):
                previous_data = data
                previous_key = tmpKey
                data_lsb = mask & data
                key_lsb = mask & tmpKey
                xor = data_lsb ^ key_lsb
                data = data >> 1
                tmpKey = tmpKey >> 1
                xor = xor << self.numBits-1
                key_lsb = key_lsb << self.numBits-1
                data = data | xor
                tmpKey = tmpKey | key_lsb
                hd_data = self.hammingDistance(previous_data, data)
                d[key].append(hd_data)

            data = ptext

        self.hammings.append(d)
        
    def generatePowerSimulationModel2(self, ptext1, ptext2):
        keys = self.binseq(self.numBits*2)
        data1 = ptext1
        data2 = ptext2
        lsbMask = 0b1
        key1Mask = 0b111111000000 #Assuming numBits == 6.
        key2Mask = 0b000000111111
        d = {}

        for key in keys:
            d[key] = []
	    key1 = (key & key1Mask) >> 6
            tmpKey1 = key1
            key2 = (key & key2Mask) 
            tmpKey2 = key2
            for i in xrange(self.numBits):
                previous_data1 = data1
                previous_data2 = data2
                previous_key1 = tmpKey1
                previous_key2 = tmpKey2
                data1_lsb = lsbMask & data1
                data2_lsb = lsbMask & data2
                key1_lsb = lsbMask & tmpKey1
                key2_lsb = lsbMask & tmpKey2
                xor1 = data1_lsb ^ key1_lsb
                xor2 = data2_lsb ^ key2_lsb
                data1 = data1 >> 1
                data2 = data2 >> 1
                tmpKey1 = tmpKey1 >> 1
                tmpKey2 = tmpKey2 >> 1 
                xor1 = xor1 << self.numBits-1
                xor2 = xor2 << self.numBits-1
                key1_lsb = key1_lsb << self.numBits-1
                key2_lsb = key2_lsb << self.numBits-1
                data1 = data1 | xor1
                data2 = data2 | xor2
                tmpKey1 = tmpKey1 | key1_lsb
                tmpKey2 = tmpKey2 | key2_lsb
                hd_data1 = self.hammingDistance(previous_data1, data1)
                hd_data2 = self.hammingDistance(previous_data2, data2)
                d[key].append(hd_data1+hd_data2)
		if self.cursor != None:
                    values = (ptext1, ptext2, key1, key2, i, hd_data1, hd_data2, hd_data1+hd_data2) 
                    self.cursor.execute('INSERT INTO hammings VALUES(?, ?, ?, ?, ?, ?, ?, ?)', values)

            data1 = ptext1
            data2 = ptext2

        self.hammings.append(d)
        
    def simulateModel1(self, ptext):
        data = ptext
        key = self.keys
        mask = 0b1
        lst = []
        
        for i in xrange(self.numBits):
            previous_data = data
            previous_key = key
            data_lsb = mask & data
            key_lsb = mask & key
            xor = data_lsb ^ key_lsb
            data = data >> 1
            key = key >> 1
            xor = xor << self.numBits-1
            key_lsb = key_lsb << self.numBits-1
            data = data | xor
            key = key | key_lsb
            hd_data = self.hammingDistance(previous_data, data)
            lst.append(hd_data)

        self.simTraces.append(lst)

    def simulateModel2(self, ptext1, ptext2):
        data1 = ptext1
        data2 = ptext2
        key1 = self.keys[0]
        key2 = self.keys[1]
        mask = 0b1
        lst = []

        for i in xrange(self.numBits):
            previous_data1 = data1
            previous_data2 = data2
            previous_key1 = key1
            previous_key2 = key2
            data1_lsb = mask & data1
            data2_lsb = mask & data2
            key1_lsb = mask & key1
            key2_lsb = mask & key2
            xor1 = data1_lsb ^ key1_lsb
            xor2 = data2_lsb ^ key2_lsb
            data1 = data1 >> 1
            data2 = data2 >> 1
            key1 = key1 >> 1
            key2 = key2 >> 1
            xor1 = xor1 << self.numBits-1
            xor2 = xor2 << self.numBits-1
            key1_lsb = key1_lsb << self.numBits-1
            key2_lsb = key2_lsb << self.numBits-1
            data1 = data1 | xor1
            data2 = data2 | xor2
            key1 = key1 | key1_lsb
            key2 = key2 | key2_lsb
            hd_data1 = self.hammingDistance(previous_data1, data1)
            hd_data2 = self.hammingDistance(previous_data2, data2)
            lst.append(hd_data1+hd_data2)

        self.simTraces.append(lst)
        
        
    def attackModel1(self, measurements):
        xlist = []
        ylist = []
        correlations = {}

        for k in xrange(self.possibleKeys):
            for i in xrange(self.numTraces):
                for r in xrange(self.numBits):
                    xlist.append(self.hammings[i][k][r])
                    ylist.append(measurements[i][r])
                    
            corr = scipy.stats.pearsonr(xlist, ylist)
            correlations[k] = corr[0]
            del xlist[:], ylist[:]
            
        return correlations    

    def attackModel2(self, measurements):
        xlist = []
        ylist = []
        correlations = {}

        for k in xrange(2**(self.numBits*2)):
            for i in xrange(self.numTraces):
                for r in xrange(self.numBits):
                    xlist.append(self.hammings[i][k][r])
                    ylist.append(measurements[i][r])
                    
            corr = scipy.stats.pearsonr(xlist, ylist)
            correlations[k] = corr[0]
            del xlist[:], ylist[:]
            
        return correlations 

    def findPeaks(self, dataFile):
        start = 34
        end = 38
	samples = []
	trace = dataFile.readline()

        for i in xrange(self.numTraces):
            self.peaks[i] = []

            if i%100 == 0 and i!=0: print 'Finding peak values for '+str(i)+'th encryption...'

            for x in xrange(self.numBits):
                t2 = end*100
                t1 = start*100

		while trace != '':
                   t = trace.split()
                   if int(t[0]) >= t1 and int(t[0]) <= t2:
                       samples.append(Sample(int(t[0]), int(t[1])))
                   elif int(t[0]) > t2:
                       break
                   trace = dataFile.readline()

                peak = max(samples, key=operator.attrgetter('ampere'))
                
                self.peaks[i].append(peak.ampere)
                
                del samples[:]
                start = start + 10
                end = end + 10
                   
            start = start + 10
            end = end + 10

        return self.peaks  
    def findKey( self, correlations, hexFormat=False ):
        lst = sorted(correlations.iteritems(), key=operator.itemgetter(1), reverse=True)
        if hexFormat == True:
            formatted = []
            for item in lst:
                key1 = hex((item[0] & 0b111111000000) >> 6)[2:]
                key2 = hex((item[0] & 0b000000111111))[2:]
                key = key1+', '+key2
	        formatted.append((key, item[1]))
	    return formatted
        else:
            return lst
    def generateMeasurementsFile(self, outFile):
        fout = open(outFile, 'w+')
	for ptext in self.peaks:
            fout.write(','.join(str(roundPeak) for roundPeak in self.peaks[ptext])+'\n')
        fout.close()
       
#######################################################################

#################################
#######  Sample Class  ##########
#################################
class Sample:
    def __init__(self, time, ampere):
        self.time = time
        self.ampere = ampere
        self.correlationValue = 0
