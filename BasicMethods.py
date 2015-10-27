RNATable = {
			'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V', 'UUC':'F', 'UUA':'L',
            'UUG':'L', 'UCU':'S', 'UCC':'S', 'CUC':'L', 'AUC':'I', 'GUC':'V',
            'CUA':'L', 'AUA':'I', 'AUG':'M', 'GUG':'V', 'CCU':'P', 'ACU':'T',
            'GCU':'A', 'CCC':'P', 'ACC':'T', 'GCC':'A', 'UCA':'S', 'CCA':'P',
            'ACA':'T', 'GCA':'A', 'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A',
            'UAU':'Y', 'CAU':'H', 'AAU':'N', 'GAU':'D', 'UAC':'Y', 'CAC':'H',
            'AAC':'N', 'GAC':'D', 'CAA':'Q', 'AAA':'K', 'GAA':'E', 'CAG':'Q',
            'AAG':'K', 'GAG':'E', 'UGU':'C', 'CGU':'R', 'AGU':'S', 'GGU':'G',
            'GGC':'G', 'CGA':'R', 'AGA':'R', 'GGA':'G', 'UGG':'W', 'UGC':'C',
            'CGG':'R', 'AGG':'R', 'GGG':'G', 'CUG':'L', 'GUA':'V', 'AGC':'S',
			'CGC':'R', 'UGA':'Stop', 'UAG':'Stop', 'UAA':'Stop'
}
#for convenience

def FASTAstr(strIn):
#takes in a list of strings in FASTA format and outputs a dictionary with the label as the key
    DNAstrings = {}
    for line in strIn:
        if '>' in line:
            label = line.lstrip('>')
            finalLabel = label.rstrip('\n')
            DNAstrings[finalLabel] = ''
        if '>' not in line:
            DNAstrings[finalLabel] = DNAstrings[finalLabel]+line.rstrip('\n')
    return DNAstrings


def gcContent(sequence):
# calculates the GC content of a given DNA sequence, returns as a percentage
    content = 0
    for letter in sequence:
        if letter is 'G' or letter is 'C':
            content+=1
    contentPer = content / len(sequence) * 100
    return contentPer
	
def profileSeqs(FASTAin):
#takes in a gene sequence and outputs a profile matrix
    profileMatrix = {}
    profileMatrix['A'] = []
    profileMatrix['T'] = []
    profileMatrix['G'] = []
    profileMatrix['C'] = []
    for word in FASTAin:
        while len(profileMatrix['A']) != len(FASTAin[word]):
            profileMatrix['A'].append(0)
            profileMatrix['T'].append(0)
            profileMatrix['G'].append(0)
            profileMatrix['C'].append(0)
        for i in range(0,len(FASTAin[word])):
            profileMatrix[FASTAin[word][i]][i] +=1
    return profileMatrix	

def CanonSeq(profile):
# returns a canonical sequence for a given profile matrix
    canonString = ''
    canonTracker = []
    for i in range(0,len(profile['A'])):
        canonTracker.append(['A',0])
    for word in profile:
        for x in range(0,len(profile[word])):
            if canonTracker[x][1] < profile[word][x]:
                canonTracker[x] = [word,profile[word][x]]
    for y in range(0,len(canonTracker)):
        canonString = canonString+canonTracker[y][0]
    return canonString

def makeAgencyList(INseqs):
#takes in a dictionary of sequences with names and generates an agency list
    agencyList = []
    for word in INseqs:
        for word2 in INseqs:
            if word != word2:
                if INseqs[word][-3:] == INseqs[word2][0:3]:
                    agencyList.append(word+ ' '+ word2)
    return agencyList

def findShortest(stringset):
#returns the shortest string in a dict of strings
    currentLength = 0
    for word in stringset:
        if currentLength > len(stringset[word]) or currentLength == 0:
            shortest = word
            currentLength = len(stringset[word])
    return shortest

	
def longestSub(stringset):
#calculates the longest common substring of a given dict of strings
    smallest = stringset[findShortest(stringset)]
    maxlen = len(smallest)
    for x in range(0,maxlen):
        sublen = maxlen-x
        for i in range(0,(maxlen-sublen)):
            testStr = smallest[i:sublen+i]
            for word in stringset:
                if testStr not in stringset[word]:
                    testStr = ''
            currentSoln = testStr
            if currentSoln != '':
                return currentSoln

def baseMatch(nucAcid,base):
#returns the complementary base for a given base
    if base.upper() == 'A':
        if nucAcid.lower() == 'rna':
                return 'U'
        elif nucAcid.lower() == 'dna':
                return 'T'
        else:
            return "Please specifcy the correct nucleic acid type"
    elif base.upper() == 'T' or base.upper() == 'U':
        return 'A'
    elif base.upper() == 'G':
        return 'C'
    elif base.upper() == 'C':
        return 'G'
    else:
        return 'Not a standard base'


def compSequence(inputSequence):
#returns the reverse complement of an input sequence
    flipSeq = inputSequence[::-1]
    compSeq = ''
    if 'U' in flipSeq:
        nucAcid = 'RNA'
    else:
        nucAcid = 'DNA'
    for letter in flipSeq:
        compSeq = compSeq+baseMatch(nucAcid,letter)
    return compSeq

def RNAtoProtein(RNA):
#takes an input RNA string and outputs the protein built from it
#also outputs whether or not translation terminated at a "stop" codon or not
    stopped = False
    position = 0
    currentCodon = ''
    proSeq = ''
    for x in range(0,len(RNA)):
        if position == 0:
            position = 1
            currentCodon = currentCodon+RNA[x]
        elif position == 1:
            position = 2
            currentCodon = currentCodon+RNA[x]
        elif position == 2:
            position = 0
            currentCodon = currentCodon+RNA[x]
            if RNATable[currentCodon]=='Stop':
                stopped = True
                break
            else:
                proSeq = proSeq + RNATable[currentCodon]
                currentCodon = ''
    return proSeq, stopped
	
def readingFrames(seqs):
#takes an input RNA sequence and returns the protein sequences based on the open reading frames
#checks the "Stopped" variable output by RNAtoProtein to verify a full open reading frame
    openFrameResults = []
    for sequ in seqs:
        frames = sequ.split('AUG')
        for x in range(1,len(frames)):
            frames[x]='AUG'+frames[x]
        for x in range(1,len(frames)):
            for y in range((x+1),len(frames)):
                frames[x]=frames[x]+frames[y]
        frames.remove(frames[0])
        for frame in frames:
            translated = RNAtoProtein(frame)
            if True in translated:
                openFrameResults.append(translated[0])
    return set(openFrameResults)