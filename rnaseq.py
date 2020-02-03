#Main function, accepts file path, returns possible amino acid sequences
def main(file_path):
    
    dna_file = open(file_path, 'r')
    # Read the sequence
    dnaseq = ""
    dna_file.readline()
    for line in dna_file:
        line = line.rstrip()
        dnaseq = dnaseq + line
    
    
    temp = input("(T)emplate Strand or (N)on-Template Strand: ")
    temp.upper()
    orientation = input("(5)-3 or (3)-5: ")

	#if temp == 1 and orientation == 0:
		#There is no change, the sequence is what we desire
    if temp == "Y" and orientation == 5:
        dnaseq = dnaseq[::-1]
    elif temp == "N" and orientation == 3:
        dnaseq = compliment(dnaseq[::-1])
    elif temp == "N" and orientation == 5:
        dnaseq = compliment(dnaseq)

    #Transcribe DNA string to RNA string and then a sequence of amino acids
    aminac = read_ORF(dnaseq.upper())
    
    return aminac

#Returns compliment of dnaseq
def compliment(dnaseq):
    #Set DNA string to uppercase
    dnaseq = dnaseq.upper()
    compliment = list(dnaseq)
    for i in range(len(dnaseq)):
        if dnaseq[i] == "A":
            compliment[i] = "T"
        elif dnaseq[i] == "T":
            compliment[i] = "A"
        elif dnaseq[i] == "C":
            compliment[i] = "G"
        elif dnaseq[i] == "G":
            compliment[i] = "C"
    return "".join(compliment)

#Returns RNA sequence of standardized DNA strand
def transcribe(dnaseq):
    #Set DNA string to uppercase
    dnaseq = dnaseq.upper()
    rnaseq = list(dnaseq)
    for i in range(len(dnaseq)):
        if dnaseq[i] == "T":
            rnaseq[i] = "U"
            
         
    return "".join(rnaseq)

#Returns list of amino acids coded by RNA sequence
def translate(rnaseq):

    aminac_dict = {'UUU':'F','UUC':'F','UUA':'F','UUG':'F','UCU':'S',
                   'UCC':'S','UCA':'S','UCG':'S','UAU':'Y','UAC':'Y',
                   'UAA':'-','UAG':'-','UGU':'C','UGC':'C','UGA':'-',
                   'UGG':'M','CUU':'L','CUC':'L','CUA':'L','CUG':'L',
                   'CCU':'P','CCC':'P','CCA':'P','CCG':'P','CAU':'H',
                   'CAC':'H','CAA':'Q','CAG':'Q','CGU':'R','CGC':'R',
                   'CGA':'R','CGG':'R','AUU':'I','AUC':'I','AUA':'I',
                   'AUG':'M','ACU':'T','ACC':'T','ACA':'T','ACG':'T',
                   'AAU':'N','AAC':'N','AAA':'K','AAG':'K','AGU':'S',
                   'AGC':'S','AGA':'R','AGG':'R','GUU':'V','GUC':'V',
                   'GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A',
                   'GCG':'A','GAU':'D','GAC':'D','GAA':'E','GAG':'E',
                   'GGU':'G','GGC':'G','GGA':'G','GGG':'G'}
    
    aminac = []
    
    for i in range(0, len(rnaseq), 3):
        aminac.append(aminac_dict[rnaseq[i:i+3]])
    
    return "".join(aminac)

#Returns a list of possible starting points
def find_STARTS(dnaseq):
    starts = []
    for i in range(len(dnaseq) - 2):
            if dnaseq[i:i+3] == 'ATG':
                starts.append(i)
    return starts

#Find the next iteration of a given codon in a distinct reading frame(start)        
def find_next(dnaseq, start, codon):   
    for i in range(start, len(dnaseq), 3):
        if dnaseq[i:i+3] == codon:
            return i
    return -1

#Finds and returns the next stop codon in the current reading frame
def find_STOP(dnaseq, start):
    stop_codons = ["TAA", "TGA", "TAG"]
    stops = []
    for stop_codon in stop_codons:
        sc_position = find_next(dnaseq, start, stop_codon)
        if sc_position != -1:
            stops.append(sc_position)
        if len(stops) > 0:
            #Return the soonest stop codon relative to the start to ensure
            #similar reading frame
            return min(stops)
    return -1

#Finds and returns any possible reading frames for coding sequences in the form
#of tuples, (starting point, stopping point)   
def find_ORF(dnaseq):
    reading_frames = []
    start_positions = find_STARTS(dnaseq)
    for start in start_positions:
        stop = find_STOP(dnaseq, start)
        if stop != -1:
            reading_frames.append((start, stop))
    return reading_frames

#Reads and processes each frame
def read_ORF(dnaseq):
    coding_sequence = []
    for start, stop in find_ORF(dnaseq):
        ORF = dnaseq[start:stop+3]
        coding_sequence.append(translate(transcribe(''.join(ORF))))
    return coding_sequence
    
    
    