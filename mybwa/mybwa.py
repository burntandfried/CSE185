#!/usr/bin/env python
# coding: utf-8

# In[13]:

import sys
if len(sys.argv) == 3:
    genomeFile = sys.argv[1]
    if '.fa' not in genomeFile and '.fasta' not in genomeFile:
        print('Reference genome must be a fasta file')
        sys.exit(1)
        
    readsFile = sys.argv[2]
    if '.fq' not in readsFile and '.fastq' not in readsFile:
        print('Reads must be in a fastq file')
        sys.exit(1)
else:
    print('Unable to run, please use the following format:')
    print('python myBwa.py genome.fa reads.fq')
    sys.exit(1)


with open(genomeFile, 'r') as file:
    genome = ''
    file.readline()
    for line in file:
        line = line.strip()
        genome = genome + line

if len(genome) == 0:
    print('Reference genome file is empty')
    sys.exit(1)
    
# In[16]:


# returns all cyclic rotations of text
def CyclicRotations(text):
    rotations = [text]
    
    curr = text[-1] + text[:len(text)-1]
    while curr != text:
        rotations.append(curr)
        curr = curr[-1] + curr[:len(text)-1]
        
    return rotations


# In[17]:


# returns the burrows-wheeler transform of text
def BWT(text):
    rotations = CyclicRotations(text)
    sortedRotations = sorted(rotations)
    
    bwt = ''
    for string in sortedRotations:
        bwt = bwt + string[-1]
    
    return bwt


# In[19]:


bwt = BWT(genome)


# In[51]:


with open(readsFile, 'r') as file:
    counter = 0
    readDict = {}
    for line in file:
        if counter < 2:
            if counter == 0:
                key = line.strip()
            elif counter == 1:
                value = line.strip()   
                readDict[key] = value
            counter = counter + 1

        else:
            counter = 0
            file.readline()

    if counter == 1:
        readDict[key] = value


# In[ ]:

if len(readDict) == 0:
    print('Reads file is empty')
    sys.exit(1)



# In[52]:


# returns all suffixes of a string
def GetSuffixes(text):
    suffixes = []
    for i in range(len(text)):
        suffixes.append(text[i:])
    return suffixes


# In[53]:


# returns the suffix array of text (sort all suffixes alphabetically, assuming $ is at beginning of alphabet)
# suffix array is the array of start indices of the suffixes in their alphabetical order
def SuffixArray(text):
    suffixes = GetSuffixes(text)
    suffixes = sorted(suffixes)
    
    suffArray = []
    for suff in suffixes:
        suffArray.append(len(text) - len(suff))
        
    return suffArray


# In[54]:


# returns a dictionary, key = sym, value is a list of the number of times the symbol occurs in the first i positions of lastColumn
def CountDict(lastColumn):
    count = {sym: [0 for i in range(len(lastColumn)+1)] for sym in lastColumn}

    for i in range(len(lastColumn)):
        sym = lastColumn[i]
        count[sym][i+1] = count[sym][i] + 1
        for j in range(i+1, len(count[sym])):
            count[sym][j] = count[sym][i+1]
        
    return count


# In[55]:


# returns a dictionary, key = symbol and value = the first position that symbol occurs in the first column of the matrix
def FirstOccurenceDict(lastColumn):
    firstOccur = {}
    firstCol = sorted(lastColumn)
    
    for i in range(len(firstCol)):
        if firstCol[i] not in firstOccur:
            firstOccur[firstCol[i]] = i
    
    return firstOccur


# In[56]:


# returns all cyclic rotations of text
def CyclicRotations(text):
    rotations = [text]
    
    curr = text[-1] + text[:len(text)-1]
    while curr != text:
        rotations.append(curr)
        curr = curr[-1] + curr[:len(text)-1]
        
    return rotations


# In[57]:


# returns the burrows-wheeler transform of text
def BWT(text):
    rotations = CyclicRotations(text)
    sortedRotations = sorted(rotations)
    
    bwt = ''
    for string in sortedRotations:
        bwt = bwt + string[-1]
    
    return bwt


# In[58]:


# faster BWTMatching implementation: returns the start indicies of the occurences of a pattern in a string using its BWT
def BetterBWTMatching(firstOccur, lastColumn, pattern, count, suffArray):
    top = 0
    bottom = len(lastColumn) - 1
    
    patternOccurences = []
    while top <= bottom:
        if pattern:
            symbol = pattern[-1] # store the last letter
            pattern = pattern[:-1] # remove the last letter
            if symbol in lastColumn[top:bottom+1]:
                top = firstOccur[symbol] + count[symbol][top]
                bottom = firstOccur[symbol] + count[symbol][bottom+1] - 1
            
            # pattern does not occur in text
            else:
                return 0
            
        # return number of occurances once pattern is empty
        else:
            for i in range(top, bottom+1):
                patternOccurences.append(suffArray[i])
            return patternOccurences
    


# In[59]:


# returns the number of times each pattern occurs in text using its BWT 
#def MatchAllPatterns(text, patterns):
def MatchAllPatterns(genome, reads):
    lastColumn = BWT(genome+'$')
    count = CountDict(lastColumn)
    firstOccur = FirstOccurenceDict(lastColumn)
    suffArray = SuffixArray(genome+'$')
    
    matchNums = []
    for header, pattern in reads.items():
        matchIndices = BetterBWTMatching(firstOccur, lastColumn, pattern, count, suffArray)
        if matchIndices != 0:
            matchIndices = sorted(matchIndices)
        matchNums.append((header, pattern, matchIndices))
    
    return matchNums


# In[60]:

def main():
    sam = MatchAllPatterns(genome, readDict)
    with open('myBwaMatching.txt', 'w') as file:
        countLines = 1
        for a in sam:
            file.write(a[0])
            file.write('\t')
            file.write(a[1])
            file.write('\t')
            countIndices = 1
            if a[2] != 0:
                file.write(' ')
                for index in a[2]:
                    file.write(str(index))
                    if countIndices < len(a[2]):
                        file.write(' ')
                    countIndices = countIndices + 1
            if countLines < len(sam):
                file.write('\n')
            countLines = countLines + 1

    with open('myBwaMatching.txt', 'r') as file:
        for line in file:
            line = line.strip()
            print(line)

main()
            
#if __name___ == "__main__":
   # main()



