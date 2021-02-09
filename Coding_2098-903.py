import re #built in packages
import itertools
from copy import copy # could use [:] in most cases instead, but this is better to maintain and more readable

# bad practice to include but done to keep everything in one file
TranslationTable = [( "UUU", "Phe", 147.1766,  0),
    ("UUC", "Phe", 147.1766, 0),
    ("UUA", "Leu", 113.1594, 0),
    ("UUG", "Leu", 113.1594, 0),
    ("UCU", "Ser", 87.0782, 0),
    ("UCC", "Ser", 87.0782, 0),
    ("UCA", "Ser", 87.0782, 0),
    ("UCG", "Ser", 87.0782, 0),
    ("UAU", "Tyr", 163.1760, 0),
    ("UAC", "Tyr", 163.1760, 0),
    ("UGU", "Cys", 103.1388, 0),
    ("UGC", "Cys", 103.1388, 0),
    ("UGG", "Trp", 186.2132, 0),
    ("CUU", "Leu", 113.1594, 0),
    ("CUC", "Leu", 113.1594, 0),
    ("CUA", "Leu", 113.1594, 0),
    ("CUG", "Leu", 113.1594, 0),
    ("CCU", "Pro", 97.1167, 0),
    ("CCC", "Pro", 97.1167, 0),
    ("CCA", "Pro", 97.1167, 0),
    ("CCG", "Pro", 97.1167, 0),
    ("CAU", "His", 137.1411, 1),
    ("CAC", "His", 137.1411, 1),
    ("CAA", "Gln", 128.1307, 0),
    ("CAG", "Gln", 128.1307, 0),
    ("CGU", "Arg", 156.1875, 1),
    ("CGC", "Arg", 156.1875, 1),
    ("CGA", "Arg", 156.1875, 1),
    ("CGG", "Arg", 156.1875, 1),
    ("AUU", "Ile", 113.1594, 0),
    ("AUC", "Ile", 113.1594, 0),
    ("AUA", "Ile", 113.1594, 0),
    ("AUG", "Met", 131.1926, 0),
    ("ACU", "Thr", 101.1051, 0),
    ("ACC", "Thr", 101.1051, 0),
    ("ACA", "Thr", 101.1051, 0),
    ("ACG", "Thr", 101.1051, 0),
    ("AAU", "Asn", 114.1038, 0),
    ("AAC", "Asn", 114.1038, 0),
    ("AAA", "Lys", 128.1741, 1),
    ("AAG", "Lys", 128.1741, 1),
    ("AGU", "Ser", 87.0782, 0),
    ("AGC", "Ser", 87.0782, 0),
    ("AGA", "Arg", 156.1875, 1),
    ("AGG", "Arg", 156.1875, 1),
    ("GUU", "Val", 99.1326, 0),
    ("GUC", "Val", 99.1326, 0),
    ("GUA", "Val", 99.1326, 0),
    ("GUG", "Val", 99.1326, 0),
    ("GCU", "Ala", 71.0788, 0),
    ("GCC", "Ala", 71.0788, 0),
    ("GCA", "Ala", 71.0788, 0),
    ("GCG", "Ala", 71.0788, 0),
    ("GAU", "Asp", 115.0886, -1),
    ("GAC", "Asp", 115.0886, -1),
    ("GAA", "Glu", 129.1155, -1),
    ("GAG", "Glu", 129.1155, -1),
    ("GGU", "Gly", 57.0519, 0),
    ("GGC", "Gly", 57.0519, 0),
    ("GGA", "Gly", 57.0519, 0),
    ("GGG", "Gly", 57.0519, 0)];

inputString = input("Please input a DNA sequence string or file path\n")
substringList = [""]
killMode = True

#read from file if a file was input, otherwise read from input
if("\\"in inputString):
    DNA_File = open(inputString, "r")
    sequence = DNA_File.read()
    DNA_File.close()
else:
    sequence = inputString

# just a function for concatenating all strings in a list, will be used later
def CombineStrings(list1):
    finalString = ""
    for string in list1:
        finalString = finalString + string
    return(finalString)

# seperate the sequence into substrings that start at promoter sequences and end at either other promoter squences or terminator sequences
substringList[0] = sequence
for string in re.findall("[CT][CT]A[ATGC][AT][CT][CT]|TATAAA|CGCGCGCGAAACGCGCGCGTTTTTTT", sequence):
    tempEndStrings = substringList[-1].split(string, 1)
    substringList[-1] = tempEndStrings[0]
    substringList.append(tempEndStrings[1])
    if(killMode):
        substringList.pop(-2)
    if(string == "TATAAA"):
        substringList[-1] = "A" + substringList[-1]
        killMode = False
    elif(string == "CGCGCGCGAAACGCGCGCGTTTTTTT"):
        substringList[-2] = substringList[-2] + "CGCGCGCGAAACGCGCGCGTTTTTTT"
        killMode = True #Killmode lets us treat the hairpin like a promoter in this function but still gets rid of all of the data after it and before the next promoter
    else:
        substringList[-1] = string[2:] + substringList[-1]
        killMode = False

# if the sequence doesnt end with a terminator, destroy the last item in the array since it isnt transcribed
x = re.findall("world$", substringList[-1])
if(not x):
    substringList.pop(-1)

#Convert all thymines to uracils
for i in range(len(substringList)):
    substringList[i] = re.sub("T", "U", substringList[i])

# find exons in each substring
exons = [""] * len(substringList)
for j in range(len(substringList)):
    intronsClipped = re.sub("GU[AG]AGU.*?CAG", "-", substringList[j]) # find all introns and replace with -
    exons[j] = intronsClipped.split("-") # split string into exons at each intron

# Find all possible exon combos
possibleExons = [[]] * len(exons)
for i in range(len(exons)):
    if(len(exons[i]) <= 1):
        possibleExons[i] = copy(exons[i])
    else:
        possibilities = list(itertools.product([0, 1], repeat=len(exons[i])-2)) # find all possible true/false combinations for the middle terms (cartesian multiplication)
        p = []
        for j in range(len(possibilities)):
            p.append(list(possibilities[j]))
            p[j].insert(0,1) # first and last term should always be true
            p[j].append(1)
            possibilities[j] = p[j]
            possibleExons[i].append(copy(exons[i])) # initialize possibleExons
        for x in range(len(possibilities)):
            for y in range(len(possibilities[x])):
                if(possibilities[x][y] == 0):
                    possibleExons[i][x][y] = ""
        #concatenate strings
        possibleExons[i] = [CombineStrings(exonPossibility) for exonPossibility in possibleExons[i]]

# find proteins for each exon combination
# split strings into ones starting with the start codon (AUG)
possibleProteins = []
for substring in possibleExons:
    for exonCombo in substring:
        possibleProteins.append(copy(exonCombo))
        for startCodon in re.findall("(?<!^)AUG", exonCombo):
            startCodonIndex = re.search("(?<!^)AUG", possibleProteins[-1]).span()[0]
            possibleProteins.append(copy(possibleProteins[-1][startCodonIndex:]))

#filter out anything added temporarily (wont start with AUG)
tempPossibleProteins = []
for possibility in possibleProteins:
    if(possibility[0:3] == "AUG"):
        tempPossibleProteins.append(possibility)
possibleProteins = copy(tempPossibleProteins)

# Take out everything starting at the stop codon (UGA, UAA, UAG)
tempPossibleProteins = []
for possibility in possibleProteins:
    for i in range(len(possibility)):
        codon = possibility[3*i:3*i+3]
        if(codon == "UGA" or codon == "UAA" or codon == "UAG"):
            tempPossibleProteins.append(possibility[:3*i])
possibleProteins = copy(tempPossibleProteins)

# remove duplicates (bad for perf with larger sets)
possibleProteins = list(set(possibleProteins))

# Convert proteins into amino acids
finalProteins = []
for possibility in possibleProteins:
    finalProteins.append(["", 0, 0])
    for i in range(len(possibility)):
        codon = possibility[3*i:3*i+3]
        for translationElement in TranslationTable:
            if(codon == translationElement[0]):
                finalProteins[-1][0] = finalProteins[-1][0] + translationElement[1] #abbreviation
                finalProteins[-1][1] = finalProteins[-1][1] + translationElement[2] #mass
                finalProteins[-1][2] = finalProteins[-1][2] + translationElement[3] #charge

# Display results
for protein in finalProteins:
    print("{0}, {1}u, {2}e".format(protein[0], protein[1], protein[2]))
input()