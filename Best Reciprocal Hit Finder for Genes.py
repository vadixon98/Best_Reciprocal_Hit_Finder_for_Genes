Python 3.13.0 (tags/v3.13.0:60403a5, Oct  7 2024, 09:38:07) [MSC v.1941 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
import sys
sys.setrecursionlimit(100000)
from humanChickenProteins import *
from blosum62 import *

def alignScore(S1, S2, gap, substitutionMatrix):
    ''' Returns the sequence alignment score for S1 and S2.'''
    if S1 == " ": return gap * len(S2)
    elif S2 == " " : return gap * len(S1)
    else:
        option1 = substitutionMatrix[(S1[0], S2[0])] + alignScore(S1[1:], S2[1:], gap, substitutionMatrix)
        option2 = gap + alignScore(S1[1:], S2, gap , substitutionMatrix)
        option3 = gap + alignScore(S1, S2[1:], gap , substitutionMatrix)
        return max(option1, option2, option3)

    

def memoAlignScore(S1, S2, gap, substitutionMatrix, memo):
    '''computes the alignment score using a a memoized version of alignScore,
memoAlignScore(S1, S2, gap, substitutionMatrix,memo), where memo is a
dictionary.'''
    # Check if result is already in memo
    if (S1, S2) in memo:
        return memo[(S1, S2)]
    
    # Base cases
    if S1 == "": return gap * len(S2)
    elif S2 == "": return gap * len(S1)
    
    # Recursive cases
    option1 = substitutionMatrix[(S1[0], S2[0])] + memoAlignScore(S1[1:], S2[1:], gap, substitutionMatrix, memo)
    option2 = gap + memoAlignScore(S1[1:], S2, gap, substitutionMatrix, memo)
    option3 = gap + memoAlignScore(S1, S2[1:], gap, substitutionMatrix, memo)
    
    # Compute max score
    best_score = max(option1, option2, option3)    
    # Store result in memo
    memo[(S1, S2)] = best_score
    
    return best_score


def allScores(geneList1, geneList2):
    '''to obtain the alignment score between all proteins in two species.
allScores takes as input two lists of genes. It then takes every protein in
the first list and calls memoAlignScore on it with every protein in the second
list. '''
    # Dictionary to store alignment scores with (gene1, gene2) as keys
    scores = {}
    memo = {}  # Memoization dictionary for memoAlignScore
    
    for gene1 in geneList1:
        for gene2 in geneList2:
            # Get protein sequences for the genes from geneD
            protein1 = geneD[gene1][3]
            protein2 = geneD[gene2][3]
            
            # Compute alignment score using memoAlignScore
            score = memoAlignScore(protein1, protein2, -9, blosum62 , memo)
            
            # Store the score in scores dictionary
            scores[(gene1, gene2)] = score
    
    return scores


def closestMatch(geneName, allScoresD):
    '''Given a gene name and a dictionary of alignment scores, closestMatch
returns the protein from the other species which is most similar
(has the highest alignment score).'''
    # Initialize variables to keep track of the highest score and best matching protein
    max_score = -float('inf')  # Use negative infinity to ensure any score will be higher
    best_match = None

    # Loop through each key-value pair in the dictionary
    for (geneA, geneB), score in allScoresD.items():
        # Check if geneName is in the current tuple (geneA, geneB)
        if geneName in (geneA, geneB):
            # Determine the other gene in the pair
            other_gene = geneB if geneA == geneName else geneA
            
            # Update max_score and best_match if this score is higher
            if score > max_score:
                max_score = score
                best_match = other_gene

    return best_match


def printBRH(geneName, allScoresD):
    ''' Prints the Best Reciprocal  Hit (BRH) for a given geneName.
If there is no best reciprocal hit it returns without printing anything. '''
    # store the best match for geneName in the other species
    best_match = closestMatch(geneName, allScoresD)
    
    # Checka if calling closestMatch on best_match returns geneName
    if best_match and closestMatch(best_match, allScoresD) == geneName:

        # Get gene information from geneD for both genes
        chromosome1, start1, a1, b1 = geneD[geneName]
        chromosome2, start2, a2, b2 = geneD[best_match]
        
        # Prints gene information for both species
        print(f"{chromosome1} {start1} {geneName} --- {chromosome2} {start2} {best_match}")



def runBRHSample():
    '''Print best reciprocal hits for sample data. First in human
    chromosome order, then in chicken chromosome order.'''
    allScoresD=allScores(sampleHumanGeneList,sampleChickenGeneList)
    print ('human --- chicken')
    for geneName in sampleHumanGeneList:
        printBRH(geneName,allScoresD)
    print
    print ('chicken --- human')
    for geneName in sampleChickenGeneList:
        printBRH(geneName,allScoresD)

        

def runBRH():
    '''Print best reciprocal hits for full data. First in human
    chromosome order, then in chicken chromosome order.'''
    allScoresD=allScores(humanGeneList,chickenGeneList)
    print ('human --- chicken')
    for geneName in humanGeneList:
        printBRH(geneName,allScoresD)
    print
    print ('chicken --- human')
    for geneName in chickenGeneList:
        printBRH(geneName,allScoresD)

        
