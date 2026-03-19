"""
Best Reciprocal Hit Finder for Genes

This script finds Best Reciprocal Hits (BRH) between genes from two species
(human and chicken) by performing sequence alignment and identifying genes
that are each other's closest matches across species.

A Best Reciprocal Hit occurs when:
- Gene A from species 1 matches most closely with Gene B from species 2
- AND Gene B from species 2 matches most closely with Gene A from species 1

This is commonly used in comparative genomics to identify orthologous genes.
"""

import sys
# Increase recursion limit for deep recursive sequence alignment calculations
sys.setrecursionlimit(100000)

# Import gene data (humanChickenProteins module contains gene lists and dictionaries)
from humanChickenProteins import *
# Import BLOSUM62 substitution matrix for protein sequence alignment scoring
from blosum62 import *

def alignScore(S1, S2, gap, substitutionMatrix):
    """
    Returns the sequence alignment score for S1 and S2.
    
    Uses dynamic programming via recursion to find the optimal global alignment score
    between two protein sequences. Considers three operations at each position:
    1. Match/mismatch (align two characters)
    2. Insert gap in S1 (skip character in S1)
    3. Insert gap in S2 (skip character in S2)
    
    Args:
        S1 (str): First protein sequence
        S2 (str): Second protein sequence
        gap (int): Penalty score for inserting a gap
        substitutionMatrix (dict): Scoring matrix (e.g., BLOSUM62) for amino acid pairs
        
    Returns:
        int: Maximum alignment score between the two sequences
    """
    # Base case: if S1 is empty, score is gap penalty times length of S2
    if S1 == "": return gap * len(S2)
    # Base case: if S2 is empty, score is gap penalty times length of S1
    elif S2 == "": return gap * len(S1)
    else:
        # Option 1: Align first characters of both sequences (match/mismatch)
        option1 = substitutionMatrix[(S1[0], S2[0])] + alignScore(S1[1:], S2[1:], gap, substitutionMatrix)
        # Option 2: Insert gap in S1, skip first character of S2
        option2 = gap + alignScore(S1[1:], S2, gap , substitutionMatrix)
        # Option 3: Insert gap in S2, skip first character of S1
        option3 = gap + alignScore(S1, S2[1:], gap , substitutionMatrix)
        # Return the maximum score among the three options
        return max(option1, option2, option3)

    

def memoAlignScore(S1, S2, gap, substitutionMatrix, memo):
    """
    Computes the alignment score using a memoized version of alignScore.
    
    This is an optimized version that caches previously computed results to avoid
    redundant calculations, significantly improving performance for large sequences.
    
    Args:
        S1 (str): First protein sequence
        S2 (str): Second protein sequence
        gap (int): Penalty score for inserting a gap
        substitutionMatrix (dict): Scoring matrix (e.g., BLOSUM62) for amino acid pairs
        memo (dict): Dictionary to cache computed alignment scores (key: (S1, S2), value: score)
        
    Returns:
        int: Maximum alignment score between the two sequences
    """
    # Check if result is already in memo - return cached value if found
    if (S1, S2) in memo:
        return memo[(S1, S2)]
    
    # Base cases: handle empty sequences
    # If S1 is empty, score is gap penalty times length of S2
    if S1 == "": return gap * len(S2)
    # If S2 is empty, score is gap penalty times length of S1
    elif S2 == "": return gap * len(S1)
    
    # Recursive cases: consider three alignment options
    # Option 1: Align first characters of both sequences (match/mismatch)
    option1 = substitutionMatrix[(S1[0], S2[0])] + memoAlignScore(S1[1:], S2[1:], gap, substitutionMatrix, memo)
    # Option 2: Insert gap in S1, skip first character of S2
    option2 = gap + memoAlignScore(S1[1:], S2, gap, substitutionMatrix, memo)
    # Option 3: Insert gap in S2, skip first character of S1
    option3 = gap + memoAlignScore(S1, S2[1:], gap, substitutionMatrix, memo)
    
    # Compute max score among the three options
    best_score = max(option1, option2, option3)
    # Store result in memo for future use (avoids recomputation)
    memo[(S1, S2)] = best_score
    
    return best_score


def allScores(geneList1, geneList2):
    """
    Obtain the alignment score between all proteins in two species.
    
    This function computes pairwise alignment scores for every combination of genes
    from the two input lists. It creates a comprehensive score matrix that can be
    used to find the best matches between species.
    
    Args:
        geneList1 (list): List of gene names from first species
        geneList2 (list): List of gene names from second species
        
    Returns:
        dict: Dictionary mapping (gene1, gene2) tuples to their alignment scores
              Key format: (gene from list1, gene from list2)
              Value: alignment score (int)
    """
    # Dictionary to store alignment scores with (gene1, gene2) as keys
    scores = {}
    # Memoization dictionary shared across all alignments for efficiency
    # This allows reuse of subproblem solutions when aligning different gene pairs
    memo = {}
    
    # Iterate through all pairs of genes from the two species
    for gene1 in geneList1:
        for gene2 in geneList2:
            # Get protein sequences for the genes from geneD dictionary
            # geneD[gene][3] contains the protein sequence (index 3)
            protein1 = geneD[gene1][3]
            protein2 = geneD[gene2][3]
            
            # Compute alignment score using memoized alignment function
            # Gap penalty: -9, Substitution matrix: BLOSUM62
            score = memoAlignScore(protein1, protein2, -9, blosum62 , memo)
            
            # Store the score in scores dictionary for later lookup
            scores[(gene1, gene2)] = score
    
    return scores


def closestMatch(geneName, allScoresD):
    """
    Find the protein from the other species which is most similar to the given gene.
    
    Given a gene name and a dictionary of alignment scores, this function identifies
    the gene from the other species that has the highest alignment score with the
    input gene. This represents the best match or closest homolog.
    
    Args:
        geneName (str): Name of the gene to find the closest match for
        allScoresD (dict): Dictionary of alignment scores with (gene1, gene2) tuples as keys
        
    Returns:
        str: Name of the gene from the other species with the highest alignment score,
             or None if no matches found
    """
    # Initialize variables to keep track of the highest score and best matching protein
    # Use negative infinity to ensure any valid score will be higher
    max_score = -float('inf')
    best_match = None

    # Loop through each key-value pair in the dictionary
    for (geneA, geneB), score in allScoresD.items():
        # Check if geneName is in the current tuple (geneA, geneB)
        # This means we're looking at a pair that involves our gene of interest
        if geneName in (geneA, geneB):
            # Determine the other gene in the pair (the one from the other species)
            other_gene = geneB if geneA == geneName else geneA
            
            # Update max_score and best_match if this score is higher than previous best
            if score > max_score:
                max_score = score
                best_match = other_gene

    return best_match


def printBRH(geneName, allScoresD):
    """
    Print the Best Reciprocal Hit (BRH) for a given geneName.
    
    A Best Reciprocal Hit occurs when:
    - Gene A's closest match is Gene B
    - AND Gene B's closest match is Gene A
    
    This indicates a strong bidirectional relationship suggesting orthology.
    If there is no best reciprocal hit, the function returns without printing anything.
    
    Args:
        geneName (str): Name of the gene to find and print BRH for
        allScoresD (dict): Dictionary of alignment scores with (gene1, gene2) tuples as keys
    """
    # Store the best match for geneName in the other species
    best_match = closestMatch(geneName, allScoresD)
    
    # Check if calling closestMatch on best_match returns geneName (reciprocal relationship)
    # This verifies that the match is bidirectional (BRH criteria)
    if best_match and closestMatch(best_match, allScoresD) == geneName:
        # Get gene information from geneD for both genes
        # geneD format: [chromosome, start_position, ?, protein_sequence]
        chromosome1, start1, a1, b1 = geneD[geneName]
        chromosome2, start2, a2, b2 = geneD[best_match]
        
        # Print gene information for both species in format:
        # chromosome start_position gene_name --- chromosome start_position gene_name
        print(f"{chromosome1} {start1} {geneName} --- {chromosome2} {start2} {best_match}")



def runBRHSample():
    """
    Print best reciprocal hits for sample data.
    
    This function processes a smaller sample dataset for testing purposes.
    Results are printed first in human chromosome order, then in chicken chromosome order.
    Uses sample gene lists (sampleHumanGeneList, sampleChickenGeneList) instead of full datasets.
    """
    # Compute all pairwise alignment scores between sample human and chicken genes
    allScoresD = allScores(sampleHumanGeneList, sampleChickenGeneList)
    
    # Print header and process human genes (in chromosome order)
    print('human --- chicken')
    for geneName in sampleHumanGeneList:
        printBRH(geneName, allScoresD)
    
    # Print separator line
    print()
    
    # Print header and process chicken genes (in chromosome order)
    print('chicken --- human')
    for geneName in sampleChickenGeneList:
        printBRH(geneName, allScoresD)

        

def runBRH():
    """
    Print best reciprocal hits for full data.
    
    This is the main function that processes the complete datasets.
    Results are printed first in human chromosome order, then in chicken chromosome order.
    Uses full gene lists (humanGeneList, chickenGeneList) - this may take considerable time
    to compute all pairwise alignments.
    """
    # Compute all pairwise alignment scores between all human and chicken genes
    # WARNING: This may take a long time for large datasets
    allScoresD = allScores(humanGeneList, chickenGeneList)
    
    # Print header and process all human genes (in chromosome order)
    print('human --- chicken')
    for geneName in humanGeneList:
        printBRH(geneName, allScoresD)
    
    # Print separator line
    print()
    
    # Print header and process all chicken genes (in chromosome order)
    print('chicken --- human')
    for geneName in chickenGeneList:
        printBRH(geneName, allScoresD)

        
