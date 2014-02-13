import sqlite3
import re
import csv
from collections import defaultdict
import math
import operator

def regexp(expr, item):
  reg = re.compile(expr)
  return reg.search(item) is not None 
  
# populate patternCounts dictionary with substrings that match both ssRegex and aaRegex
def collectPatterns(ssRegex, aaRegex, ssString, aaString, patternCounts, minPatternLength):
  matchesList = [] # create a list of tuples of start and end positions of matching substrings
  for ssMatch in re.finditer(ssRegex, ssString): # match on secondary structure
    matchesList.append((ssMatch.start(),ssMatch.end()))
  for aaMatch in re.finditer(aaRegex, aaString): # match on AA identity
    # now look if the aaMatch bounds and ssMatch bounds overlap in any way
    for matchBoundsTuple in matchesList:
      ssBegins = matchBoundsTuple[0];
      ssEnds = matchBoundsTuple[1];
      # check if the substrings are overlapping in any way 
      if( not( (aaMatch.end() < ssBegins) or (ssEnds < aaMatch.start()) ) ):
        # if so, find the largest beginning and the smallest ending AKA the most conservative overlap
        beginningIndex = max(ssBegins, aaMatch.start())
        endingIndex = min(ssEnds, aaMatch.end())
        # check to be sure the substring is large enough
        if ( (endingIndex - beginningIndex) >= minPatternLength ):
          # and add the appropriate substrings to the pattern dict
          patternTuple = (ssString[beginningIndex:endingIndex], aaString[beginningIndex:endingIndex])
          patternCounts[patternTuple] += 1
          
# compile a count of the amino acids  founds in a patternCounts dictionary
def countAAs(patternCounts, aaCounts):
  for patternTuple in patternCounts.keys(): # go thru all pattern tuples
    aaString = patternTuple[1]
    for index in range(len(aaString)):
     aaCounts[aaString[index]] += 1

def calculateProportions(inDict, outDict):
  # go through the whole dictionary 
  for key in inDict.keys():
    # add the proportion of each character to outDict
    outDict[key] = float(1.0*inDict[key]/sum(inDict.values()))

# connection to db for reading data
readConn = sqlite3.connect('EPITOPES.sqlite')
readConn.create_function("REGEXP", 2, regexp)
readCursor = readConn.cursor()

epitope_pattern_counts = defaultdict(int) # counts of arbitrary patterns in the epitope strings
epitope_AA_counts = defaultdict(int) # counts of AA's in all epitopes
epitope_AA_proportions = defaultdict(int) # counts of SS's in all epitopes

source_pattern_counts = defaultdict(int) # counts of arbitrary pattern in the source strings
source_AA_counts = defaultdict(int) # counts of AA's in all source sequences
source_AA_proportions = defaultdict(int) # counts of SS's in all source sequences

readCursor.execute("SELECT * FROM COMBINED_NONREDUNDANT") # select from all known binders

# -------------------------------------------
# -------------------------------------------
# SET THE REGEX HERE ------------------------
regex_ss_character = ur"[G,H,I]+" # beta sheet secondary structure
regex_aa_character = ur"([F,L,I,M,V][R,K,D,E,N,Q])+" # pattern of polar/nonpolar amino acids

# go through each row in the database
for row in readCursor:
  # get the epitope sequence
  epitope_seq = row[1]
  epitope_str = row[2]

  collectPatterns(regex_ss_character, regex_aa_character, epitope_str, epitope_seq, epitope_pattern_counts, 3)
  
  # get the source sequence
  source_seq = row[4]
  source_str = row[5]
  
  collectPatterns(regex_ss_character, regex_aa_character, source_str, source_seq, source_pattern_counts, 3)

countAAs(epitope_pattern_counts,epitope_AA_counts)
calculateProportions(epitope_AA_counts, epitope_AA_proportions)

print "Epitope Amino Acid Proportions:"
# print epitope_pattern_counts
print epitope_AA_proportions

countAAs(source_pattern_counts,source_AA_counts)
calculateProportions(source_AA_counts, source_AA_proportions)

# print "Source Amino Acid Proportions:"
# print source_pattern_counts
# print source_AA_proportions