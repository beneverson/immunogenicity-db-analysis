import sqlite3
import re
import csv
from collections import defaultdict
import math
import operator

def regexp(expr, item):
  reg = re.compile(expr)
  return reg.search(item) is not None 

def calculateProportions(inDict, outDict):
  # go through the whole dictionary 
  for key in inDict.keys():
    # add the proportion of each character to outDict
    outDict[key] = float(1.0*inDict[key]/sum(inDict.values()))

# return the standard error of the difference between proportion 1 and 2
def calculateSE(p1, n1, p2, n2):
  f1 = (p1*(1-p1))/n1
  f2 = (p2*(1-p2))/n2
  return math.sqrt(f1+f2)

# return the Z-value for the given arguments
def calculateZ(p1, n1, p2, n2):
  return (p1-p2)/calculateSE(p1,n1,p2,n2)

# calculate the z-values for the given dictionaries
def calculateZvals(epiCountDict, epiPropDict, srcCountDict, srcPropDict, zDict):
  for key in epiCountDict: # go through the dict of counts
    if key in srcPropDict.keys():
      n2 = sum(epiCountDict.values()) 
      p2 = epiPropDict[key]
      n1 = sum(srcCountDict.values())
      p1 = srcPropDict[key]
      zDict[key] = calculateZ(p1, n1, p2, n2)

# populate aa_list by filtering ss_string for ssElement, and adding the
# aa element at the corresponding location 
def collectAAsFromSSElements(ssElement, ss_string, aa_string, aa_list):
  regex = ssElement
  for match in re.finditer(regex, ss_string):
    aa_list[aa_string[match.start():match.end()]] += 1
  
# connection to db for reading data
readConn = sqlite3.connect('EPITOPES.sqlite')
readConn.create_function("REGEXP", 2, regexp)
readCursor = readConn.cursor()

epitope_AA_counts = defaultdict(int) # counts of AA's in all epitopes
epitope_SS_counts = defaultdict(int) # counts of SS's in all epitopes
epitope_AA_proportions = {} # proportions
epitope_SS_proportions = {} # proportions

source_AA_counts = defaultdict(int) # counts of AA's in all source sequences
source_SS_counts = defaultdict(int) # counts of SS's in all source sequences
source_AA_proportions = {} # proportions
source_SS_proportions = {} # proportions

AA_zvalues = {} # z values of all AA's between epitope and source
SS_zvalues = {} # z values of all SS's between epitope and source

readCursor.execute("SELECT * FROM combined")

# -------------------------------------------
# -------------------------------------------
# SET THE REGEX HERE ------------------------
regex_ss_character = ur"[T,S,C,*]"

# go through each row in the database
for row in readCursor:
  # get the epitope sequence
  epitope_seq = row[1]
  epitope_str = row[2]
  
  # collect counts of helical AA's into epitope_AA_counts
  collectAAsFromSSElements(regex_ss_character, epitope_str, epitope_seq, epitope_AA_counts) 
    
  # get the source sequence
  source_seq = row[5]
  source_str = row[6]
  
  # collect counts of helical AA's into epitope_AA_counts
  collectAAsFromSSElements(regex_ss_character, source_str, source_seq, source_AA_counts)

# compute the 'proportions' dictionaries
calculateProportions(epitope_AA_counts, epitope_AA_proportions)
# calculateProportions(epitope_SS_counts, epitope_SS_proportions)
calculateProportions(source_AA_counts, source_AA_proportions)
# calculateProportions(source_SS_counts, source_SS_proportions)

# compute the zValue dictionaries
calculateZvals(epitope_AA_counts, epitope_AA_proportions, source_AA_counts, source_AA_proportions, AA_zvalues)
# calculateZvals(epitope_SS_counts, epitope_SS_proportions, source_SS_counts, source_SS_proportions, SS_zvalues)

# print out the statistically significant differences 
print("-----------------------------------STRUCTURAL-AMINO-ACIDS---------------------------------------------------")
print "SS ELEMENT IS: " + regex_ss_character;

#for aaEpiKey in epitope_AA_proportions:
  #if epitope_AA_proportions[aaEpiKey] > 0.001:
    # print aaEpiKey + " has a proportion of " + str(epitope_AA_proportions[aaEpiKey]) + " in epitopes."
    
#for aaSrcKey in source_AA_proportions:
  #if source_AA_proportions[aaSrcKey] > 0.001:
    # print aaSrcKey + " has a proportion of " + str(source_AA_proportions[aaSrcKey]) + " in source sequences."


print("-----------------------------------SIGNIFICANCE-------------------------------------------------------------")
for aaKey in AA_zvalues:
  if AA_zvalues[aaKey] > 1.96 or AA_zvalues[aaKey] < -1.96:
    print aaKey + " has a z-value of " + str(AA_zvalues[aaKey])
