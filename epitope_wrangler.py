import sqlite3
import re

def regexp(expr, item):
  if item is None:
    return False
  reg = re.compile(expr)
  return reg.search(item) is not None
  
def appendToVals(row, vals):
  rowLength = len(row)
  for x in range(0, rowLength):
    vals.append(row[x])
  
# connection for writing to the database
writeConn = sqlite3.connect('EPITOPES.sqlite')
writeConn.create_function("REGEXP", 2, regexp)
writeCursor = writeConn.cursor()

# connection to db for checking for duplicates
dupConn = sqlite3.connect('EPITOPES.sqlite')
dupConn.create_function("REGEXP", 2, regexp)
dupCursor = dupConn.cursor()

# connection to db for reading data
readConn = sqlite3.connect('EPITOPES.sqlite')
readConn.create_function("REGEXP", 2, regexp)
readCursor = readConn.cursor()
readCursor.execute("SELECT * FROM COMBINED_REDUNDANT")

# instantiate a list to hold all the nonredundant data we find
nonredundant_rows = []

for readRow in readCursor:
  epitope_seq = readRow[1] # readRow[1] is "epitope_seq"
  pdb_chain = readRow[3] # readRow[3] is "pdb_chain"
  
  # print "epitope_seq is " + epitope_seq
  # print "pdb_chain is " + pdb_chain
  
  # query the database for all epitopes that contain the same epitope sequence as the current row, as well as the same pdb chain
  queryVals = []
  queryVals.append(epitope_seq)
  queryVals.append(epitope_seq)
  queryVals.append(pdb_chain)
  dupCursor.execute("SELECT * FROM COMBINED_REDUNDANT WHERE \"epitope_seq\" REGEXP ? AND \"epitope_seq\" IS NOT ? AND \"pdb_chain\" IS ?", queryVals)
  
  # find the longest epitope sequence out of those returned
  longest_seq = epitope_seq
  longest_pdb = pdb_chain
  writeVals = []
  appendToVals(readRow, writeVals)
  for dupRow in dupCursor:
    current_seq = dupRow[1]
    current_pdb = dupRow[3]
    # print "current_seq is " + current_seq
    if current_seq.__len__() > longest_seq.__len__():
      # reset the longest_seq string
      # print "Found a longer sequence"
      longest_seq = current_seq
      longest_pdb = current_pdb
      # reset the writeVals list
      del writeVals[:]
      appendToVals(dupRow, writeVals)
  
  # writeVals now contains the longest distinct row, to be written to the non-redundant database
  # we'll hold all the nonredundant rows in a list before writing them to sql
  
  # first handle the case where nonredundant_rows is empty
  if not nonredundant_rows:
    nonredundant_rows.append(writeVals)

  # go through all of nonredundant_rows
  wasFound = False
  for nrVals in nonredundant_rows:
    # if writeVals is the same as nrVals
    if set(nrVals) == set(writeVals):
      # set our flag to 'True'
      wasFound = True
      
  # if writeVals was not found in nonredundant_rows
  if not wasFound:
    # append writeVals to nonredundant_rows
    nonredundant_rows.append(writeVals)
    
# when we're done, insert the contents of nonredundant_rows into COMBINED_NONREDUNDANT table
for newWriteValue in nonredundant_rows:
  writeCursor.execute("INSERT INTO COMBINED_NONREDUNDANT VALUES(?,?,?,?,?,?)", newWriteValue)

writeConn.commit()


  
      
    

		
	

  


  
    
  
  