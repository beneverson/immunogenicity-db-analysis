import sqlite3
import re

def regexp(expr, item):
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

for readRow in readCursor:
  epitope_seq = readRow[1] # readRow[1] is "epitope_seq" in MHCPEP
  source_seq = readRow[5] # readRow[5] is "source_seq" in MHCPEP
  
  # query the database for all epitopes that contain the same epitope sequence as the current row, as well as the same source sequence
  queryVals = []
  queryVals.append(epitope_seq)
  queryVals.append(epitope_seq)
  queryVals.append(source_seq)
  dupCursor.execute("SELECT * FROM combined WHERE \"epitope_seq\" REGEXP ? AND \"epitope_seq\" IS NOT ? AND \"source_seq\" IS ?", queryVals)
  
  # find the longest epitope sequence out of those returned
  longest_seq = epitope_seq
  writeVals = []
  appendToVals(readRow, writeVals)
  for dupRow in dupCursor:
    current_seq = dupRow[1]
    if current_seq.__len__() > longest_seq.__len__():
      # reset the longest_seq string
      longest_seq = current_seq
      # reset the writeVals list
      del writeVals[:]
      appendToVals(dupRow, writeVals)
  
  print 'longest_seq among duplicates is %s' % longest_seq
  
  # write the entry corresponding to the longest row of all duplicates into a new database
  writeCursor.execute("INSERT INTO combined2 values(?,?,?,?,?,?,?)", writeVals)

writeConn.commit()
  
      
    

		
	

  


  
    
  
  