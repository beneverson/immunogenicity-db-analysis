The files in this directory are being used as part of a project to determine statistics on the space of known immunogenic regions in proteins with known structures. 

XLSX files: 

Initial reformats of various databases of experimentally-determined epitope sequences in proteins. 

.py files:

Various scripts for analyzing sqlite3 databases of large numbers of peptide sequences. In particular, db_regex_analyzer.py will determine statistics on epitopes whose string of secondary structure codes (DSSP) match a given regular expression input. 

EPITOPES.sqlite: 

The master database of experimentally-validated epitope sequences, secondary structure strings, PDB acession codes, and parent protein sequences. 
