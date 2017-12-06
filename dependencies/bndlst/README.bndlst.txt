README.bndlst                                                 JMW - 4/9/99

The program "bndlst" lists connectivity information for atoms in a PDB
file. It puts the information out in a tabular form, separated by
colons for easy parsing with programs such as awk. A record is output
for each atom, followed by a record for each atom it is connected to,
terminating in an EOL record signaling the end of the list.

By default, only the H-bonds and covalent nearest bonded neighbors are
listed. Alternatively, by using the -neighbors, -rad# and -cut# it is
possible to show non-bonded non-hbonds, and not show close bonded atoms.

Example

   bndlst 1xyzH > 1xyzH.bonds


To build an executable for bndlst:
1) cp makefile.xxx makefile  ## where xxx is the operating system type
2) --modify makefile, if neccessary-- 
3) make                      ## compiles and links the program

Mike Word
mike.word@duke.edu

Richardson Lab
Biochemistry Department
Duke University
Durham, NC USA 27710
