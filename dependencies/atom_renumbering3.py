# -*- coding: utf-8 -*-
# This function renumbers the ATOM and HETATM records following the output from reduce such that they can be treated by FIRST

import sys

if (len(sys.argv) > 1):
    print (sys.argv[1])

    with open(sys.argv[1], 'r') as fin:
        lines = fin.readlines()
        
    with open(sys.argv[1], 'w') as fout:
        number=0
        for line in lines:
            if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
                if line[16] == 'A' or line[16] == ' ':
                    number = number + 1
                    whitespace = '     '
                    atom_number = whitespace[0:5-len(str(number))] + str(number)
                    fout.write(line.replace(line[6:11], atom_number, 1)) #only want to replace first occurance or can change atom coords!!
            else:
                fout.write(line)




                


"""

subprocess.call(pg.settings.REDUCE + ' -Quiet -BUILD -DB ' +
                pg.settings.FIRST_ROOT + '/lib/het_dict.txt ' 
                + '/home/maxhodges/proteinnetwork/proteingraph3/testing/shorttest.pdb' + ' > ' + '/home/maxhodges/proteinnetwork/proteingraph3/testing/shorttest.pdb'[0:-4] 
                + '_H_temp.pdb', shell=True)

"""
