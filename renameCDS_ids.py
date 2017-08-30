__author__="ycherukuri@hudsonalpha.org"
__date__ ="Created:  05/1/2017"
#======================================================================#
#script to rename duplicate read names in fasta file
#python renameCDS_ids.py fastafile
#======================================================================#
from sys import argv
import subprocess 
import re
from os.path import join
from sys import stderr
#======================================================================#
def iterFASTA(file_handle):
    while True:
        line = file_handle.readline()
        #skip blank lines
        if line == "": return
        ####
        if line[0] == ">": break
    ####
    while True:
        splitID = line.split(None)
        #split the header into id and description
        id      = splitID[0]
        id      = id[1:].strip()
        desp    = ''.join(splitID[1:])
        #create an empty list to store sequence lines of a read
        seq_list = []
        line = file_handle.readline()
        while True:
            # end of line
            if not line: break
            ####
            # break loop if header of next read is reached
            if line[0] == ">": break
            ####
            # append sequence lines to list
            seq_list.append(line.strip())
            line = file_handle.readline()
        ####
        #yield each read record
        yield {'id': id, 'desp': desp, 'seq': "".join(seq_list)}
        # stop iteration
        if not line:return
####
#======================================================================#
def writeFASTA(recordDict,outfile):
    outfile.write(">%s\n"recordDict['id'])
    outfile.write("%s\n"recordDict['seq'])
####
#=======================================================================#
def real_main():
    fastafile   = "%s"%(argv[1])
    outfilename = "%s"%(argv[2])
    system("cat %s | grep ">" | sed 's/'>'//g' | sort | uniq -c | awk '{if($1>1)print $0}' > DuplicateIDs.dat "%(fastafile))
	if(True):
	    ReadIdDict = {}
	    tmpDict = {}
	    outfile = open("%s.IDrenamed.fasta"%(outfilename),'w')
	    for line in open("DuplicateIDs.dat",'r'):
	         splitLine = line.split(None)
	         count     = int(splitLine[0])
	         ReadId    = splitLine[1].strip()
	         ReadIdDict[ReadId] = count
	    ####
	    for record in iterFASTA(open("%s"%fastafile,'r')):
	        if  (">%s"%record.id in ReadIdDict.keys()): 
	             try   : tmpDict[record['id']]+=1
	             except: 
	                 tmpDict[record['id']] = 1
                     newrecord = {'id': "%s_%i"%(record['id'], ), 'desp': '' , 'seq':record['seq']  }
                     writeFASTA(newrecord,outfile)
	         ####
	        else :writeFASTA(newrecord,outfile)
	    #### 
	    outfile.close()
####          
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
