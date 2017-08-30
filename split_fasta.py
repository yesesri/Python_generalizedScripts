#!/usr/bin/python 
__author__="ycherukuri"
__date__ ="8/31/16"
#==============================================================#
#Install subprocess ; opt parser and make sure it is in the python path
#how to run the code
### python split_fasta.py  fastafile   no. of Files
#==============================================================#
import subprocess
import optparse
from hagsc_lib import  writeFASTA
#==============================================================#
def iterFASTA (file_handle):
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
#==============================================================#
def writeFASTA(batch,outfile):
    for record in batch : outfile.write(">%s\n%s\n%s\n"%(record['id'],record['seq'],record['desp']))
    outfile.close()
#==============================================================#
def real_main():
    readFile = "%s"%(argv[1])
    nFiles   = int( argv[2] )
    ext      = readFile[-2:]
    if(ext == "gz")    : p  = subprocess.Popen('cat |gunzip -c %s'%readFile,shell = True, stdout = subprocess.PIPE)
    elif(ext == "bz2") : p  = subprocess.Popen('cat |bunzip2 -c %s'%readFile,shell = True, stdout = subprocess.PIPE)
    else :               p  = subprocess.Popen('cat %s'%readFile,shell = True, stdout = subprocess.PIPE)
    totalbasecount = 0
    for record in iterFASTA(p.stdout): totalbasecount+=len(record.seq)
    print totalbasecount
    p.poll()
    split_cutoff    = totalbasecount/20
    recount_totalbase = 0  
    base_count = 0   
    file_count = 1  
    batch = []
    proc = subprocess.Popen('cat wheatgrass.Scaffolds.fasta.gz|gunzip -c',shell = True, stdout = subprocess.PIPE)
    for record in iterFASTA(proc.stdout):
        base_count+= len(record)
        batch.append(record)
        if(base_count >= split_cutoff):
           recount_totalbase+=base_count
           writeFASTA(batch,open("split_targetfile_%i.fasta"%(file_count), 'w'))
           print("file : split_targetfile_%i.fasta, basecount:%i " %(file_count,base_count))
           file_count+=1
           base_count = 0   
           batch = []
        ####
    ####
    proc.poll()
    recount_totalbase+=base_count
    writeFASTA(batch,gzip.open("split_targetfile_%i.fasta"%(file_count), 'w'))
    print("file : split_targetfile_%i.fasta, basecount:%i " %(file_count,base_count))
    #check the files with *fasta extension before running the shell script
    subprocess.call('sh blat_batchrun.sh',shell = True)
    oh = open('mergedfile.blat','w')
    for each_split_file in range(1,7):
        oh_tmp = open("split_targetfile_%i.blat"%each_split_file,'r')
        for n in xrange(5): oh_tmp.readline()
        for line in oh_tmp: oh.write( line )
        oh_tmp.close()
    ####
    oh.close()
#==============================================================#
if ( __name__ == '__main__' ):
    real_main()



 