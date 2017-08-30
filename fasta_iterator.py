__author__="ycherukuri"
__date__ ="8/4/16"

#==============================================================#
def fasta_iterator(file_handle):
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
        # do something to each record            
        print "next record" 
        # stop iteration                                                   
        if not line:return                                                      
        ####
    ####
#==============================================================#
def real_main():
    for record in fasta_iterator(open("test.fasta",'r')):
            print record
    ####
   
#==============================================================
if ( __name__ == '__main__' ):
    real_main()



