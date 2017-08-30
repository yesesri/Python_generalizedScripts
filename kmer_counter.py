
#!/usr/bin/python
_author__="ycherukuri"
__date__ ="8/12/16"
#==============================================================#
#Install subprocess ; opt parser and make sure it is in the python path
#===================================================================#
from string import maketrans
import optparse
import subprocess
import re 
import math
import os.path
find_longest_string = re.compile(r'(1+)').finditer
trans_table = maketrans("ATGC","TACG")
#==================================================================#
def fastq_iterator(file_handle):
    while True:
        seq_id    = file_handle.stdout.readline()
        seq       = file_handle.stdout.readline()
        qual_id   = file_handle.stdout.readline()
        qual      = file_handle.stdout.readline()
        if not qual: return
        if qual_id.strip() == '+':
            yield {'ID': seq_id.strip(), 'seq': seq.strip(), 'qual': qual.strip()}
        else:
            print "inappropriate Input file",seq_id
            break
        ####
    ####
#===================================================================#
def get_phred_score_table(phred_offset):
    return dict((chr(char + phred_offset),char) for char in range(0, 44))
####
#===================================================================#
def get_score_table(window_size,cut_off):
    score_dict = {}
    for n in range((2*window_size),(43*window_size)):
        if n >= cut_off:
           score_dict[n] = '1'
        else:
           score_dict[n] = '0'
        ####
    return score_dict
####
#======================================================================#
def fastq_trimmer(record,window_size,phred_score_table,score_table):
    #Ascii to phred scale conversion
    try:
      converted_qual = [phred_score_table[each_char] for each_char in record['qual']]
    except KeyError:
        print "Invalid character encountered"
        print record['qual']
        assert False
    ####
    # appending binary scores based on average window score
    join_score     = "".join(score_table[sum(converted_qual[n:(n + window_size)])]
                             for n in range(len(record['seq']) - window_size + 1))
    #finding the quality fragment of a read
    long_string    = [(len(substring.group()), substring.span()[0]) for substring in find_longest_string(join_score)]
    long_string    =  sorted(long_string,reverse= True,key= lambda x:x[0])
    try:
       right_most_str_len,right_start_pos = sorted(long_string,reverse= True,key= lambda x:x[1])[0]
       left_most_str_len,start_pos = sorted(long_string,key= lambda x:x[1])[0]
    except IndexError:
        print "Low quality read"
        return False
    end_pos   =  right_start_pos + right_most_str_len + window_size -1
    #trim sequences and yield output
    trimed_seq        = record['seq'][start_pos:end_pos]
    trimed_qual_score = record['qual'][start_pos:end_pos]
    return {'ID': record['ID'], 'seq': trimed_seq.strip(), 'qual': trimed_qual_score.strip()}
####    
#===================================================================#   
def kmer_iter(kmer_size,trimmed_seq):
    for position in range(len(trimmed_seq) - kmer_size+1):
        each_kmer = trimmed_seq[position:(position + kmer_size)]
        #excluding the kmers having "N" bases
        if (each_kmer.count("N") == 0):
           yield sorted((each_kmer,each_kmer.translate(trans_table)[::-1]))[0]
    ####
####
#==================================================================#
def real_main():
    #Addiding commndline options
    usage  = "usage: %prog [input_file][kmer_size][options]"
    parser = optparse.OptionParser(usage)
    #scanning windowsize for fastq_trimmer  
    window_size = 20
    parser.add_option( '-w', \
                       "--windowsize", \
                       type    = "int", \
                       help    = "set the scanning windowsize for fastq_trimmer ; default: %i"%window_size, \
                       default = window_size )
    #minqual_score of bases to be included in the trimmed read
    qual_min = 20
    parser.add_option( '-q', \
                       "--minimumquality", \
                       type    = "int", \
                       help    = "set minqual_score of bases to be included in the trimmed read ; default: %i"%qual_min, \
                       default = qual_min)
    #standard deviation to vary color range in plot
    sd = 1.0
    parser.add_option( '-s', \
                       "--sd", \
                       type    = "int", \
                       help    = "set standarddeviation, to vary color range in the plot ; default: %f"%sd, \
                       default = sd)    
    #user defined prefix for outputfiles                  
    outputfile = "kmer"
    parser.add_option( '-o', \
                   "--outputfile", \
                   type    = "str", \
                   help    = "give the prefix for outputfiles, example: 24_mer ; default: %s"%outputfile, \
                   default = outputfile)                                
    # Parsing the arguments
    (options, args) = parser.parse_args()
    # Checking the user input
    if ( len(args) < 2 ):
        parser.error( "Incorrect number of arguments.  " + \
                      "View usage using --help option." )
    else:
        # Pulling the fileNames
        input_file = str(args[0])
        #check if inputfile exist or not
        if ( not os.path.isfile(input_file) ): parser.error( '%s can not be found'%input_file)
        kmer_size  = int(args[1])
    #determining the phred_scale
    #initializing an empty list to capture the numeric score of each character in qual_score line of each read
    character_list = []
    #initialize lop value to loop till 100 reads
    loop_counter   = 1
    for record in fastq_iterator(subprocess.Popen('cat %s| bunzip2 -c'%input_file,shell=True, stdout=subprocess.PIPE)):
        for each_char in record['qual']:
            character_list.append(ord(each_char))
            loop_counter+=1 
            #break the loop when 100th read is reached
            if(loop_counter == 100): break 
        ####
    ####
    #unique set of numeric qual_scores 
    character_set  = set(character_list)
    #capture the min_qual score
    min_qual_score = min(character_set)
    #capture the max qual score
    max_qual_score = max(character_set)
    #check the range of qual scores and determine the phred_scale
    if (min_qual_score >= 34 and max_qual_score <= 74): phred_offset = 33
    elif(min_qual_score >= 66 and max_qual_score <= 104): phred_offset = 64
    else: print("phred offset cannot be determined") 
    #calculating cutoff   
    cut_off              = options.windowsize*options.minimumquality
    #phredscoretable : converts phred ascii symbols to numeric values
    phred_score_table    = get_phred_score_table(phred_offset)
    #score table append 1 or 0 to the first base of the window, based on the average qual_score of bases in the window
    score_table          = get_score_table(options.windowsize,cut_off)
    #initializing empty list to measure readlength of each  read before trimming
    readlen_list         = []
    #initializing empty list to measure readlength of each trimed read
    trimmed_readlen_list = []
    #initializing empty dictionary for counting frequency of each kmer
    kmer_dict            =    {}
    for record in fastq_iterator(subprocess.Popen('cat %s| bunzip2 -c'%input_file,shell=True, stdout=subprocess.PIPE)):
        #measuring the readlength of read before trimming
        readlen_list.append(len(record['seq']))
        #each read is trimmed and parsed to kmer_iter
        trimmed_record  = fastq_trimmer(record,options.windowsize,phred_score_table,score_table)
        if not trimmed_record: continue
        #meauring the readlenth of each trimmed read
        trimmed_readlen_list.append(len(trimmed_record['seq']))
        #count of each kmer on both the strands of DNA
        for lexi_kmer in kmer_iter(kmer_size,trimmed_record['seq']):
            try:
                kmer_dict[lexi_kmer] +=1
            except KeyError:
               kmer_dict[lexi_kmer] = 1
            ####
        ####
   ####
   #calculating stats of reads before trimming
    average_read_length        = float(sum(readlen_list))/len(readlen_list)
    read_len_set               = set(readlen_list)
   #calculating the stats of trimmed reads
    average_trimmedread_length = float(sum(trimmed_readlen_list))/len(trimmed_readlen_list)
    trimmedread_len_set        = set(trimmed_readlen_list)
   #number of kmers in each fre bin are counted
    freq_dict = {}
    for kmer, kmer_count in kmer_dict.iteritems():
        try: 
           freq_dict[kmer_count] +=1
        except KeyError:
           freq_dict[kmer_count] = 1
        ####
    ####
    #eliminating the kmers whose frequency counts are 1 i.e. sequencing errors
    del freq_dict[1]
    #the base count of kmers in each freq_bin are counted
    base_count_dict = {}
    for freq,ncount in freq_dict.iteritems():
        base_count_dict[freq] = float(freq*ncount*kmer_size)
    ####
    #calculation of  "percent total sequence" of kmers in each fre_bin 
    total_base_count       = float(sum(base_count_dict.values()))
    total_seq_percent_dict = {}
    for freq,ncount in base_count_dict.iteritems():
         total_seq_percent_dict[freq] = 100.*float(ncount)/float(total_base_count)
    ####
    #average G and C base count of kmers in each freqbin is calculated
    GC_sum_count_dict ={}
    for kmer,kmer_count in kmer_dict.iteritems():
        try:
           GC_sum_count_dict[kmer_count] += kmer.count('G')+kmer.count('C')
        ####
        except KeyError:
            GC_sum_count_dict[kmer_count] = kmer.count('G')+kmer.count('C')
        ####
    ####
    #eliminating the GC_count of kmers whose frequency is 1
    del GC_sum_count_dict[1]
    #GC percent calculation
    GC_percent_dict  = dict((freq,(100.*(float(GC_sum_count_dict[freq])/float(freq_dict[freq]))/float(kmer_size))) for freq in GC_sum_count_dict)
    #calculating the mean of "GC_percent" (discrete random variable mean calculation)
    mean_calc_dict   = {}
    SD_calc_dict     = {}
    for freq in GC_percent_dict.iterkeys():
        mean_calc_dict[freq] =  (float(total_seq_percent_dict[freq])*float(GC_percent_dict[freq]))/100.
        #calculating Standard deviation 
        SD_calc_dict[freq] = float(mean_calc_dict[freq])*float(GC_percent_dict[freq])
    ####
    mean = sum(mean_calc_dict.values())
    SD   = math.sqrt(sum(SD_calc_dict.values()) - (mean*mean))
    #Calculating cumulative of  "percent total sequence"    
    cumulative_dict  = {}
    cum_sum          = 0
    sorted(total_seq_percent_dict)
    for freq,seq_percent in total_seq_percent_dict.iteritems():
        cum_sum               = cum_sum+seq_percent
        cumulative_dict[freq] = cum_sum
    ####
    #writing the final output matrix to file
    file = open("%s.dat"%options.outputfile,'w')
    for freq in total_seq_percent_dict.iterkeys():
        file.write('%i %i %2E %2E %2E \n' %(freq,freq_dict[freq],total_seq_percent_dict[freq],cumulative_dict[freq],GC_percent_dict[freq]))
    ####
    file.close()
    #writing mean and SD to file
    file = open("%s_mean_sd.txt"%options.outputfile,'w')
    file.write('%.2f' %(mean))
    file.write("\n")
    file.write('%.2f' %(SD))
    file.close()
    #defining the SD percintile to Zscore conversion table
    SD_Zscore_table_dict = {0.5:0.43,1.0:1,1.5:1.31,2.0:1.96,2.5:2.17,3.0 : 2.58}
    upperlimit = float(mean + (SD_Zscore_table_dict[options.sd]*SD))
    lowerlimit = float(mean - (SD_Zscore_table_dict[options.sd]*SD))
    #creating a gnu text file to generate gnuplot
    #####
    file = open("%s.txt"%options.outputfile,"w")
    file.write('set terminal png font arial 14 size 1000,800')
    file.write("\n")
    file.write('set output "%s.png"'%options.outputfile) 
    file.write("\n")
    file.write('set title "%s_plot"'%options.outputfile) 
    file.write("\n")
    file.write('set logscale xy 10')
    file.write("\n")
    file.write('set xtics auto')
    file.write("\n")
    file.write('set ytics auto')
    file.write("\n")
    file.write('set y2tics auto')
    file.write("\n")
    file.write('set xlabel "kmer frequency"')
    file.write("\n")
    file.write('set ylabel "Percent Total Sequence"')
    file.write("\n")
    file.write('set y2label "Cumulative"')
    file.write("\n")
    file.write('set boxwidth')
    file.write("\n")
    file.write('set style fill solid 1')
    file.write("\n")
    file.write("set palette rgb 33,13,10")
    file.write("\n")
    file.write("set cbrange[%.2f:%.2f]" %(upperlimit,lowerlimit))
    file.write("\n")
    file.write('plot "%s.dat" using 1:3:5  title "" with boxes palette axes x1y1, "%s.dat" using 1:4  title "" with line lt -1 lw -2 axes x1y2'%(options.outputfile,options.outputfile))
    file.close()
    ######
    #plot data using gnuplot
    subprocess.call("gnuplot   %s.txt"%options.outputfile,shell = True) 
    #writing the stas of raw and trimmed reads   to outputfile
    #########
    file = open("%s_inputfile_statistics"%options.outputfile,'w')
    file.write("Phred_scale is : %i"%phred_offset)
    file.write("\n")
    file.write("===========================   Raw Read Statistics  ======================== ")
    file.write("\n")
    file.write("average_read_length: %.2f"%average_read_length)
    file.write("\n")
    file.write("Minimum_read_length: %i"%min(read_len_set))
    file.write("\n")
    file.write("Maximum_read_length: %i"%max(read_len_set))
    file.write("\n")
    file.write("Number of reads: %i"%len(readlen_list))
    file.write("\n")
    file.write("===========================   Trimmed Read Statistics  ======================== ")
    file.write("\n")
    file.write("average_read_length: %.2f"%average_trimmedread_length)
    file.write("\n")
    file.write("Minimum_read_length: %i"%min(trimmedread_len_set))
    file.write("\n")
    file.write("Maximum_read_length: %i"%max(trimmedread_len_set))
    file.write("\n")
    file.write("Number of reads: %i"%len(trimmed_readlen_list))
    file.write("\n")
    file.close() 
    ######
    #comments to the user
    print "kmer counting done!!!!"
    print"files generated are located in current directory"
    print "the stats of the input file are written to '%s_inputfile_statistics'"%options.outputfile 
    print "data file is named as '%s.dat'"%options.outputfile 
    print "gnuscript file is named as '%s.txt'"%options.outputfile
    print "plot file is named as '%s.png' "%options.outputfile
####
#==================================================================#
if ( __name__ == '__main__' ):
    real_main()
