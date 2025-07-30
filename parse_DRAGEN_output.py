#!/usr/bin/python
_author__            = "m189786[yesesri_cherukuri]"
__Start_date__       = "11/30/2021"
__Lastmodified__     = "12/1/2021"
#===================================================================#
#modules  to be imported
from string import maketrans
import optparse
import subprocess
import os.path
from sys import stderr
import textwrap
from os import system
import pandas as pd
import sys
sys.path.insert(0,"/research/bsi/projects/staff_analysis/m189786/python_packages/")
import xlrd
#==================================================================#
def parse_PlateMaps(dir) :
    map_dict = {}
    for each_excel_file in os.listdir(dir) :
        plate = each_excel_file.strip().split("_")[2].replace(".xlsx","")
        wb    = xlrd.open_workbook("%s/%s"%(dir,each_excel_file.strip()))
        sheet = wb.sheet_by_index(0)
        for i in range(1,sheet.nrows) :
            row = sheet.row_values(i)
            sample = str(row[3]).encode('utf-8').split(".")[0]
            if(sample in map_dict.keys()) :
                stderr.write("%s : \n duplicate samples in plateMaps\n"%sample)
                continue
            assign_n = str(row[0]).encode('utf-8').split(".")[0]
            well   = str(row[1]).encode('utf-8').split(".")[0]
            map_dict[sample] = (plate,assign_n ,well )
        #####
    return(map_dict)
    ####
#==================================================================#
def parse_ClinicalData(file) :
    map_dict = {}
    wb = xlrd.open_workbook(file)
    sheet = wb.sheet_by_index(0)
    for i in range(1, sheet.nrows):
        temp_list = []
        for e in sheet.row_values(i) : temp_list.append(str(e).encode('utf-8'))
        sample = temp_list[0].split(".")[0]
        line = ""
        for i in range(1,len(temp_list)) :  line+="%s\t"%temp_list[i]
        map_dict[sample] = line
    ######
    return (map_dict)
#==================================================================#
def parse_CDC_VC(file) :
    map_dict = {}
    wb = xlrd.open_workbook(file)
    sheet = wb.sheet_by_index(0)
    for i in range(1, sheet.nrows):
        temp_list = []
        for e in sheet.row_values(i) : temp_list.append(e.encode('utf-8'))
        k = temp_list[1]
        line = "%s\t%s\t%s\t%s\t%s\t%s\t"%(temp_list[3],temp_list[4],temp_list[5],temp_list[6],temp_list[7],temp_list[8])
        map_dict[k] = line
    ######
    return (map_dict)
#==================================================================#
def real_main():
    # Addiding commndline options
    usage = "usage: %prog  [options]"
    parser = optparse.OptionParser(usage)
    # reference file name
    parser.add_option('-b',
                      "--basePath",
                      type="str",
                      help="location to current run analyses folder, where per sample dragen out folders were "
                           "located.  \
                           \n Example : /research/bsi/projects/PI/tertiary/Cerhan_James_cerhan/s301447.COVID-19"
                           "/processing/bs/output/211109_A00220_0264_BHJJHLDRXY/analyses/")

    parser.add_option('-r',
                      "--run_label",
                      type="str",
                      help="label the current run.  \
                           \n Example :NovaSeq23_211109")

    parser.add_option('-p',
                      "--plateMap",
                      type="str",
                      help="location to recently downloaded  platemap files from onedrive/sharepoint")
    parser.add_option('-c',
                      "--clinical_data_file",
                      type="str",
                      help="location to recently downloaded  clinical data files from onedrive/sharepoint")

    parser.add_option('-v',
                      "--CDC_VC_file",
                      type="str",
                      help="location to CDC variant classification file\n"
                           "location : /research/bsi/projects/PI/tertiary/Cerhan_James_cerhan/s301447.COVID-19"
                           "/CDC_variant_classifications/current/CDC_SARS_CoV-2_Variant_Classifications_as_of_10.06.2021.xlsx")
    parser.add_option('-s',
                      "--spikeN_count_file",
                      type="str",
                      help="location to spikeN_count file generated using spike_missing.sh\n")
    parser.add_option('-n',
                      "--contamination_check_file",
                      type="str",
                      help="location to contamination check file generated using contamination_check.pl\n")

    parser.add_option('-o',
                      "--out_dir",
                      type="str",
                      help="path where out files should be written to")

    # Parsing the arguments
    (options, args) = parser.parse_args()
    ##################################################
    #outfiles
    log_f            =  open("%s/%s.log.txt"%(options.out_dir,options.run_label),'w')
    DRAGEN_summary_f =  open("%s/%s.Illumina.Dragen_Summary.txt"%(options.out_dir,options.run_label),'w')

    ####################################################
    #Parse user Inputs
    basePath           = options.basePath
    PlateMaps_dict     = parse_PlateMaps(options.plateMap)
    clinical_map_dict  = parse_ClinicalData(options.clinical_data_file)
    CDC_VC_dict        = parse_CDC_VC(options.CDC_VC_file)
    spikeN_count_dict  = {}
    for line in open(options.spikeN_count_file) :
        x = line.strip().split()
        spikeN_count_dict[x[0]] = int(x[1])
    ######
    # location to NRCA table
    file = "/research/bsi/projects/PI/tertiary/Cerhan_James_cerhan/s301447.COVID-19" \
           "/m189786_Illumina_DRAGEN_summary_scripts/NRCA_table.txt"
    NRCA_map_dict = {}
    ####
    for line in open(file) :
        x    = line.strip().split()
        line = ""
        if(x[0].count("#") > 0 ) : continue
        line+="%s\t%s\t%s\t"%(x[0],x[1],x[2])
        #zip is key
        NRCA_map_dict[x[0]] = line
    #####
    zip_county_map_dict = {}
    for line in open("/research/bsi/projects/PI/tertiary/Cerhan_James_cerhan/s301447.COVID-19"
                     "/m189786_Illumina_DRAGEN_summary_scripts/Zip_county_062021.txt") :
        x = line.strip().split("\t")
        if(line.count("res_ratio") > 0 ) : continue
        # one zip code can have multiple county number's , pick the one with high residential ratio.
        try :
            zip_county_map_dict[x[0]].append( (x[1],float(x[4])) )
        except KeyError :
            zip_county_map_dict[x[0]] = [ ( x[1], float(x[4]) ) ]
    #####
    contamination_check_dict = {}
    for line in open(options.contamination_check_file) :
        if(line.count("Run") > 0 ) : continue
        x = line.strip().split("\t")
        try  :
            sample = x[1]
        except IndexError :
            continue
        line = ""
        for i in range(2, len(x)): line+="%s\t" % x[i]
        contamination_check_dict[sample] = line
    #####
    ########################################################
    #HARDCODED
    header    = "Run\tRunID1\tRunID2\tSample\tPlate\tDetected-HumanControl\tDetected-SARS-CoV2\tfracKmersCovered" \
                "-HumanControl\t" \
                "fracKmersCovered-SARS-CoV2\tfracTargetsCovered-Human Control\t" \
                "fracTargetsCovered-SARS-CoV2\tnKmerMatch-HumanControl\tnKmerMatch-SARS-CoV2\t" \
                "nTargetsDetected-HumanControl\tnTargetsDetected-SARS-CoV2\tnKmer-HumanControl\t" \
                "nKmer-SARS-CoV2\tnTarget-HumanControl\tnTarget-SARS-CoV2\t" \
                "lineage\tprobability\tconflict\tambiguity_score\tscorpio_call\tscorpio_support\t" \
                "scorpio_conflict\tversion\tpangolin_version\tpangoLEARN_version\tpango_version\tstatus\tnote\t" \
                "clade\tqc.overallScore\tqc.overallStatus\ttotalGaps\ttotalInsertions\ttotalMissing\ttotalMutations" \
                "\ttotalNonACGTNs\ttotalPcrPrimerChanges\tsubstitutions\tdeletions\tinsertions\tmissing\tnonACGTNs" \
                "\tpcrPrimerChanges\taaSubstitutions\ttotalAminoacidSubstitutions\taaDeletions\ttotalAminoacidDeletions\t" \
                "alignmentEnd\talignmentScore\talignmentStart\tqc.missingData.missingDataThreshold\tqc.missingData.score\t" \
                "qc.missingData.status\tqc.missingData.totalMissing\tqc.mixedSites.mixedSitesThreshold\tqc.mixedSites.score\t" \
                "qc.mixedSites.status\tqc.mixedSites.totalMixedSites\tqc.privateMutations.cutoff\tqc.privateMutations.excess\t" \
                "qc.privateMutations.score\tqc.privateMutations.status\tqc.privateMutations.total\t" \
                "qc.snpClusters.clusteredSNPs\tqc.snpClusters.score\tqc.snpClusters.status\tqc.snpClusters.totalSNPs" \
                "\terrors\tCOLLECTION.DT\tCOLLECT.CENTER.NAME\tCLINIC.NAME\tMRN\tAGE\tGENDER\tPT.ETHNICITY\tPT.RACE" \
                "\tPT.CITY\tPT.COUNTY\tPT.STATE\tZIP\tPT.COUNTRY\tSource_Soft.Order\tTarget.1.Value..ATTR." \
                "\tGROUP" \
                ".TEST.ID\tX.Days.Post.Pos.PCR.Dx\tTesting.Platform\tRegeneron.Set..In.Out.\tAssignNumber\tWell\t" \
                "Run_Sample\tSampleID0\tSampleID1\tSampleID2\tSampleID3\tRun_Plate\tSequencingCompletionDate\t" \
                "CDC_VC\tCDC_VC_clade\tCDC_VC_WHO.label\tCDC_VC_First.identified\tCDC_VC_Spike.Protein.Substitutions" \
                "\tCDC_VC_Attributes\tSpikeN_count\tcountyfips\tcountyname\tcountystate\tMatched.Variants\tMatched" \
                ".HighQual.Variants\tHet.HQ.Var.Count\tPercent.Het\tTraining.Set.Var.Match.Count\tContamination.Status\tSpikeMissing\tQuality\n"
    ########################################################
    DRAGEN_summary_f.write(header)
    for each_sample_folder  in os.listdir(options.basePath) :
        each_sample_folder = each_sample_folder.strip()
        sample = each_sample_folder.split("_")[0].replace("-","_")
        try:
            plate    = PlateMaps_dict[sample.replace("-","_")][0]
            assign_n = PlateMaps_dict[sample.replace("-","_")][1]
            well     = PlateMaps_dict[sample.replace("-","_")][2]
        except KeyError:
            plate     = "NA"
            assign_n  = "NA"
            well      = "NA"
        #####
        try :
            spikeN_count = spikeN_count_dict[sample]
        except KeyError :
            spikeN_count = "NA"
        #####
        try :
            contamination_check_line = contamination_check_dict[sample]
        except KeyError :
            contamination_check_line = "NA\tNA\tNA\tNA\tNA\tNA\t"
        #####
        runID = "%s"%options.run_label
        runID_1 = runID.split("_")[0]
        runID_2 = runID.split("_")[1]
        Info_line = ""
        Info_line+="%s\t%s\t%s\t%s\t%s\t"%(runID,runID_1,runID_2,sample,plate)
        ####
        # Parse DRAGEN outfile
        #Summary file
        summary_f = "%s/%s/Summary.tsv"%(basePath,each_sample_folder)
        if (os.path.isfile(summary_f)):
            n=0
            for line in open(summary_f) :
                n+=1
                if(n==1) : continue
                else :
                    x = line.strip().split("\t")
                    for i in range(1, len(x)): Info_line+="%s\t"%x[i]
        else :
            stderr.write("No Sumaary.tsv file %s\n"%each_sample_folder)
            continue
        ##### end of summary file.
        #Lineage file
        CDC_VC_line = "NA\tNA\tNA\tNA\tNA\tNA\t"
        lineage_f = "%s/%s/consensus/%s.lineage_report.csv" % (basePath, each_sample_folder.strip(), sample)
        if( os.path.isfile(lineage_f) ) :
            n=0
            for line in open(lineage_f):
                n += 1
                if (n == 1):continue
                else:
                    x = line.strip().split(",")
                    lineage   = x[1]
                    Info_line+="%s\tNA\t"%(lineage)
                    for i in range(2,len(x)):Info_line+="%s\t"%x[i]
                    try :
                        CDC_VC_line = CDC_VC_dict[lineage]
                    except KeyError :
                        CDC_VC_line = "NA\tNA\tNA\tNA\tNA\tNA\t"
        else :
            log_f.write("No lineage file %s.lineage_report.csv\n"%sample)
            Info_line+="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
         ##### end of lineage file
        #clade file
        clade_f = "%s/%s/consensus/%s.nextclade.tsv" % (basePath, each_sample_folder.strip(), sample)
        if(os.path.isfile(clade_f)) :
            n = 0
            for line in open(clade_f):
                n += 1
                if (n == 1): continue
                else:
                    x =  line.strip().split("\t")
                    for i in range(1 , len(x)): Info_line+="%s\t"%x[i]
                    if(len(x) < 41 ) : Info_line+="None\t"
        else:
            log_f.write("No clade file %s.nextclade.tsv\n"%sample)
            Info_line+="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" \
                       "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
        ##### end of clade file
        ##### End of Parsing DRAGEN files for current sample
        # Merge Clinical Information
        try:
            clinical_Info = clinical_map_dict[sample]
        except KeyError:
            clinical_Info = "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"
        #####
        Info_line+="%s"%(clinical_Info)
        Info_line+="%s\t%s\t"%(assign_n,well)
        # additional sample Infor
        Run_sample = "%s_%s"%(runID, sample)
        SampleID0 = runID
        SampleID1 = sample
        if (len(sample.split("_")) > 1):
            SampleID2 = sample.split("_")[0]
            SampleID3 = sample.split("_")[1]
        else:
            SampleID2 = "NA"
            SampleID3 = "NA"
        Run_Plate = "%s_%s" % (runID, plate)
        SequencingCompletionDate = runID.split("_")[1]
        Info_line+="%s\t%s\t%s\t%s\t%s\t%s\t%s\t"%(Run_sample,SampleID0,SampleID1,SampleID2,SampleID3,Run_Plate,
                                   SequencingCompletionDate)
        # Merge CDC varaint classification.
        Info_line+="%s"%(CDC_VC_line)
        # Merge SpikeN.
        Info_line+="%s\t"%(spikeN_count)
        # Merge NRCA table.
        try :
            zip_code   = clinical_Info.split("\t")[11].split(".")[0]
            temp_list  = sorted( zip_county_map_dict[zip_code] , key = lambda x: x[1] )
            county_num =  temp_list[-1][0]
            NRCA_line  =  NRCA_map_dict[county_num]
        except KeyError :
            NRCA_line = "NA\tNA\tNA\t"
        #######
        Info_line+="%s"%(NRCA_line)
        # Merge Contamination data.
        Info_line+="%s"%contamination_check_line
        #SpineNmissing and Quality
        #Threshold defined in Pete's script /research/bsi/projects/PI/tertiary/Cerhan_James_cerhan/s301447.COVID-19/processing/bs/summary/COVID_DataClean_DataQC_1011-11_threechanges_and retitle_2021__mforge.Rmd
        qc_overallStatus = Info_line.split("\t")[34]

        Quality = "NA"
        try :
            totalMissing     = int(Info_line.split("\t")[37])
            spikeN_count     = int(spikeN_count)

            if (qc_overallStatus == 'good' and spikeN_count == 0): Quality = "PASS"
            if (qc_overallStatus == 'good' and spikeN_count == 0 and totalMissing > 0): Quality = "High"

            if (qc_overallStatus != 'bad' and spikeN_count > 0 and spikeN_count < 168): Quality = "Medium"
            if (qc_overallStatus != 'bad' and qc_overallStatus != 'good' and spikeN_count < 168): Quality = "Medium"

            if (qc_overallStatus != 'bad' and spikeN_count >= 168): Quality = "Low"

            if (qc_overallStatus == 'bad' ): Quality =  "Fail"

            if ( qc_overallStatus != 'bad' and qc_overallStatus != 'good' and qc_overallStatus != 'mediocre'):  Quality = "Fail"

        except ValueError :
            Quality = "Fail"
        #####
        Info_line+="%s\t%s"%(spikeN_count, Quality)
        #write merged line to outfile
        if not ( len(Info_line.split("\t")) == len(header.split("\t"))) :
            print (  len(Info_line.split("\t")) )
            stderr.write("Total number of columns not correct , check")
            #continue
        #####
        DRAGEN_summary_f.write("%s\n" % Info_line)
    #### End of parsing data and merging for all samples.

##### end of main function
#==================================================================#
if ( __name__ == '__main__' ):
    real_main()
