#!/usr/bin/env python

import math
import pandas as pd
from optparse import OptionParser
import os
import sys



class input_data_pretreatment:
    def __init__( self , parameters_pretreatment ):
        self.chr_list = [ ]
        self.chr2maxCoordinate_dict = { }
        self.chr2lengthSetting_dict = { }
        self.chr2lines_dict = { }
        self.chr2outputfile_dict = { }
        if parameters_pretreatment['single_or_paired_setting'] == 's':
            self.do_pretreatment( parameters_pretreatment )
        elif parameters_pretreatment['single_or_paired_setting'] == 'p':
            self.do_pretreatment_for_paired_end( parameters_pretreatment )
        self.write_down_summary( parameters_pretreatment )

    def do_pretreatment( self , parameters_pretreatment ):
        chr2outputTEXT_dict = { }
        count_raw_line = 0 
        print('Read raw data: ' , parameters_pretreatment['Raw_input_file'] )
        print('Path to record preliminary data: ' , parameters_pretreatment['Output_path'] )
        inputTEXT = open( parameters_pretreatment['Raw_input_file'] , 'r' ) 
        while True:
            try:
                original_line = next( inputTEXT )
                count_raw_line = count_raw_line + 1  
                sys.stdout.write('\r    Raw line: '+str(count_raw_line) )
                line2list = original_line.split('\n')[0].split('\r')[0].split('\t')
                if line2list[ 5 ] == '+':
                    right_coordinate = int( line2list[1] ) + 150
                elif line2list[ 5 ] == '-':
                    right_coordinate = int( line2list[2] )
                chromosome = line2list[ 0 ]
                if chromosome not in self.chr2outputfile_dict.keys():
                    self.chr2outputfile_dict[ chromosome ] = parameters_pretreatment['Output_path']+chromosome+'.bed'
                    chr2outputTEXT_dict[ chromosome ] = open( self.chr2outputfile_dict[ chromosome ] , 'w' )
                    self.chr2maxCoordinate_dict[ chromosome ] = 0
                    self.chr2lines_dict[ chromosome ] = 0
                    self.chr_list.append( chromosome )
                elif chromosome in self.chr2outputfile_dict.keys():
                    pass
                chr2outputTEXT_dict[ chromosome ].write( original_line )
                self.chr2lines_dict[ chromosome ] = self.chr2lines_dict[ chromosome ] + 1
                if right_coordinate > self.chr2maxCoordinate_dict[ chromosome ]:
                    self.chr2maxCoordinate_dict[ chromosome ] = right_coordinate
                else:
                    pass
                sys.stdout.flush()
            except StopIteration:
                break
        inputTEXT.close() 
        for chromosome in self.chr_list:
            chr2outputTEXT_dict[ chromosome ].close()
            self.chr2lengthSetting_dict[ chromosome ] = self.chr2maxCoordinate_dict[ chromosome ]

    def do_pretreatment_for_paired_end( self , parameters_pretreatment ):
        chr2outputTEXT_dict = { }
        count_raw_line = 0 
        print('Read raw data: ' , parameters_pretreatment['Raw_input_file'] )
        print('Path to record preliminary data: ' , parameters_pretreatment['Output_path'] )
        inputTEXT = open( parameters_pretreatment['Raw_input_file'] , 'r' ) 
        while True:
            try:
                original_line = next( inputTEXT )
                count_raw_line = count_raw_line + 1  
                sys.stdout.write('\r    Raw line: '+str(count_raw_line) )
                line2list = original_line.split('\n')[0].split('\r')[0].split('\t')
                right_coordinate = int( line2list[2] )
                chromosome = line2list[ 0 ]
                if chromosome not in self.chr2outputfile_dict.keys():
                    self.chr2outputfile_dict[ chromosome ] = parameters_pretreatment['Output_path']+chromosome+'.bed'
                    chr2outputTEXT_dict[ chromosome ] = open( self.chr2outputfile_dict[ chromosome ] , 'w' )
                    self.chr2maxCoordinate_dict[ chromosome ] = 0
                    self.chr2lines_dict[ chromosome ] = 0
                    self.chr_list.append( chromosome )
                elif chromosome in self.chr2outputfile_dict.keys():
                    pass
 
                chr2outputTEXT_dict[ chromosome ].write( original_line )
                self.chr2lines_dict[ chromosome ] = self.chr2lines_dict[ chromosome ] + 1
                if right_coordinate > self.chr2maxCoordinate_dict[ chromosome ]:
                    self.chr2maxCoordinate_dict[ chromosome ] = right_coordinate
                else:
                    pass
                sys.stdout.flush()
            except StopIteration:
                break
        inputTEXT.close() 
        for chromosome in self.chr_list:
            chr2outputTEXT_dict[ chromosome ].close()
            self.chr2lengthSetting_dict[ chromosome ] = self.chr2maxCoordinate_dict[ chromosome ]

    def write_down_summary( self , parameters_pretreatment ):
        print('Write down summary: ' , parameters_pretreatment['OutputFile_summary'] )
        outputTEXT = open( parameters_pretreatment['OutputFile_summary'] , 'w' )
        outputTEXT.write( ('\t').join( [ 'Chromosome' , 'Tags' , 'Max_coordinate' , 'Recorded_chromosome_length' ] ) + '\n' )
        for chromosome in sorted( self.chr_list ):
            outputTEXT.write( ('\t').join( [ chromosome , str(self.chr2lines_dict[ chromosome ]) , str(self.chr2maxCoordinate_dict[ chromosome ]) , str(self.chr2lengthSetting_dict[ chromosome ]) ] ) + '\n' )
        outputTEXT.close()


class NucDetect:
    def __init__(self,FilesParameters,ConvolutionParameters,threshold):
        self.inputfile=FilesParameters['NucPosition_in']
        self.celltype=self.inputfile.split('/')[-1].split('Nucleosome')[0]
        self.chrlength=chromosome_length[self.chromosome]
        self.chrRecordingListLength=0
        self.outputfile_wiggle=self.outputfile_path+self.inputfile.split('/')[-1][0]+'_'+self.chromosome+FilesParameters['outputfile_suffix']+'.like_wig'
        self.outputfile_nucleosome=self.outputfile_path+self.inputfile.split('/')[-1][0]+'_'+self.chromosome+FilesParameters['outputfile_suffix']+'.like_bed'
        self.ConvolutionParameters=ConvolutionParameters
        self.threshold=threshold
        self.score_list=[]
        self.tag_list = [ ]
        self.Gaussian_list=[]
        self.FDoG_list=[]
        self.LoG_list=[]
        self.TDoG_list=[]
        self.score_table=[]  
        self.secondary_LoG_list=[]

    def convolution(x,Gaussian,First_Derivative_of_Gaussian,LoG,Third_Derivative_of_Gaussian,sigma,times_of_sigma):  
        y=0
        FDoG_y=0
        LoG_y=0
        TDoG_y=0
        if x>=sigma*times_of_sigma and x<=len(self.score_list)-sigma*times_of_sigma-1:
            for n in range(x-sigma*times_of_sigma,x+sigma*times_of_sigma+1):
                y=y+self.score_list[n]*Gaussian[n-(x-sigma*times_of_sigma)]
                FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[n-(x-sigma*times_of_sigma)]
                LoG_y=LoG_y+self.score_list[n]*LoG[n-(x-sigma*times_of_sigma)]
                TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[n-(x-sigma*times_of_sigma)]
        elif x<sigma*times_of_sigma:
            for n in range(0,x+sigma*times_of_sigma+1):
                y=y+self.score_list[n]*Gaussian[sigma*times_of_sigma-x+n]
                FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
                LoG_y=LoG_y+self.score_list[n]*LoG[sigma*times_of_sigma-x+n]
                TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
        elif x>len(self.score_list)-sigma*times_of_sigma-1:
            for n in range(x-sigma*times_of_sigma,len(self.score_list)):
                y=y+self.score_list[n]*Gaussian[sigma*times_of_sigma-x+n]
                FDoG_y=FDoG_y+self.score_list[n]*First_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
                LoG_y=LoG_y+self.score_list[n]*LoG[sigma*times_of_sigma-x+n]
                TDoG_y=TDoG_y+self.score_list[n]*Third_Derivative_of_Gaussian[sigma*times_of_sigma-x+n]
        return y,FDoG_y,LoG_y,TDoG_y


    def extremum_detection(self):
        print('Detecting extremum for peaks and valleys ......')
        self.peak_dict={}
        self.peak_details_dict={}
        one_peak=['left_valley_ending_coordinate','max_left_beginning_coordinate','max_right_ending_coordinate','right_valley_beginning_coordinate',
                  'left_valley_ending_index','max_left_beginning_index','max_right_ending_index','right_valley_beginning_index',
                  'serial_number']
        one_peak_finished='none'
        holding_on=0
        above_or_below=0
        serial_number=0
        for t in range(1,len(self.score_table)):
            if self.FDoG_list[t]>0:    
                if above_or_below==0:  
                    one_peak[0]=self.score_table[t][0]
                    one_peak[4]=t
                elif above_or_below==1:    
                    one_peak[1]='max_left_beginning_coordinate'
                    one_peak[5]='max_left_beginning_index'
                elif above_or_below==-1:    
                    if type(one_peak[3])!=int or type(one_peak[7])!=int:
                        one_peak[3]=self.score_table[t-1][0]
                        one_peak[7]=t-1
                    one_peak_finished=one_peak[:]
                    one_peak=['left_valley_ending_coordinate','max_left_beginning_coordinate','max_right_ending_coordinate','right_valley_beginning_coordinate','left_valley_ending_index','max_left_beginning_index','max_right_ending_index','right_valley_beginning_index','serial_number']
                    one_peak[0]=self.score_table[t][0]
                    one_peak[4]=t
                above_or_below=1
                holding_on=0
            elif self.FDoG_list[t]<0:    
                if above_or_below==-1:
                    one_peak[3]='right_valley_beginning_coordinate'
                    one_peak[7]='right_valley_beginning_index'
                elif above_or_below==1:
                    one_peak[2]=self.score_table[t-holding_on][0]
                    one_peak[6]=t-holding_on
                    if type(one_peak[1])!=int or type(one_peak[5])!=int:
                        one_peak[1]=self.score_table[t-max(1,holding_on)][0]
                        one_peak[5]=t-max(1,holding_on)
                above_or_below=-1
                holding_on=0
            elif self.FDoG_list[t]==0:    
                holding_on=holding_on+1
                if above_or_below==1 and self.FDoG_list[t-1]>0:
                    one_peak[1]=self.score_table[t][0]
                    one_peak[5]=t
                elif above_or_below==-1 and self.FDoG_list[t-1]<0:
                    one_peak[3]=self.score_table[t-1][0]
                    one_peak[7]=t-1
            if t==len(self.score_table)-1:
                if self.FDoG_list[t]<0:                                     
                    if above_or_below==-1:                                  
                        if type(one_peak[3])!=int or type(one_peak[7])!=int:  
                            one_peak[3]=self.score_table[t][0]                
                            one_peak[7]=t                                    
                one_peak_finished=one_peak[:]
            if one_peak_finished!='none' and type(one_peak_finished)==list and type(one_peak_finished[0])==int and type(one_peak_finished[1])==int and type(one_peak_finished[2])==int and type(one_peak_finished[3])==int and type(one_peak_finished[4])==int and type(one_peak_finished[5])==int and type(one_peak_finished[6])==int and type(one_peak_finished[7])==int:
                serial_number=serial_number+1
                one_peak_finished[8]=serial_number
                self.peak_dict[serial_number]=one_peak_finished 
                for p in range(one_peak_finished[4],one_peak_finished[7]+1):
                    self.peak_details_dict[p]=one_peak_finished    
                one_peak_finished='none'
        print('\tFinished ==> Peak numbers:  ',len(self.peak_dict))


    def inflection_pairs_detection(self):
        print('Detecting inflection pairs for nucleosome candidates ......')
        self.inflection_pairs_list=[]
        one_pair_of_inflection=['unknown','unknown','n1','n2','nindex']
        holding_on=0
        above_or_below=0
        nindex=0
        for t in range(1,len(self.score_table)):
            if self.score_table[t][3]>0:
                holding_on=0
                if above_or_below==1:
                    pass
                elif above_or_below!=1:
                    if type(one_pair_of_inflection[0])==int and type(one_pair_of_inflection[2])==int:
                        one_pair_of_inflection[1]=self.score_table[t-1][0]
                        one_pair_of_inflection[3]=t-1
                    else:
                        pass
                    above_or_below=1
            elif self.score_table[t][3]<0:
                holding_on=0
                if above_or_below==-1:
                    pass
                elif above_or_below!=-1:
                    one_pair_of_inflection[0]=self.score_table[t][0]
                    one_pair_of_inflection[2]=t
                    one_pair_of_inflection[1]='unknown'
                    one_pair_of_inflection[3]='n2'
                    above_or_below=-1
            elif self.score_table[t][3]==0:
                holding_on=holding_on+1
            if type(one_pair_of_inflection[0])==int and type(one_pair_of_inflection[1])==int and type(one_pair_of_inflection[2])==int and type(one_pair_of_inflection[3])==int:
                nindex=nindex+1
                one_pair_of_inflection[4]=nindex
                self.inflection_pairs_list.append(one_pair_of_inflection)
                one_pair_of_inflection=['unknown','unknown','n1','n2','nindex']
        print('\tFinished ==> Number of nucleosome candidates:  ',len(self.inflection_pairs_list))


    def inflection_pairs_midpoint_detection(self):
        print('Detecting midpoint of every inflection pairs ......')
        midpoints_candidates_list=[]   
        midpoints_candidates_dict={}    
        midpoints_candidates_details_dict={}   
        midpoints_candidates_index=0
        one_midpoint=['unknown','unknown','m1','m2','mindex']
        holding_on=0
        above_or_below=0
        for t in range(1,len(self.score_table)):
            if self.TDoG_list[t]>0:
                if above_or_below==1:
                    pass
                elif above_or_below==-1:
                    one_midpoint[1]=self.score_table[t][0]
                    one_midpoint[3]=t
                    one_midpoint[0]=self.score_table[t-max(1,holding_on)][0]
                    one_midpoint[2]=t-max(1,holding_on)
                above_or_below=1
                holding_on=0
            elif self.TDoG_list[t]<0:
                above_or_below=-1
                holding_on=0
            elif self.TDoG_list[t]==0:
                holding_on=holding_on+1
            if type(one_midpoint[0])==int and type(one_midpoint[1])==int and type(one_midpoint[2])==int and type(one_midpoint[3])==int:
                midpoints_candidates_index=midpoints_candidates_index+1
                one_midpoint[4]=midpoints_candidates_index
                midpoints_candidates_list.append(one_midpoint)    ### Recording
                midpoints_candidates_dict[midpoints_candidates_index]=one_midpoint   ### Recording
                for m in range(one_midpoint[2],one_midpoint[3]+1):
                    midpoints_candidates_details_dict[m]=one_midpoint    ### Recording
                one_midpoint=['unknown','unknown','m1','m2','mindex']
        print('...... The detection for all midpoint candidates of every inflection pairs is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        self.inflection_pairs_list_with_midpoint=[]
        self.inflection_pairs_dict_with_midpoint={}
        real_midpoints_index=0
        for a_pair_of_inflection in self.inflection_pairs_list:
            a_pair_of_inflection_with_midpoint=a_pair_of_inflection[:]
            mapped_midpoint_mindex=[]
            for point in range(a_pair_of_inflection[2],a_pair_of_inflection[3]+1):
                if point in midpoints_candidates_details_dict.keys():
                    if midpoints_candidates_details_dict[point][4] not in mapped_midpoint_mindex:
                        mapped_midpoint_mindex.append(midpoints_candidates_details_dict[point][4])
                    else:
                        pass
            if len(mapped_midpoint_mindex)>0:
                real_midpoints_index=real_midpoints_index+len(mapped_midpoint_mindex)
                for pointt in mapped_midpoint_mindex:
                    a_pair_of_inflection_with_midpoint.append(midpoints_candidates_dict[pointt])    ### 给一对拐点注明其“中点”
            else:
                pass
            self.inflection_pairs_list_with_midpoint.append(a_pair_of_inflection_with_midpoint)
            self.inflection_pairs_dict_with_midpoint[a_pair_of_inflection_with_midpoint[4]]=a_pair_of_inflection_with_midpoint
        print('\tChecking ==> Total inflection pairs with midpoint:  ',len(self.inflection_pairs_list_with_midpoint))

    def preliminary_nucleosome_position(self):
        print('Beginning preliminary main nucleosome and shoulder positioning ......')
        integrated_pairs_counting=0    
        self.unintegrated_inflection_pairs_list=[]
        self.nucleosome_list=[]
        self.nucleosome_dict={}
        self.shoulder_list=[]
        self.shoulder_dict={}
        for one_pair_of_inflection in self.inflection_pairs_list_with_midpoint:
            peak_checking=0    
            valley_left_checking=0    
            valley_right_checking=0    
            peak_index_list=[]
            for n in range(one_pair_of_inflection[2],one_pair_of_inflection[3]+1):
                if n in self.peak_details_dict.keys():
                    if self.peak_details_dict[n][8] not in peak_index_list:
                        peak_index_list.append(self.peak_details_dict[n][8])
                    else:
                        pass
            if len(peak_index_list)==1:
                peak_index=peak_index_list[0]
                self.peak_dict[peak_index].append(one_pair_of_inflection)    
                integrated_pairs_counting=integrated_pairs_counting+1
                max_beginning=self.peak_dict[peak_index][5]
                max_beginning_coordinate=self.peak_dict[peak_index][1]
                max_ending=self.peak_dict[peak_index][6]
                max_ending_coordinate=self.peak_dict[peak_index][2]
                peak_beginning=self.peak_dict[peak_index][4]
                peak_beginning_coordinate=self.peak_dict[peak_index][0]
                peak_ending=self.peak_dict[peak_index][7]
                peak_ending_coordinate=self.peak_dict[peak_index][3]
                if (self.peak_dict[peak_index][5]>=one_pair_of_inflection[2] and self.peak_dict[peak_index][5]<=one_pair_of_inflection[3]) or (self.peak_dict[peak_index][6]>=one_pair_of_inflection[2] and self.peak_dict[peak_index][6]<=one_pair_of_inflection[3]):    ####peak_dict[peak_index][5]>=one_pair_of_inflection[2] and peak_dict[peak_index][6]<=one_pair_of_inflection[3]:
                    peak_checking=1    
                if one_pair_of_inflection[2]>=self.peak_dict[peak_index][4]:
                    valley_left_checking=1    
                if one_pair_of_inflection[3]<=self.peak_dict[peak_index][7]:
                    valley_right_checking=1   
            else:
                self.unintegrated_inflection_pairs_list.append(one_pair_of_inflection)
            if peak_checking==1 and valley_left_checking==1 and valley_right_checking==1:   
                self.nucleosome_list.append([one_pair_of_inflection[0],one_pair_of_inflection[1],one_pair_of_inflection[2],one_pair_of_inflection[3],one_pair_of_inflection[4],'max(original_score_fragment)','area',[peak_beginning_coordinate,max_beginning_coordinate,max_ending_coordinate,peak_ending_coordinate,peak_beginning,max_beginning,max_ending,peak_ending,peak_index]])
                self.nucleosome_dict[peak_index]=[one_pair_of_inflection[0],one_pair_of_inflection[1],one_pair_of_inflection[2],one_pair_of_inflection[3],one_pair_of_inflection[4],'max(original_score_fragment)','area',[peak_beginning_coordinate,max_beginning_coordinate,max_ending_coordinate,peak_ending_coordinate,peak_beginning,max_beginning,max_ending,peak_ending,peak_index]]
            elif peak_checking==0 and valley_left_checking==1 and valley_right_checking==1:   
                self.shoulder_list.append([one_pair_of_inflection[0],one_pair_of_inflection[1],one_pair_of_inflection[2],one_pair_of_inflection[3],one_pair_of_inflection[4],'max(original_score_fragment)','area',[peak_beginning_coordinate,max_beginning_coordinate,max_ending_coordinate,peak_ending_coordinate,peak_beginning,max_beginning,max_ending,peak_ending,peak_index]])
                self.shoulder_dict[one_pair_of_inflection[4]]=[one_pair_of_inflection[0],one_pair_of_inflection[1],one_pair_of_inflection[2],one_pair_of_inflection[3],one_pair_of_inflection[4],'max(original_score_fragment)','area',[peak_beginning_coordinate,max_beginning_coordinate,max_ending_coordinate,peak_ending_coordinate,peak_beginning,max_beginning,max_ending,peak_ending,peak_index]]
            else:
                pass
        print('\tChecking ==> Total main nucleosome:  ',len(self.nucleosome_list))
        print('\tChecking ==> Number of shoulders:  ',len(self.shoulder_list))


def count_chromosome_length( input_file ):
    max_coordinate = 0
    inputTEXT = open( input_file , 'r' )
    while True:
        try:
            line2list = next( inputTEXT ).split('\n')[0].split('\r')[0].split('\t')
            if line2list[ 5 ] == '+':
                right_coordinate = int( line2list[1] ) + 150
            elif line2list[ 5 ] == '-':
                right_coordinate = int( line2list[2] )
            if right_coordinate > max_coordinate:
                max_coordinate = right_coordinate
            else:
                pass
        except StopIteration:
            break
    inputTEXT.close()
    return max_coordinate

def Merge_overall_results( parameters_Merge_results ):
    outputTEXT = open( parameters_Merge_results['OutputFile'] , 'w' )
    total_nucleosome____dict = { }
    headline2 = ''
    for chromosome in sorted( parameters_Merge_results['Chromosome_list'] ):
        InputFile = parameters_Merge_results['Input_Files'] + '_' + chromosome + '.like_bed'
        inputTEXT = open( InputFile , 'r' ) 
        headline1 = next( inputTEXT ) 
        outputTEXT.write( headline1 ) 
        headline2 = next( inputTEXT ) 
        inputTEXT.close() 
    outputTEXT.write( headline2 )    
    for chromosome in sorted( parameters_Merge_results['Chromosome_list'] ):
        InputFile = parameters_Merge_results['Input_Files'] + '_' + chromosome + '.like_bed'
        print('Read from: ' , InputFile )
        inputTEXT = open( InputFile , 'r' ) 
        next( inputTEXT )
        next( inputTEXT ) 
        while True:
            try:
                original_line = next( inputTEXT )
                outputTEXT.write( original_line )    
            except StopIteration:
                break
        inputTEXT.close() 
    outputTEXT.close()
    print('...... Collecting nucleosome results ')
    print('\n')


FilesParameters={'NucPosition_in':'none',
                 'Chosen_Chromosome_Abbreviation':'none',
                 'Chosen_Chromosome_Length':'none',
                 'outputfile_like_bed':'',
                 'outputfile_like_wig':'',
                 }

if __name__=='__main__':
    
    print('\n\nProgram to detect nucleosomes from MNase-seq data.')
    print('\nThe program should run in Python3.\n')
    usage='usage: NucDetect.py [options]'
    parser = OptionParser( usage = '%prog  -i /path/inputfile  -o /path/outputfile  -c chromosome_name  -l chromosome_length  --s_p single_or_paired_end_data' , version = '%prog Version:1.2.1' )
    parser.add_option('-i', '--input',      dest='input_file',                        type='string', help='"/path/filename" Path of input ')
    parser.add_option('-o', '--output',     dest='output_file',                       type='string', help='"/path/filename" gather the detected nucleosomes on every chromosome. Note that a path "/path/filename/" or "/path/filename_[ChromosomeName]/" will be built to record the preliminary and intermediate data.')
    parser.add_option('-c', '--chrname',    dest='chromosome_name',                   type='string', help='Specify the name (or abbreviation) of the chromosome,')
    parser.add_option('-l', '--chrlength',  dest='chromosome_length',                 type='int',    help='The length of the chromosome specified by parameter "-c" or "--chrname". ')
       parser.add_option('--pe_max',           dest='superior_limit_of_paired_end_tags', type='int',    default=200 , help='The superior limit of the length of paired-end tags.')
    parser.add_option('--pe_min',           dest='inferior_limit_of_paired_end_tags', type='int',    default=100 , help='The inferior limit of the length of paired-end tags. ')
    (options, args) = parser.parse_args( )
    if options.input_file == None:
        parser.error('-h for help or provide the input file name!')
    else:
        pass
    if options.output_file == None:
        parser.error('-h for help or provide the output file name!')
    else:
        pass
    ConvolutionParameters=Gaussian_profile(PreliminaryConvolutionParameters)

    if options.chromosome_name != None:
        Intermediate_Path = options.output_file+'_'+options.chromosome_name+'/' 
        if os.path.exists( Intermediate_Path ):
            pass
        else:
            os.makedirs( Intermediate_Path )
        parameters_pretreatment_for_single_chromosome = { 'Raw_input_file': options.input_file,
                                                          'Output_path':    Intermediate_Path,
                                                          'OutputFile_summary': Intermediate_Path+'InputData_Summary.txt',
                                                          'chromosome':   options.chromosome_name ,
                                                          'single_or_paired_setting': options.single_or_paired_end,
                                                          }
        if options.chromosome_length == None:
            parameters_pretreatment_for_single_chromosome[ 'length' ] = 'none'
        elif options.chromosome_length != None:
            parameters_pretreatment_for_single_chromosome[ 'length' ] = options.chromosome_length
        InputData_Summary_for_single_chromosome = input_data_pretreatment_for_single_chromosome( parameters_pretreatment_for_single_chromosome ) 
        if ( options.chromosome_name in InputData_Summary_for_single_chromosome.chr_recording_list ) and ( InputData_Summary_for_single_chromosome.chromosome_lines > 0 ):
            FilesParameters[ 'Chosen_Chromosome_Abbreviation' ] = options.chromosome_name
            FilesParameters[ 'NucPosition_in' ] = InputData_Summary_for_single_chromosome.OutputFile_Tags
            FilesParameters[ 'outputfile_like_bed' ] = options.output_file + '_' + options.chromosome_name + '.like_bed'
            FilesParameters[ 'outputfile_like_wig' ] = options.output_file + '_' + options.chromosome_name + '.like_wig'
            if options.chromosome_length == None:
                FilesParameters[ 'Chosen_Chromosome_Length' ] = InputData_Summary_for_single_chromosome.chromosome_lengthSetting
            elif options.chromosome_length != None:
                FilesParameters[ 'Chosen_Chromosome_Length' ] = options.chromosome_length
            Positioning_improvement( FilesParameters , ConvolutionParameters , threshold )
        else:
            print('Please provide correct chromosome name, or choose chromosome among: ' , (',').join(sorted(InputData_Summary_for_single_chromosome.chr_recording_list)) )
    elif options.chromosome_name == None:
        Intermediate_Path = options.output_file+'/' 
        if os.path.exists( Intermediate_Path ):
            pass
        else:
            os.makedirs( Intermediate_Path )
        print('Do data pretreatment ......')
        parameters_pretreatment = { 'Raw_input_file': options.input_file,
                                    'Output_path':    Intermediate_Path,
                                    'OutputFile_summary': Intermediate_Path+'InputData_Summary.txt',
                                    'single_or_paired_setting': options.single_or_paired_end,
                                    }
        InputData_Summary = input_data_pretreatment( parameters_pretreatment )
        for each_chromosome in sorted(InputData_Summary.chr_list):
            FilesParameters[ 'Chosen_Chromosome_Abbreviation' ] = each_chromosome
            FilesParameters[ 'NucPosition_in' ] = InputData_Summary.chr2outputfile_dict[ each_chromosome ]
            FilesParameters[ 'outputfile_like_bed' ] = options.output_file + '_' + each_chromosome + '.like_bed'
            FilesParameters[ 'outputfile_like_wig' ] = options.output_file + '_' + each_chromosome + '.like_wig'
            FilesParameters[ 'Chosen_Chromosome_Length' ] = InputData_Summary.chr2lengthSetting_dict[ each_chromosome ]
            Positioning_improvement( FilesParameters , ConvolutionParameters , threshold )
        parameters_Merge_results = { 'Input_Files': options.output_file,
                                     'Chromosome_list': InputData_Summary.chr_list,
                                     'OutputFile':  options.output_file + '_Gathering.like_bed',
                                     }
        Merge_overall_results( parameters_Merge_results )

print('\nNucDetect Finished,good luck\n\n')



