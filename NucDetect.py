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



    def collect_and_sort_shoulders(self):
        print('Collecting and sorting the shoulders with main nucleosomes in every peak region ......')
        peak_nucleosome_shoulder_simple_collection={}
        for one_shoulder in self.shoulder_list:
            inflection_pair_index=one_shoulder[4]    
            peak_index=one_shoulder[-1][-1]    
            if peak_index in self.nucleosome_dict.keys():
                if peak_index not in peak_nucleosome_shoulder_simple_collection.keys():
                    peak_nucleosome_shoulder_simple_collection[peak_index]=[]
                    peak_nucleosome_shoulder_simple_collection[peak_index].append(inflection_pair_index)
                elif peak_index in peak_nucleosome_shoulder_simple_collection.keys():
                    peak_nucleosome_shoulder_simple_collection[peak_index].append(inflection_pair_index)
            elif peak_index not in self.nucleosome_dict.keys():
                pass

        for peak_index in sorted( peak_nucleosome_shoulder_simple_collection.keys() ):
            self.shoulders_sorted_with_main_nucleosomes[peak_index]={}   
            main_nucleosome_in_the_peak=self.nucleosome_dict[peak_index]
            left_part=[]
            right_part=[]
            for inflection_pair_index in peak_nucleosome_shoulder_simple_collection[peak_index]:
                if self.inflection_pairs_dict_with_midpoint[inflection_pair_index][0]>self.nucleosome_dict[peak_index][1]:
                    right_part.append(inflection_pair_index)
                elif self.inflection_pairs_dict_with_midpoint[inflection_pair_index][1]<self.nucleosome_dict[peak_index][0]:
                    left_part.append(inflection_pair_index)
            if len(left_part)==1:
                self.shoulders_sorted_with_main_nucleosomes[peak_index][-1]=left_part[0]
            elif len(left_part)>1:
                for i in range(0,len(left_part)-1):
                    for j in range(i+1,len(left_part)):
                        if self.inflection_pairs_dict_with_midpoint[left_part[i]][1]<self.inflection_pairs_dict_with_midpoint[left_part[j]][0]:
                            temp=left_part[i]
                            left_part[i]=left_part[j]
                            left_part[j]=temp
                        else:
                            pass
                for k in range(len(left_part)):
                    self.shoulders_sorted_with_main_nucleosomes[peak_index][(k+1)*(-1)]=left_part[k]
            if len(right_part)==1:
                self.shoulders_sorted_with_main_nucleosomes[peak_index][1]=right_part[0]
            elif len(right_part)>1:
                for i in range(0,len(right_part)-1):
                    for j in range(i+1,len(right_part)):
                        if self.inflection_pairs_dict_with_midpoint[right_part[i]][0]>self.inflection_pairs_dict_with_midpoint[right_part[j]][1]:
                            temp=right_part[i]
                            right_part[i]=right_part[j]
                            right_part[j]=temp
                        else:
                            pass
                for k in range(len(right_part)):
                    self.shoulders_sorted_with_main_nucleosomes[peak_index][(k+1)]=right_part[k]
        


    def precise_nucleosome_position(self):
        def Pearson_Correlation_Coefficient(Original,Smoothed):
            meanO=sum(Original)/len(Original)
            meanS=sum(Smoothed)/len(Smoothed)
            sumOO=0
            sumSS=0
            sumOS=0
            for i in range(len(Original)):
                sumOO=sumOO+(Original[i]-meanO)*(Original[i]-meanO)
                sumSS=sumSS+(Smoothed[i]-meanS)*(Smoothed[i]-meanS)
                sumOS=sumOS+(Original[i]-meanO)*(Smoothed[i]-meanS)
            pcc_result=(sumOS+0.00000001)/math.sqrt((sumOO+0.00000001)*(sumSS+0.00000001))
            return pcc_result

        def concavity_and_convexity(Smoothed):
            BG=min(Smoothed)
            k=(Smoothed[-1]-Smoothed[0])/len(Smoothed)
            y_list=[]
            for i in range(len(Smoothed)):
                y=(i-0)*k+Smoothed[0]
                y_list.append(y)
            ratio=(sum([s-BG for s in Smoothed])+0.0000001)/(sum([t-BG for t in y_list])+0.0000001)
            return ratio

        def height_ratio(smallONE,bigONE):
            small_original_score_fragment=[]
            for nnp in range(smallONE[2],smallONE[3]+1):
                small_original_score_fragment.append(self.score_list[nnp])
            smallHEIGHT=max(small_original_score_fragment)
            big_original_score_fragment=[]
            for nnp in range(bigONE[2],bigONE[3]+1):
                big_original_score_fragment.append(self.score_list[nnp])
            bigHEIGHT=max(big_original_score_fragment)
            ratio=smallHEIGHT/bigHEIGHT
            return ratio

        def preparation(self,fragment_beginning_index,fragment_ending_index,gap_length,LorR,peak_index,LeftONE,RightONE):
            PCC_score=[]
            CaC_score=[]
            if LorR=='left':
                for midpoint_index in fragment_beginning_index:
                    Original=self.score_list[midpoint_index:fragment_ending_index+1]
                    Smoothed=self.Gaussian_list[midpoint_index:fragment_ending_index+1]
                    pcc_result=Pearson_Correlation_Coefficient(Original,Smoothed)
                    PCC_score.append(pcc_result)
                    cac_result=concavity_and_convexity(Smoothed)
                    CaC_score.append(cac_result)
                    height_ratio_score=height_ratio(LeftONE,RightONE)    
            elif LorR=='right':
                for midpoint_index in fragment_ending_index:
                    Original=self.score_list[fragment_beginning_index:midpoint_index+1]
                    Smoothed=self.Gaussian_list[fragment_beginning_index:midpoint_index+1]
                    pcc_result=Pearson_Correlation_Coefficient(Original,Smoothed)
                    PCC_score.append(pcc_result)
                    cac_result=concavity_and_convexity(Smoothed)
                    CaC_score.append(cac_result)
                    height_ratio_score=height_ratio(RightONE,LeftONE)  
            if min(PCC_score)>self.threshold['PCC_independence_1'] and max(CaC_score)>self.threshold['concavity_and_convexity_1'] and gap_length<self.threshold['gap_length_1'] and height_ratio_score>=self.threshold['height_ratio_1']:
                relationship='shifting'
            elif min(PCC_score)>self.threshold['PCC_independence_2'] and max(CaC_score)>self.threshold['concavity_and_convexity_2'] and gap_length<self.threshold['gap_length_2'] and height_ratio_score>=self.threshold['height_ratio_2']:
                relationship='shifting'
            elif min(PCC_score)>self.threshold['PCC_independence_3'] and max(CaC_score)>self.threshold['concavity_and_convexity_3'] and gap_length<self.threshold['gap_length_3'] and height_ratio_score>=self.threshold['height_ratio_3']:
                relationship='shifting'
            elif min(PCC_score)>self.threshold['PCC_independence_4'] and max(CaC_score)>self.threshold['concavity_and_convexity_4'] and gap_length<self.threshold['gap_length_4'] and height_ratio_score>=self.threshold['height_ratio_4']:
                relationship='shifting'
            elif min(PCC_score)>self.threshold['PCC_independence_5'] and max(CaC_score)>self.threshold['concavity_and_convexity_5'] and gap_length<self.threshold['gap_length_5'] and height_ratio_score>=self.threshold['height_ratio_5']:
                relationship='shifting'
            elif min(PCC_score)>self.threshold['PCC_independence_6'] and max(CaC_score)>self.threshold['concavity_and_convexity_6'] and gap_length<self.threshold['gap_length_6'] and height_ratio_score>=self.threshold['height_ratio_6']:
                relationship='shifting'
            else:
                relationship='independent'
            return relationship
            
        def relation_determination(self,peak_index,shoulders_in_a_peak):
            shoulders_in_a_peak_relation={}
            for position in sorted(shoulders_in_a_peak.keys()):
                inflection_pair_index=shoulders_in_a_peak[position]
                if len(self.inflection_pairs_dict_with_midpoint[inflection_pair_index])>5:
                    midpoints_list_for_pair=self.inflection_pairs_dict_with_midpoint[inflection_pair_index][5:]
                else:
                    midpoints_list_for_pair=[self.inflection_pairs_dict_with_midpoint[inflection_pair_index][:]]
                midpoints_index_list_for_pair=[p[2] for p in midpoints_list_for_pair]+[p[3] for p in midpoints_list_for_pair]
                if position<-1:
                    fragment_beginning_index = midpoints_index_list_for_pair
                    fragment_ending_index = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position+1]][2]-1
                    gap_length = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position+1]][2]-self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]
                    LeftONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]]
                    RightONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position+1]]
                    independent_or_shifting = preparation(self,fragment_beginning_index,fragment_ending_index,gap_length,'left',peak_index,LeftONE,RightONE)
                elif position==-1:
                    fragment_beginning_index = midpoints_index_list_for_pair
                    fragment_ending_index = self.nucleosome_dict[peak_index][2]-1
                    gap_length = self.nucleosome_dict[peak_index][2]-self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]
                    LeftONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]]
                    RightONE = self.nucleosome_dict[peak_index]
                    independent_or_shifting = preparation(self,fragment_beginning_index,fragment_ending_index,gap_length,'left',peak_index,LeftONE,RightONE)
                elif position==1:
                    fragment_beginning_index = self.nucleosome_dict[peak_index][3]+1
                    fragment_ending_index = midpoints_index_list_for_pair
                    gap_length = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]-self.nucleosome_dict[peak_index][3]
                    LeftONE = self.nucleosome_dict[peak_index]
                    RightONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]]
                    independent_or_shifting = preparation(self,fragment_beginning_index,fragment_ending_index,gap_length,'right',peak_index,LeftONE,RightONE)
                elif position>1:
                    fragment_beginning_index = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position-1]][3]+1
                    fragment_ending_index = midpoints_index_list_for_pair
                    gap_length = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]-self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position-1]][3]
                    LeftONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position-1]]
                    RightONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]]
                    independent_or_shifting = preparation(self,fragment_beginning_index,fragment_ending_index,gap_length,'right',peak_index,LeftONE,RightONE)
                shoulders_in_a_peak_relation[position]=independent_or_shifting
            return shoulders_in_a_peak_relation    

        self.more_independent_shoulders=[]
        shifting_counting=0
        independent_counting=0
        shoulder_candidates_counting=0
        main_nucleosome_integrated_with_shoulder_counting=0
        for peak_index in sorted(self.shoulders_sorted_with_main_nucleosomes.keys()):
            shoulders_in_a_peak=self.shoulders_sorted_with_main_nucleosomes[peak_index]
            shoulders_in_a_peak_relation=relation_determination(self,peak_index,shoulders_in_a_peak)
            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
            main_nucleosome_has_been_treated='no'  
            for position in sorted(shoulders_in_a_peak_relation.keys()):
                shoulder_candidates_counting=shoulder_candidates_counting+1
                if position<-1:
                    if shoulders_in_a_peak_relation[position]=='independent':
                        independent_counting=independent_counting+1    ### 记录independent
                        if one_independent_shoulder['left_index']=='none':
                            one_independent_shoulder['left_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]    
                        elif one_independent_shoulder['left_index']!='none':
                            pass
                        one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]   
                     
                        self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                        one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                    elif shoulders_in_a_peak_relation[position]=='shifting':
                        shifting_counting=shifting_counting+1    
                        if one_independent_shoulder['left_index']=='none':
                            one_independent_shoulder['left_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]
                        elif one_independent_shoulder['left_index']!='none':
                            pass
                elif position==-1:
                    if shoulders_in_a_peak_relation[position]=='independent':
                        independent_counting=independent_counting+1    
                        if one_independent_shoulder['left_index']=='none':
                            one_independent_shoulder['left_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]    
                        elif one_independent_shoulder['left_index']!='none':
                            pass
                        one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]   
                        self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                        one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}   
                    elif shoulders_in_a_peak_relation[position]=='shifting':
                        shifting_counting=shifting_counting+1    
                        if one_independent_shoulder['left_index']=='none':
                            one_independent_shoulder['left_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]
                        elif one_independent_shoulder['left_index']!='none':
                            pass
                        if max(shoulders_in_a_peak_relation.keys())==-1:  。
                            one_independent_shoulder['right_index']=self.nucleosome_dict[peak_index][3]    
                            self.nucleosome_dict[peak_index][2]=one_independent_shoulder['left_index']
                            self.nucleosome_dict[peak_index][3]=one_independent_shoulder['right_index']
                            self.nucleosome_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                            self.nucleosome_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                            main_nucleosome_integrated_with_shoulder_counting=main_nucleosome_integrated_with_shoulder_counting+1
                            self.nucleosome_dict[peak_index].append('Main_nucleosome:integrated')
                            main_nucleosome_has_been_treated='yes'    
                            one_independent_shoulders={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif max(shoulders_in_a_peak_relation.keys())>=1:
                            pass
                elif position==1:
                    if shoulders_in_a_peak_relation[position]=='independent':
                        independent_counting=independent_counting+1  
                        if one_independent_shoulder['left_index']!='none':   
                            one_independent_shoulder['right_index']=self.nucleosome_dict[peak_index][3]  
                            self.nucleosome_dict={peak_index:[pairBC, pairEC, pairBI, pairEI, pair_index, 'height', 'area', [peakBC, maxBC, maxEC, peakEC, peakBI, maxBI, maxEI, peakEI, peak_index]], ...}
                            self.nucleosome_dict[peak_index][2]=one_independent_shoulder['left_index']
                            self.nucleosome_dict[peak_index][3]=one_independent_shoulder['right_index']
                            self.nucleosome_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                            self.nucleosome_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                            main_nucleosome_integrated_with_shoulder_counting=main_nucleosome_integrated_with_shoulder_counting+1
                            self.nucleosome_dict[peak_index].append('Main_nucleosome:integrated')
                            main_nucleosome_has_been_treated='yes'    
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif one_independent_shoulder['left_index']=='none':   
                        
                            main_nucleosome_has_been_treated='yes'
                          
                        one_independent_shoulder['left_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]
                        if max(shoulders_in_a_peak_relation.keys())==position:    
                            one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]    
                            self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif max(shoulders_in_a_peak_relation.keys())>position:
                            pass
                    elif shoulders_in_a_peak_relation[position]=='shifting':
                        shifting_counting=shifting_counting+1   
                        if one_independent_shoulder['left_index']!='none':    
                            pass
                        elif one_independent_shoulder['left_index']=='none':    
                            one_independent_shoulder['left_index']=self.nucleosome_dict[peak_index][2]    
                        if max(shoulders_in_a_peak_relation.keys())==position:    
                            one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]    
                            self.nucleosome_dict[peak_index][2]=one_independent_shoulder['left_index']
                            self.nucleosome_dict[peak_index][3]=one_independent_shoulder['right_index']
                            self.nucleosome_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                            self.nucleosome_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                            main_nucleosome_integrated_with_shoulder_counting=main_nucleosome_integrated_with_shoulder_counting+1
                            self.nucleosome_dict[peak_index].append('Main_nucleosome:integrated')
                            main_nucleosome_has_been_treated='yes'    
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif max(shoulders_in_a_peak_relation.keys())>position:
                            pass
                elif position>1:
                    if shoulders_in_a_peak_relation[position]=='independent':
                        independent_counting=independent_counting+1    
                        if one_independent_shoulder['left_index']!='none':    
                            one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position-1]][3]    
                            if main_nucleosome_has_been_treated=='yes':    
                                self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                            elif main_nucleosome_has_been_treated=='no':    
                                self.nucleosome_dict[peak_index][2]=one_independent_shoulder['left_index']
                                self.nucleosome_dict[peak_index][3]=one_independent_shoulder['right_index']
                                self.nucleosome_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                                self.nucleosome_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                                main_nucleosome_integrated_with_shoulder_counting=main_nucleosome_integrated_with_shoulder_counting+1
                                self.nucleosome_dict[peak_index].append('Main_nucleosome:integrated')
                                main_nucleosome_has_been_treated='yes'    
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}   
                        elif one_independent_shoulder['left_index']=='none':
                            pass
                      
                        one_independent_shoulder['left_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][2]
                        if max(shoulders_in_a_peak_relation.keys())==position:    
                            one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]    
                            self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif max(shoulders_in_a_peak_relation.keys())>position:
                            pass
                    elif shoulders_in_a_peak_relation[position]=='shifting':
                        shifting_counting=shifting_counting+1  
                        if one_independent_shoulder['left_index']!='none':
                            pass
                        elif one_independent_shoulder['left_index']=='none':
                            pass
                        if max(shoulders_in_a_peak_relation.keys())==position:   
                            one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]    
                            if main_nucleosome_has_been_treated=='yes':   
                                self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                            elif main_nucleosome_has_been_treated=='no':   
                                self.nucleosome_dict[peak_index][2]=one_independent_shoulder['left_index']
                                self.nucleosome_dict[peak_index][3]=one_independent_shoulder['right_index']
                                self.nucleosome_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                                self.nucleosome_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                                main_nucleosome_integrated_with_shoulder_counting=main_nucleosome_integrated_with_shoulder_counting+1
                                self.nucleosome_dict[peak_index].append('Main_nucleosome:integrated')
                                main_nucleosome_has_been_treated='yes'   
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'} 
                        elif max(shoulders_in_a_peak_relation.keys())>position:
                            pass

       
        self.nucleosome_integrated=[]
        shoulder_counting=0
        shoulder_counting_suplimit=len(self.more_independent_shoulders)-1
        peak_index=1
        if len( self.nucleosome_dict.keys() ) > 0:
            peak_index_suplimit = sorted( list( self.nucleosome_dict.keys() ) )[-1]  
        else:
            peak_index_suplimit = 0
        total_shoulders=0
        total_main_nucleosomes=0
        isolated_main_nucleosome_counting=0
        if len( self.more_independent_shoulders ) > 0:
            for bp in range(len(self.score_table)):
                if self.more_independent_shoulders[shoulder_counting][0]==self.score_table[bp][0]:
                    self.nucleosome_integrated.append([self.more_independent_shoulders[shoulder_counting][0],self.more_independent_shoulders[shoulder_counting][1],self.more_independent_shoulders[shoulder_counting][2],self.more_independent_shoulders[shoulder_counting][3],'height','area',[self.peak_dict[peak_index][0],self.peak_dict[peak_index][1],self.peak_dict[peak_index][2],self.peak_dict[peak_index][3],self.peak_dict[peak_index][4],self.peak_dict[peak_index][5],self.peak_dict[peak_index][6],self.peak_dict[peak_index][7],peak_index]])
                    self.nucleosome_integrated[-1].append('Independent_shoulder')
                    total_shoulders=total_shoulders+1   
                    if shoulder_counting<shoulder_counting_suplimit:
                        shoulder_counting=shoulder_counting+1
                    else:
                        pass
                else:
                    while peak_index not in self.nucleosome_dict.keys():
                        peak_index=peak_index+1
                    if self.nucleosome_dict[peak_index][0]==self.score_table[bp][0]:
                        self.nucleosome_integrated.append(self.nucleosome_dict[peak_index])
                        self.nucleosome_integrated[-1].pop(4)
                        if self.nucleosome_integrated[-1][-1]=='Main_nucleosome:integrated':
                            pass
                        elif self.nucleosome_integrated[-1][-1]!='Main_nucleosome:integrated':
                            self.nucleosome_integrated[-1].append('Main_nucleosome:isolated')
                            isolated_main_nucleosome_counting=isolated_main_nucleosome_counting+1
                        total_main_nucleosomes=total_main_nucleosomes+1
                        if peak_index<peak_index_suplimit:
                            peak_index=peak_index+1
                        else:
                            pass
        else:
            for bp in range(len(self.score_table)):
                while peak_index not in self.nucleosome_dict.keys():
                    peak_index=peak_index+1
                if self.nucleosome_dict[peak_index][0]==self.score_table[bp][0]:
                    self.nucleosome_integrated.append(self.nucleosome_dict[peak_index])
                    self.nucleosome_integrated[-1].pop(4)
                    if self.nucleosome_integrated[-1][-1]=='Main_nucleosome:integrated':
                        pass
                    elif self.nucleosome_integrated[-1][-1]!='Main_nucleosome:integrated':
                        self.nucleosome_integrated[-1].append('Main_nucleosome:isolated')
                        isolated_main_nucleosome_counting=isolated_main_nucleosome_counting+1
                    total_main_nucleosomes=total_main_nucleosomes+1
                    if peak_index<peak_index_suplimit:
                        peak_index=peak_index+1
                    else:
                        pass
        print('\tChecking ==> ',total_main_nucleosomes,'+',total_shoulders,'=',total_main_nucleosomes+total_shoulders)
        self.final_nucleosome=[]
        for i in range(len(self.nucleosome_integrated)):    
            self.final_nucleosome.append(self.nucleosome_integrated[i][0:6]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][7]])
        print('...... Precise nucleosome positioning is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))


            
    def adjust_border(self):
        print('Adjust the border of every inflection pair with secondary inflection detection ......')
        def influence_coefficient(self,ii,influence_from_left_or_right):
            if influence_from_left_or_right=='left':
                influence_index=ii-1
                distance=self.nucleosome_integrated[ii][2]-self.nucleosome_integrated[influence_index][3]    
                if self.nucleosome_integrated[ii][-1]=='Main_nucleosome:integrated' or self.nucleosome_integrated[ii][-1]=='Main_nucleosome:isolated':
                    center_1=self.nucleosome_integrated[ii][6][5]
                elif self.nucleosome_integrated[ii][-1]=='Independent_shoulder':
                    center_1=(self.nucleosome_integrated[ii][2]+self.nucleosome_integrated[ii][3])/2
                if self.nucleosome_integrated[influence_index][-1]=='Main_nucleosome:integrated' or self.nucleosome_integrated[influence_index][-1]=='Main_nucleosome:isolated':
                    center_2=self.nucleosome_integrated[influence_index][6][6]
                elif self.nucleosome_integrated[influence_index][-1]=='Independent_shoulder':
                    center_2=(self.nucleosome_integrated[influence_index][2]+self.nucleosome_integrated[influence_index][3])/2
                distance_center=center_1-center_2    
            elif influence_from_left_or_right=='right':
                influence_index=ii+1
                distance=self.nucleosome_integrated[influence_index][2]-self.nucleosome_integrated[ii][3]   
                if self.nucleosome_integrated[influence_index][-1]=='Main_nucleosome:integrated' or self.nucleosome_integrated[influence_index][-1]=='Main_nucleosome:isolated':
                    center_1=self.nucleosome_integrated[influence_index][6][5]
                elif self.nucleosome_integrated[influence_index][-1]=='Independent_shoulder':
                    center_1=(self.nucleosome_integrated[influence_index][2]+self.nucleosome_integrated[influence_index][3])/2
                if self.nucleosome_integrated[ii][-1]=='Main_nucleosome:integrated' or self.nucleosome_integrated[ii][-1]=='Main_nucleosome:isolated':
                    center_2=self.nucleosome_integrated[ii][6][6]
                elif self.nucleosome_integrated[ii][-1]=='Independent_shoulder':
                    center_2=(self.nucleosome_integrated[ii][2]+self.nucleosome_integrated[ii][3])/2
                distance_center=center_1-center_2    
            temp_height_self=[]
            temp_smoothedheight_self=[]
            for nnp in range(self.nucleosome_integrated[ii][2],self.nucleosome_integrated[ii][3]+1):
                temp_height_self.append(self.score_table[nnp][1])
                temp_smoothedheight_self.append(self.score_table[nnp][2])
            height_self=max(temp_height_self)
            area_self=sum(temp_height_self)
            smoothedheight_self=max(temp_smoothedheight_self)
            smoothedarea_self=sum(temp_smoothedheight_self)
            temp_height_influence=[]
            temp_smoothedheight_influence=[]
            for nnp in range(self.nucleosome_integrated[influence_index][2],self.nucleosome_integrated[influence_index][3]+1):
                temp_height_influence.append(self.score_table[nnp][1])
                temp_smoothedheight_influence.append(self.score_table[nnp][2])
            height_influence=max(temp_height_influence)
            area_influence=sum(temp_height_influence)
            smoothedheight_influence=max(temp_smoothedheight_influence)
            smoothedarea_influence=sum(temp_smoothedheight_influence)
            coefficient=(height_influence/height_self)*10/distance
            coefficient_absolute=height_influence*10/distance
            if height_influence<height_self or smoothedheight_influence<smoothedheight_self:
                effect='no'
            elif distance>self.threshold['influence_coefficient_distance'] or distance_center>self.threshold['influence_coefficient_distance_center']:
                effect='no'
            else:
                if coefficient<self.threshold['influence_coefficient_cutoff'] or coefficient_absolute<self.threshold['influence_coefficient_absolute_cutoff']:
                    effect='no'
                else:
                    effect='yes'
            return [effect,coefficient]

        influence_list=[]
        for i in range(len(self.nucleosome_integrated)):
            influence_list.append([])
            peak_index=self.nucleosome_integrated[i][-2]
            if i==0:
                right_effect=influence_coefficient(self,i,'right')
                if right_effect[0]=='yes':
                    influence_list[-1].append('right')
                else:
                    pass
            elif i>0 and i<len(self.nucleosome_integrated)-1:
                right_effect=influence_coefficient(self,i,'right')
                left_effect=influence_coefficient(self,i,'left')
                if right_effect[0]=='yes':
                    influence_list[-1].append('right')
                else:
                    pass
                if left_effect[0]=='yes':
                    influence_list[-1].append('left')
                else:
                    pass
            elif i==len(self.nucleosome_integrated)-1:
                left_effect=influence_coefficient(self,i,'left')
                if left_effect[0]=='yes':
                    influence_list[-1].append('left')
                else:
                    pass

        if len( self.nucleosome_dict.keys() ) > 0:
            peak_index_suplimit = sorted( list(self.nucleosome_dict.keys()) )[-1] 
        else:
            peak_index_suplimit = 0
        for j in range(len(self.nucleosome_integrated)):
            peak_index=self.nucleosome_integrated[j][-2][-1]
            adjacent_influence=influence_list[j]
            secondary_canditates=[]
            if peak_index==1:
                peak_index_list=[peak_index,peak_index+1]
            elif peak_index>1 and peak_index<peak_index_suplimit:
                peak_index_list=[peak_index-1,peak_index,peak_index+1]
            elif peak_index==peak_index_suplimit:
                peak_index_list=[peak_index-1,peak_index]
            for peak_index_number in peak_index_list:
                if peak_index_number in self.secondary_inflection_pairs_dict_with_peak.keys():
                    for candidate in self.secondary_inflection_pairs_dict_with_peak[peak_index_number]:
                        if candidate not in secondary_canditates:
                            secondary_canditates.append(candidate)
                        else:
                            pass
                else:
                    pass
            temp_left_index=0
            temp_left=secondary_canditates[0][2]-self.nucleosome_integrated[j][2]
            temp_right_index=len(secondary_canditates)-1
            temp_right=secondary_canditates[len(secondary_canditates)-1][3]-self.nucleosome_integrated[j][3]
            for k in range(len(secondary_canditates)):
                if abs(secondary_canditates[k][2]-self.nucleosome_integrated[j][2])<abs(temp_left):
                    temp_left_index=k
                    temp_left=secondary_canditates[k][2]-self.nucleosome_integrated[j][2]
                    if 'left' in adjacent_influence:
                        if secondary_canditates[k][2]<self.nucleosome_integrated[j][2] and secondary_canditates[k][2]>(self.nucleosome_integrated[j-1][3]+self.nucleosome_integrated[j][2])/2:
                            self.nucleosome_integrated[j][0]=secondary_canditates[temp_left_index][0]
                            self.nucleosome_integrated[j][2]=secondary_canditates[temp_left_index][2]
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
                if abs(secondary_canditates[len(secondary_canditates)-1-k][3]-self.nucleosome_integrated[j][3])<abs(temp_right):
                    temp_right_index=len(secondary_canditates)-1-k
                    temp_right=secondary_canditates[len(secondary_canditates)-1-k][3]-self.nucleosome_integrated[j][3]
                    if 'right' in adjacent_influence:
                        if secondary_canditates[len(secondary_canditates)-1-k][3]>self.nucleosome_integrated[j][3] and secondary_canditates[len(secondary_canditates)-1-k][3]<(self.nucleosome_integrated[j][3]+self.nucleosome_integrated[j+1][2])/2:
                            self.nucleosome_integrated[j][1]=secondary_canditates[temp_right_index][1]
                            self.nucleosome_integrated[j][3]=secondary_canditates[temp_right_index][3]
                        else:
                            pass
                    else:
                        pass
                else:
                    pass

        for i in range(len(self.nucleosome_integrated)):    
            self.final_nucleosome.append(self.nucleosome_integrated[i][0:6]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][7]])
        print('...... Border adjustment is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))



    def filter_results(self):            

        print('Filter and annotate the nucleosomes ......')
        inteStream_Downal_property_list=[]
        for n in range(len(self.nucleosome_integrated)-1):
            one_inteStream_Downal={'Distance_between_center':'none',  'InteStream_Downal_length':'none',  'Original_InteStream_Downal_Depth':'none',  'Smoothed_InteStream_Downal_Depth':'none',
                          'Left_OriginalHeight':'none',   'Left_OriginalAUC':'none',   'Left_SmoothedHeight':'none',   'Left_SmoothedAUC':'none',
                          'Right_OriginalHeight':'none',  'Right_OriginalAUC':'none',  'Right_SmoothedHeight':'none',  'Right_SmoothedAUC':'none',
                          'Smoothed_InteStream_Downal_Depth_percentage':'none',
                          'HeightRatio':'none',
                          'Merge':'no',
                          'Noise_on_left_or_right':'none'}
            one_inteStream_Downal['Distance_between_center']=(self.nucleosome_integrated[n+1][0]+self.nucleosome_integrated[n+1][1])/2-(self.nucleosome_integrated[n][0]+self.nucleosome_integrated[n][1])/2
            one_inteStream_Downal['InteStream_Downal_length']=self.nucleosome_integrated[n+1][0]-self.nucleosome_integrated[n][1]  
            InteStream_Downal_original_score_fragment=[]
            InteStream_Downal_smoothed_score_fragment=[]
            for nnq in range(self.nucleosome_integrated[n][3]+1,self.nucleosome_integrated[n+1][2]):
                InteStream_Downal_original_score_fragment.append(self.score_list[nnq])
                InteStream_Downal_smoothed_score_fragment.append(self.Gaussian_list[nnq])
            one_inteStream_Downal['Original_InteStream_Downal_Depth']=min(InteStream_Downal_original_score_fragment)
            one_inteStream_Downal['Smoothed_InteStream_Downal_Depth']=min(InteStream_Downal_smoothed_score_fragment)

            Left_original_score_fragment=[]
            Left_smoothed_score_fragment=[]
            for nnp in range(self.nucleosome_integrated[n][2],self.nucleosome_integrated[n][3]+1):
                Left_original_score_fragment.append(self.score_list[nnp])
                Left_smoothed_score_fragment.append(self.Gaussian_list[nnp])
            one_inteStream_Downal['Left_OriginalHeight']=max(Left_original_score_fragment)
            one_inteStream_Downal['Left_OriginalAUC']=sum(Left_original_score_fragment)
            one_inteStream_Downal['Left_SmoothedHeight']=max(Left_smoothed_score_fragment)
            one_inteStream_Downal['Left_SmoothedAUC']=sum(Left_smoothed_score_fragment)

            Right_original_score_fragment=[]
            Right_smoothed_score_fragment=[]
            for nnp in range(self.nucleosome_integrated[n+1][2],self.nucleosome_integrated[n+1][3]+1):
                Right_original_score_fragment.append(self.score_list[nnp])
                Right_smoothed_score_fragment.append(self.Gaussian_list[nnp])
            one_inteStream_Downal['Right_OriginalHeight']=max(Right_original_score_fragment)
            one_inteStream_Downal['Right_OriginalAUC']=sum(Right_original_score_fragment)
            one_inteStream_Downal['Right_SmoothedHeight']=max(Right_smoothed_score_fragment)
            one_inteStream_Downal['Right_SmoothedAUC']=sum(Right_smoothed_score_fragment)
            one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']=2*one_inteStream_Downal['Smoothed_InteStream_Downal_Depth']/(one_inteStream_Downal['Left_SmoothedHeight']+one_inteStream_Downal['Right_SmoothedHeight'])    ### 关键
            one_inteStream_Downal['HeightRatio']=min(one_inteStream_Downal['Left_OriginalHeight'],one_inteStream_Downal['Right_OriginalHeight'])/max(one_inteStream_Downal['Left_OriginalHeight'],one_inteStream_Downal['Right_OriginalHeight'])    ### 关键
 
            if one_inteStream_Downal['Distance_between_center']<self.threshold['merging_center_distance']:
                if max(one_inteStream_Downal['Left_OriginalHeight'],one_inteStream_Downal['Right_OriginalHeight'])<self.threshold['merging_height_watershed']:
                    if one_inteStream_Downal['HeightRatio']>self.threshold['merging_height_ratio_1'] and one_inteStream_Downal['InteStream_Downal_length']<=self.threshold['merging_gap_1'] and one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']>self.threshold['merging_percentage_1']:
                        one_inteStream_Downal['Merge']='yes'
                    elif one_inteStream_Downal['HeightRatio']>self.threshold['merging_height_ratio_3'] and one_inteStream_Downal['InteStream_Downal_length']<=self.threshold['merging_gap_3'] and one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']>self.threshold['merging_percentage_3']:
                        one_inteStream_Downal['Merge']='yes'
                    elif one_inteStream_Downal['HeightRatio']>self.threshold['merging_height_ratio_4'] and one_inteStream_Downal['InteStream_Downal_length']<=self.threshold['merging_gap_4'] and one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']>self.threshold['merging_percentage_4']:
                        one_inteStream_Downal['Merge']='yes'
                    else:
                        one_inteStream_Downal['Merge']='no'
                elif max(one_inteStream_Downal['Left_OriginalHeight'],one_inteStream_Downal['Right_OriginalHeight'])>=self.threshold['merging_height_watershed']:
                    if one_inteStream_Downal['HeightRatio']>self.threshold['merging_height_ratio_2'] and one_inteStream_Downal['InteStream_Downal_length']<=self.threshold['merging_gap_2'] and one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']>self.threshold['merging_percentage_2']:
                        one_inteStream_Downal['Merge']='yes'
                    elif one_inteStream_Downal['HeightRatio']>self.threshold['merging_height_ratio_3'] and one_inteStream_Downal['InteStream_Downal_length']<=self.threshold['merging_gap_3'] and one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']>self.threshold['merging_percentage_3']:
                        one_inteStream_Downal['Merge']='yes'
                    elif one_inteStream_Downal['HeightRatio']>self.threshold['merging_height_ratio_4'] and one_inteStream_Downal['InteStream_Downal_length']<=self.threshold['merging_gap_4'] and one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']>self.threshold['merging_percentage_4']:
                        one_inteStream_Downal['Merge']='yes'
                    else:
                        one_inteStream_Downal['Merge']='no'
            else:
                one_inteStream_Downal['Merge']='no'
            inteStream_Downal_property_list.append(one_inteStream_Downal)

        print('\tChecking ==> Number of inteStream_Downal:  ',len(inteStream_Downal_property_list))
        for n in range(len(self.nucleosome_integrated)-2):
            if inteStream_Downal_property_list[n]['Right_OriginalHeight']==inteStream_Downal_property_list[n+1]['Left_OriginalHeight'] and inteStream_Downal_property_list[n]['Right_OriginalAUC']==inteStream_Downal_property_list[n+1]['Left_OriginalAUC'] and inteStream_Downal_property_list[n]['Right_SmoothedHeight']==inteStream_Downal_property_list[n+1]['Left_SmoothedHeight'] and inteStream_Downal_property_list[n]['Right_SmoothedAUC']==inteStream_Downal_property_list[n+1]['Left_SmoothedAUC']:
                pass
            else:
                print('ERROR in inteStream_Downal mapping:  ',n,' and ',n+1)


        def merge_and_annotate(self,inteStream_Downal_property_list):
            merging_choice_list=['none']*len(self.nucleosome_integrated)
            for i in range(len(self.nucleosome_integrated)):
                if i==0:
                    if inteStream_Downal_property_list[i]['Merge']=='yes':
                        merging_choice_list[i]='right'
                    else:
                        pass
                elif i>0 and i<len(self.nucleosome_integrated)-1:
                    if inteStream_Downal_property_list[i]['Merge']=='yes' and inteStream_Downal_property_list[i-1]['Merge']=='no':
                        merging_choice_list[i]='right'
                    elif inteStream_Downal_property_list[i]['Merge']=='no' and inteStream_Downal_property_list[i-1]['Merge']=='yes':
                        merging_choice_list[i]='left'
                    elif inteStream_Downal_property_list[i]['Merge']=='yes' and inteStream_Downal_property_list[i-1]['Merge']=='yes':
                        if inteStream_Downal_property_list[i]['Smoothed_InteStream_Downal_Depth']>inteStream_Downal_property_list[i-1]['Smoothed_InteStream_Downal_Depth']:
                            merging_choice_list[i]='right'
                        elif inteStream_Downal_property_list[i]['Smoothed_InteStream_Downal_Depth']<inteStream_Downal_property_list[i-1]['Smoothed_InteStream_Downal_Depth']:
                            merging_choice_list[i]='left'
                        elif inteStream_Downal_property_list[i]['Smoothed_InteStream_Downal_Depth']==inteStream_Downal_property_list[i-1]['Smoothed_InteStream_Downal_Depth']:
                                pass
                    elif inteStream_Downal_property_list[i]['Merge']=='no' and inteStream_Downal_property_list[i-1]['Merge']=='no':
                        pass
                elif i==len(self.nucleosome_integrated)-1:
                    if inteStream_Downal_property_list[i-1]['Merge']=='yes':
                        merging_choice_list[i]='left'
                    else:
                        pass
 
            merging_counting=0
            unmerging_counting=0
            merging_switch='off'
            final_nucleosome_temp_recorder=[]    
            one_nucleosome_for_final_nucleosome_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
            ### self.nucleosome_integrated=[[pairBC, pairEC, pairBI, pairEI, 'height', 'area', [peakBC, maxBC, maxEC, peakEC, peakBI, maxBI, maxEI, peakEI, peak_index], main/shoulder], ...]
            for i in range(len(self.nucleosome_integrated)):
                if merging_choice_list[i]=='right':
                    if merging_switch=='off':    
                        one_nucleosome_for_final_nucleosome_temp_recorder[0]=self.nucleosome_integrated[i][0]
                        one_nucleosome_for_final_nucleosome_temp_recorder[2]=self.nucleosome_integrated[i][2]
                        merging_switch='on'
                    elif merging_switch=='on':   
                        final_nucleosome_temp_recorder.append(self.nucleosome_integrated[i-1][0:6]+[self.nucleosome_integrated[i-1][6]]+[self.nucleosome_integrated[i-1][6]]+[self.nucleosome_integrated[i-1][7]])    ### 记录前一个
                        unmerging_counting=unmerging_counting+1
                        one_nucleosome_for_final_nucleosome_temp_recorder[0]=self.nucleosome_integrated[i][0]
                        one_nucleosome_for_final_nucleosome_temp_recorder[2]=self.nucleosome_integrated[i][2]

                elif merging_choice_list[i]=='none':
                    if merging_switch=='off':
                        final_nucleosome_temp_recorder.append(self.nucleosome_integrated[i][0:6]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][7]])    ### 记录这个
                        unmerging_counting=unmerging_counting+1
                    elif merging_switch=='on':
                        final_nucleosome_temp_recorder.append(self.nucleosome_integrated[i-1][0:6]+[self.nucleosome_integrated[i-1][6]]+[self.nucleosome_integrated[i-1][6]]+[self.nucleosome_integrated[i-1][7]])    ### 记录前一个
                        final_nucleosome_temp_recorder.append(self.nucleosome_integrated[i][0:6]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][7]])    ### 记录这个
                        unmerging_counting=unmerging_counting+2
                        merging_switch='off'
                        one_nucleosome_for_final_nucleosome_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
                elif merging_choice_list[i]=='left':
                    if merging_switch=='off':
                        final_nucleosome_temp_recorder.append(self.nucleosome_integrated[i][0:6]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][7]])    ### 记录这个
                        unmerging_counting=unmerging_counting+1
                    elif merging_switch=='on':
                        one_nucleosome_for_final_nucleosome_temp_recorder[1]=self.nucleosome_integrated[i][1]
                        one_nucleosome_for_final_nucleosome_temp_recorder[3]=self.nucleosome_integrated[i][3]
                        one_nucleosome_for_final_nucleosome_temp_recorder[7]=self.nucleosome_integrated[i][6]
                        if type(one_nucleosome_for_final_nucleosome_temp_recorder[0])==int and type(one_nucleosome_for_final_nucleosome_temp_recorder[1])==int and type(one_nucleosome_for_final_nucleosome_temp_recorder[2])==int and type(one_nucleosome_for_final_nucleosome_temp_recorder[3])==int and len(one_nucleosome_for_final_nucleosome_temp_recorder[6])==9 and len(one_nucleosome_for_final_nucleosome_temp_recorder[7])==9:
                            if ('nucleosome' in one_nucleosome_for_final_nucleosome_temp_recorder[8]) or ('nucleosome' in self.nucleosome_integrated[i][7]):
                                one_nucleosome_for_final_nucleosome_temp_recorder[8]='Main_nucleosome:integrated:doublet'
                            else:
                                one_nucleosome_for_final_nucleosome_temp_recorder[8]='Shoulder:doublet'
                            final_nucleosome_temp_recorder.append(one_nucleosome_for_final_nucleosome_temp_recorder)   
                            merging_counting=merging_counting+1
                            merging_switch='off'
                            one_nucleosome_for_final_nucleosome_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
                        else:
                            merging_switch='off'
                            one_nucleosome_for_final_nucleosome_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
                            print('Error in merging!')
            print('\tChecking ==> ',merging_counting,'+',unmerging_counting,'=',merging_counting+unmerging_counting)
            final_nucleosome_temp_recorder_noise_discarded=[]    
            discarded_noise=0
            for i in range(len(final_nucleosome_temp_recorder)):
                original_score_fragment=[]
                LoG_fragment=[]
                negative_counting=0
                negative_counting_switch='off'
                negative_counting_list=[]
                LoG_Sigma3_length=0
                for nnp in range(final_nucleosome_temp_recorder[i][2],final_nucleosome_temp_recorder[i][3]+1):
                    original_score_fragment.append(self.score_table[nnp][1])
                    LoG_fragment.append(self.score_table[nnp][3])
                    if self.score_table[nnp][3]<0:
                        if nnp<final_nucleosome_temp_recorder[i][3]:
                            negative_counting_switch='on'
                            negative_counting=negative_counting+1
                        elif nnp==final_nucleosome_temp_recorder[i][3]:
                            negative_counting=negative_counting+1
                            negative_counting_list.append(negative_counting)
                            negative_counting=0
                            negative_counting_switch='off'
                    elif self.score_table[nnp][3]>=0:
                        if negative_counting_switch=='on':
                            negative_counting_list.append(negative_counting)
                            negative_counting=0
                            negative_counting_switch='off'
                        elif negative_counting_switch=='off':
                            pass
                LoG_Sigma3_length=max(negative_counting_list)
                if sum(original_score_fragment)<self.threshold['discarded_noise_selfAUC_1'] and LoG_Sigma3_length<=self.threshold['discarded_noise_selfLength_1'] and sum(LoG_fragment)/len(LoG_fragment)>self.threshold['discarded_noise_LoG_average_1']:
                    discarded_noise=discarded_noise+1
                elif len(original_score_fragment)<=self.threshold['discarded_noise_selfLength_2'] and sum(LoG_fragment)/len(LoG_fragment)>self.threshold['discarded_noise_LoG_average_2']:
                    discarded_noise=discarded_noise+1
                elif max(original_score_fragment)<=self.threshold['discarded_noise_selfHeight_3'] and sum(LoG_fragment)/len(LoG_fragment)>self.threshold['discarded_noise_LoG_average_3']:
                    discarded_noise=discarded_noise+1
                elif LoG_Sigma3_length<=self.threshold['discarded_noise_RealSelfLength_4'] and sum(LoG_fragment)/len(LoG_fragment)>self.threshold['discarded_noise_LoG_average_4']:
                    discarded_noise=discarded_noise+1
                elif LoG_Sigma3_length<=self.threshold['discarded_noise_RealSelfLength_5'] and sum(LoG_fragment)/len(LoG_fragment)>self.threshold['discarded_noise_LoG_average_5']:
                    discarded_noise=discarded_noise+1
                elif LoG_Sigma3_length<=self.threshold['discarded_noise_RealSelfLength_6'] and sum(LoG_fragment)/len(LoG_fragment)>self.threshold['discarded_noise_LoG_average_6']:
                    discarded_noise=discarded_noise+1
                else:
                    final_nucleosome_temp_recorder_noise_discarded.append(final_nucleosome_temp_recorder[i])
            print('\tChecking ==> Discarded nucleosomes as noise:  ',discarded_noise)
            print('\tChecking ==> Final remaining nucleosomes:  ',len(final_nucleosome_temp_recorder_noise_discarded))
            return final_nucleosome_temp_recorder_noise_discarded

        if threshold['filter_switch']=='off':
            self.final_nucleosome=[]
            for i in range(len(self.nucleosome_integrated)):    
                self.final_nucleosome.append(self.nucleosome_integrated[i][0:6]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][6]]+[self.nucleosome_integrated[i][7]])
        elif threshold['filter_switch']=='on':
            self.final_nucleosome=[]
            self.final_nucleosome=merge_and_annotate(self,inteStream_Downal_property_list)
        print('\tChecking ==> Nucleosomes after filtered:  ',len(self.final_nucleosome))
        print('...... Results are filtered and annotated.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
            
    def Test_significance( self ):
        self.Log_Pvalue_Peak_dict = { }
        self.Log_Pvalue_Valley_dict = { }
        ########
        print('Test significance ......')
        next_done = 'no'
        beginning_coordinate = 1
        ending_coordinate = len( self.tag_list ) * 10 - 9
        for i in range( len( self.final_nucleosome ) ):
            sys.stdout.write( '\r    Nucleosome: '+str(i+1) )
            ### nucleosome=[[pairBC, pairEC, pairBI, pairEI, 'height', 'area', [peakBC,maxBC,maxEC,peakEC,peakBI,maxBI,maxEI,peakEI,BeginningPeakIndex], [peakBC,maxBC,maxEC,peakEC,peakBI,maxBI,maxEI,peakEI,EndingPeakIndex], main/shoulder], ...]
            LeftInflection  = self.final_nucleosome[i][0]    #### 峰的起点坐标
            RightInflection = self.final_nucleosome[i][1]    #### 峰的终点坐标
            Nucleosome_mid = ( LeftInflection + RightInflection ) // 2 // 10 * 10 + 1
            Nuc_BG_left_1000  = max( beginning_coordinate  , Nucleosome_mid - 1000 )
            Nuc_BG_right_1000 = min( Nucleosome_mid + 1000 , ending_coordinate     )
            Nuc_BG_left_5000  = max( beginning_coordinate  , Nucleosome_mid - 5000 )
            Nuc_BG_right_5000 = min( Nucleosome_mid + 5000 , ending_coordinate     )
            Nuc_BG_left_10000  = max( beginning_coordinate   , Nucleosome_mid - 10000 )
            Nuc_BG_right_10000 = min( Nucleosome_mid + 10000 , ending_coordinate      )
            ####
            LeftInflection_index  = ( LeftInflection - beginning_coordinate ) // 10
            RightInflection_index = ( RightInflection - beginning_coordinate ) // 10
            Nuc_BG_left_1000_index  = ( Nuc_BG_left_1000 - beginning_coordinate ) // 10
            Nuc_BG_right_1000_index = ( Nuc_BG_right_1000 - beginning_coordinate ) // 10
            Nuc_BG_left_5000_index  = ( Nuc_BG_left_5000 - beginning_coordinate ) // 10
            Nuc_BG_right_5000_index = ( Nuc_BG_right_5000 - beginning_coordinate ) // 10
            Nuc_BG_left_10000_index  = ( Nuc_BG_left_10000 - beginning_coordinate ) // 10
            Nuc_BG_right_10000_index = ( Nuc_BG_right_10000 - beginning_coordinate ) // 10
            BG_1000_for_peak  = sum( [ self.tag_list[index] for index in range( Nuc_BG_left_1000_index , Nuc_BG_right_1000_index+1 ) ] ) / (Nuc_BG_right_1000_index-Nuc_BG_left_1000_index+1) * (RightInflection_index-LeftInflection_index+1)
            BG_5000_for_peak  = sum( [ self.tag_list[index] for index in range( Nuc_BG_left_5000_index , Nuc_BG_right_5000_index+1 ) ] ) / (Nuc_BG_right_5000_index-Nuc_BG_left_5000_index+1) * (RightInflection_index-LeftInflection_index+1)
            BG_10000_for_peak = sum( [ self.tag_list[index] for index in range( Nuc_BG_left_10000_index , Nuc_BG_right_10000_index+1 ) ] ) / (Nuc_BG_right_10000_index-Nuc_BG_left_10000_index+1) * (RightInflection_index-LeftInflection_index+1)
            BG_for_peak = max( BG_1000_for_peak , BG_5000_for_peak , BG_10000_for_peak )
            foreground_for_peak = sum( [ self.tag_list[index] for index in range( LeftInflection_index , RightInflection_index+1 ) ] )
            Pvalue_for_peak , Score_for_peak = Poisson_test.greater_fast( BG_for_peak , foreground_for_peak )

            if i == 0:
                Upstream_end = self.final_nucleosome[i][0] - 10    
                Upstream_start = max( Upstream_end - 300 , beginning_coordinate )   
                Upstream_mid = ( Upstream_start + Upstream_end ) // 2 // 10 * 10 + 1
                Stream_Up_left_1000  = max( beginning_coordinate  , Upstream_mid - 1000 )
                Stream_Up_right_1000 = min( Upstream_mid + 1000 , ending_coordinate     )
                Stream_Up_left_5000  = max( beginning_coordinate  , Upstream_mid - 5000 )
                Stream_Up_right_5000 = min( Upstream_mid + 5000 , ending_coordinate     )
                Stream_Up_left_10000  = max( beginning_coordinate  , Upstream_mid - 10000 )
                Stream_Up_right_10000 = min( Upstream_mid + 10000 , ending_coordinate     )
                ####
                Upstream_start_index = ( Upstream_start - beginning_coordinate ) // 10
                Upstream_end_index   = ( Upstream_end - beginning_coordinate ) // 10
                Stream_Up_left_1000_index  = ( Stream_Up_left_1000 - beginning_coordinate ) // 10
                Stream_Up_right_1000_index = ( Stream_Up_right_1000 - beginning_coordinate ) // 10
                Stream_Up_left_5000_index  = ( Stream_Up_left_5000 - beginning_coordinate ) // 10
                Stream_Up_right_5000_index = ( Stream_Up_right_5000 - beginning_coordinate ) // 10
                Stream_Up_left_10000_index  = ( Stream_Up_left_10000 - beginning_coordinate ) // 10
                Stream_Up_right_10000_index = ( Stream_Up_right_10000 - beginning_coordinate ) // 10
                BG_1000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_1000_index , Stream_Up_right_1000_index+1 ) ] ) / (Stream_Up_right_1000_index-Stream_Up_left_1000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                BG_5000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_5000_index , Stream_Up_right_5000_index+1 ) ] ) / (Stream_Up_right_5000_index-Stream_Up_left_5000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                BG_10000_for_Up = sum( [ self.tag_list[index] for index in range( Stream_Up_left_10000_index , Stream_Up_right_10000_index+1 ) ] ) / (Stream_Up_right_10000_index-Stream_Up_left_10000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                BG_for_Up = min( BG_1000_for_Up , BG_5000_for_Up , BG_10000_for_Up )
                if BG_for_Up < 0:
                    BG_for_Up = 0
                foreground_for_Up = sum( [ self.tag_list[index] for index in range( Upstream_start_index , Upstream_end_index+1 ) ] )
                Pvalue_for_Up , Score_for_Up = Poisson_test.less_fast( BG_for_Up , foreground_for_Up )
                ########
                Stream_Down_start = self.final_nucleosome[i][1] + 10   
                if len( self.final_nucleosome ) == 1:
                    Stream_Down_end = min( Stream_Down_start + 300 , ending_coordinate )    
                elif len( self.final_nucleosome ) > 1:
                    if self.final_nucleosome[i+1][0] - Stream_Down_start <= 1000:
                        Stream_Down_end = self.final_nucleosome[i+1][0] - 10    
                        next_done = 'yes'
                    else:
                        Stream_Down_end = Stream_Down_start + 1000    
                        next_done = 'no'
                Stream_Down_mid = ( Stream_Down_start + Stream_Down_end ) // 2 // 10 * 10 + 1
                Stream_Down_left_1000  = max( beginning_coordinate   , Stream_Down_mid - 1000 )
                Stream_Down_right_1000 = min( Stream_Down_mid + 1000 , ending_coordinate      )
                Stream_Down_left_5000  = max( beginning_coordinate   , Stream_Down_mid - 5000 )
                Stream_Down_right_5000 = min( Stream_Down_mid + 5000 , ending_coordinate      )
                Stream_Down_left_10000  = max( beginning_coordinate   , Stream_Down_mid - 10000 )
                Stream_Down_right_10000 = min( Stream_Down_mid + 10000 , ending_coordinate      )
                Stream_Down_start_index = ( Stream_Down_start - beginning_coordinate ) // 10
                Stream_Down_end_index   = ( Stream_Down_end - beginning_coordinate ) // 10
                Stream_Down_left_1000_index  = ( Stream_Down_left_1000 - beginning_coordinate ) // 10
                Stream_Down_right_1000_index = ( Stream_Down_right_1000 - beginning_coordinate ) // 10
                Stream_Down_left_5000_index  = ( Stream_Down_left_5000 - beginning_coordinate ) // 10
                Stream_Down_right_5000_index = ( Stream_Down_right_5000 - beginning_coordinate ) // 10
                Stream_Down_left_10000_index  = ( Stream_Down_left_10000 - beginning_coordinate ) // 10
                Stream_Down_right_10000_index = ( Stream_Down_right_10000 - beginning_coordinate ) // 10
                BG_1000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_1000_index , Stream_Down_right_1000_index+1 ) ] ) / (Stream_Down_right_1000_index-Stream_Down_left_1000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_5000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_5000_index , Stream_Down_right_5000_index+1 ) ] ) / (Stream_Down_right_5000_index-Stream_Down_left_5000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_10000_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_left_10000_index , Stream_Down_right_10000_index+1 ) ] ) / (Stream_Down_right_10000_index-Stream_Down_left_10000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_for_Stream_Down = min( BG_1000_for_Stream_Down , BG_5000_for_Stream_Down , BG_10000_for_Stream_Down )
                if BG_for_Stream_Down < 0:
                    BG_for_Stream_Down = 0
                foreground_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_start_index , Stream_Down_end_index+1 ) ] )
                Pvalue_for_Stream_Down , Score_for_Stream_Down = Poisson_test.less_fast( BG_for_Stream_Down , foreground_for_Stream_Down )
            elif ( i > 0 ) and ( i < len( self.final_nucleosome ) - 1 ):
                if next_done == 'yes':
                    Pvalue_for_Up , Score_for_Up = Pvalue_for_Stream_Down , Score_for_Stream_Down
                    next_done = 'no'
                elif next_done == 'no':
                    Upstream_end = self.final_nucleosome[i][0] - 10    
                    Upstream_start = max( Upstream_end - 1000 , beginning_coordinate )    
                    Upstream_mid = ( Upstream_start + Upstream_end ) // 2 // 10 * 10 + 1
                    Stream_Up_left_1000  = max( beginning_coordinate  , Upstream_mid - 1000 )
                    Stream_Up_right_1000 = min( Upstream_mid + 1000 , ending_coordinate     )
                    Stream_Up_left_5000  = max( beginning_coordinate  , Upstream_mid - 5000 )
                    Stream_Up_right_5000 = min( Upstream_mid + 5000 , ending_coordinate     )
                    Stream_Up_left_10000  = max( beginning_coordinate  , Upstream_mid - 10000 )
                    Stream_Up_right_10000 = min( Upstream_mid + 10000 , ending_coordinate     )
                    ####
                    Upstream_start_index = ( Upstream_start - beginning_coordinate ) // 10
                    Upstream_end_index   = ( Upstream_end - beginning_coordinate ) // 10
                    Stream_Up_left_1000_index  = ( Stream_Up_left_1000 - beginning_coordinate ) // 10
                    Stream_Up_right_1000_index = ( Stream_Up_right_1000 - beginning_coordinate ) // 10
                    Stream_Up_left_5000_index  = ( Stream_Up_left_5000 - beginning_coordinate ) // 10
                    Stream_Up_right_5000_index = ( Stream_Up_right_5000 - beginning_coordinate ) // 10
                    Stream_Up_left_10000_index  = ( Stream_Up_left_10000 - beginning_coordinate ) // 10
                    Stream_Up_right_10000_index = ( Stream_Up_right_10000 - beginning_coordinate ) // 10
                    BG_1000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_1000_index , Stream_Up_right_1000_index+1 ) ] ) / (Stream_Up_right_1000_index-Stream_Up_left_1000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_5000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_5000_index , Stream_Up_right_5000_index+1 ) ] ) / (Stream_Up_right_5000_index-Stream_Up_left_5000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_10000_for_Up = sum( [ self.tag_list[index] for index in range( Stream_Up_left_10000_index , Stream_Up_right_10000_index+1 ) ] ) / (Stream_Up_right_10000_index-Stream_Up_left_10000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_for_Up = min( BG_1000_for_Up , BG_5000_for_Up , BG_10000_for_Up )
                    if BG_for_Up < 0:
                        BG_for_Up = 0
                    foreground_for_Up = sum( [ self.tag_list[index] for index in range( Upstream_start_index , Upstream_end_index+1 ) ] )
                    Pvalue_for_Up , Score_for_Up = Poisson_test.less_fast( BG_for_Up , foreground_for_Up )
                Stream_Down_start = self.final_nucleosome[i][1] + 10   
                if self.final_nucleosome[i+1][0] - Stream_Down_start <= 1000:
                    Stream_Down_end = self.final_nucleosome[i+1][0] - 10    
                    next_done = 'yes'
                else:
                    Stream_Down_end = Stream_Down_start + 1000    
                    next_done = 'no'
                Stream_Down_mid = ( Stream_Down_start + Stream_Down_end ) // 2 // 10 * 10 + 1
                Stream_Down_left_1000  = max( beginning_coordinate   , Stream_Down_mid - 1000 )
                Stream_Down_right_1000 = min( Stream_Down_mid + 1000 , ending_coordinate      )
                Stream_Down_left_5000  = max( beginning_coordinate   , Stream_Down_mid - 5000 )
                Stream_Down_right_5000 = min( Stream_Down_mid + 5000 , ending_coordinate      )
                Stream_Down_left_10000  = max( beginning_coordinate   , Stream_Down_mid - 10000 )
                Stream_Down_right_10000 = min( Stream_Down_mid + 10000 , ending_coordinate      )
                ####
                Stream_Down_start_index = ( Stream_Down_start - beginning_coordinate ) // 10
                Stream_Down_end_index   = ( Stream_Down_end - beginning_coordinate ) // 10
                Stream_Down_left_1000_index  = ( Stream_Down_left_1000 - beginning_coordinate ) // 10
                Stream_Down_right_1000_index = ( Stream_Down_right_1000 - beginning_coordinate ) // 10
                Stream_Down_left_5000_index  = ( Stream_Down_left_5000 - beginning_coordinate ) // 10
                Stream_Down_right_5000_index = ( Stream_Down_right_5000 - beginning_coordinate ) // 10
                Stream_Down_left_10000_index  = ( Stream_Down_left_10000 - beginning_coordinate ) // 10
                Stream_Down_right_10000_index = ( Stream_Down_right_10000 - beginning_coordinate ) // 10
                BG_1000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_1000_index , Stream_Down_right_1000_index+1 ) ] ) / (Stream_Down_right_1000_index-Stream_Down_left_1000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_5000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_5000_index , Stream_Down_right_5000_index+1 ) ] ) / (Stream_Down_right_5000_index-Stream_Down_left_5000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_10000_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_left_10000_index , Stream_Down_right_10000_index+1 ) ] ) / (Stream_Down_right_10000_index-Stream_Down_left_10000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_for_Stream_Down = min( BG_1000_for_Stream_Down , BG_5000_for_Stream_Down , BG_10000_for_Stream_Down )
                if BG_for_Stream_Down < 0:
                    BG_for_Stream_Down = 0
                foreground_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_start_index , Stream_Down_end_index+1 ) ] )
                Pvalue_for_Stream_Down , Score_for_Stream_Down = Poisson_test.less_fast( BG_for_Stream_Down , foreground_for_Stream_Down )
            elif i == len( self.final_nucleosome ) - 1:
                if next_done == 'yes':
                    Pvalue_for_Up , Score_for_Up = Pvalue_for_Stream_Down , Score_for_Stream_Down
                    next_done = 'no'
                elif next_done == 'no':
                    Upstream_end = self.final_nucleosome[i][0] - 10   
                    Upstream_start = max( Upstream_end - 1000 , beginning_coordinate )    
                    Upstream_mid = ( Upstream_start + Upstream_end ) // 2 // 10 * 10 + 1
                    Stream_Up_left_1000  = max( beginning_coordinate  , Upstream_mid - 1000 )
                    Stream_Up_right_1000 = min( Upstream_mid + 1000 , ending_coordinate     )
                    Stream_Up_left_5000  = max( beginning_coordinate  , Upstream_mid - 5000 )
                    Stream_Up_right_5000 = min( Upstream_mid + 5000 , ending_coordinate     )
                    Stream_Up_left_10000  = max( beginning_coordinate  , Upstream_mid - 10000 )
                    Stream_Up_right_10000 = min( Upstream_mid + 10000 , ending_coordinate     )
                    Upstream_start_index = ( Upstream_start - beginning_coordinate ) // 10
                    Upstream_end_index   = ( Upstream_end - beginning_coordinate ) // 10
                    Stream_Up_left_1000_index  = ( Stream_Up_left_1000 - beginning_coordinate ) // 10
                    Stream_Up_right_1000_index = ( Stream_Up_right_1000 - beginning_coordinate ) // 10
                    Stream_Up_left_5000_index  = ( Stream_Up_left_5000 - beginning_coordinate ) // 10
                    Stream_Up_right_5000_index = ( Stream_Up_right_5000 - beginning_coordinate ) // 10
                    Stream_Up_left_10000_index  = ( Stream_Up_left_10000 - beginning_coordinate ) // 10
                    Stream_Up_right_10000_index = ( Stream_Up_right_10000 - beginning_coordinate ) // 10
                    BG_1000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_1000_index , Stream_Up_right_1000_index+1 ) ] ) / (Stream_Up_right_1000_index-Stream_Up_left_1000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_5000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_5000_index , Stream_Up_right_5000_index+1 ) ] ) / (Stream_Up_right_5000_index-Stream_Up_left_5000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_10000_for_Up = sum( [ self.tag_list[index] for index in range( Stream_Up_left_10000_index , Stream_Up_right_10000_index+1 ) ] ) / (Stream_Up_right_10000_index-Stream_Up_left_10000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_for_Up = min( BG_1000_for_Up , BG_5000_for_Up , BG_10000_for_Up )
                    if BG_for_Up < 0:
                        BG_for_Up = 0
                    foreground_for_Up = sum( [ self.tag_list[index] for index in range( Upstream_start_index , Upstream_end_index+1 ) ] )
                    Pvalue_for_Up , Score_for_Up = Poisson_test.less_fast( BG_for_Up , foreground_for_Up )
                Stream_Down_start = self.final_nucleosome[i][1] + 10    
                Stream_Down_end = min( Stream_Down_start + 300 , ending_coordinate )    
                Stream_Down_mid = ( Stream_Down_start + Stream_Down_end ) // 2 // 10 * 10 + 1
                Stream_Down_left_1000  = max( beginning_coordinate   , Stream_Down_mid - 1000 )
                Stream_Down_right_1000 = min( Stream_Down_mid + 1000 , ending_coordinate      )
                Stream_Down_left_5000  = max( beginning_coordinate   , Stream_Down_mid - 5000 )
                Stream_Down_right_5000 = min( Stream_Down_mid + 5000 , ending_coordinate      )
                Stream_Down_left_10000  = max( beginning_coordinate   , Stream_Down_mid - 10000 )
                Stream_Down_right_10000 = min( Stream_Down_mid + 10000 , ending_coordinate      )
                Stream_Down_start_index = ( Stream_Down_start - beginning_coordinate ) // 10
                Stream_Down_end_index   = ( Stream_Down_end - beginning_coordinate ) // 10
                Stream_Down_left_1000_index  = ( Stream_Down_left_1000 - beginning_coordinate ) // 10
                Stream_Down_right_1000_index = ( Stream_Down_right_1000 - beginning_coordinate ) // 10
                Stream_Down_left_5000_index  = ( Stream_Down_left_5000 - beginning_coordinate ) // 10
                Stream_Down_right_5000_index = ( Stream_Down_right_5000 - beginning_coordinate ) // 10
                Stream_Down_left_10000_index  = ( Stream_Down_left_10000 - beginning_coordinate ) // 10
                Stream_Down_right_10000_index = ( Stream_Down_right_10000 - beginning_coordinate ) // 10
                BG_1000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_1000_index , Stream_Down_right_1000_index+1 ) ] ) / (Stream_Down_right_1000_index-Stream_Down_left_1000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_5000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_5000_index , Stream_Down_right_5000_index+1 ) ] ) / (Stream_Down_right_5000_index-Stream_Down_left_5000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_10000_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_left_10000_index , Stream_Down_right_10000_index+1 ) ] ) / (Stream_Down_right_10000_index-Stream_Down_left_10000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_for_Stream_Down = min( BG_1000_for_Stream_Down , BG_5000_for_Stream_Down , BG_10000_for_Stream_Down )
                if BG_for_Stream_Down < 0:
                    BG_for_Stream_Down = 0
                foreground_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_start_index , Stream_Down_end_index+1 ) ] )
                Pvalue_for_Stream_Down , Score_for_Stream_Down = Poisson_test.less_fast( BG_for_Stream_Down , foreground_for_Stream_Down )
            else:
                print('\tError ==> Nucleosome: ' , i )
            Score_for_V = ( Score_for_Up + Score_for_Stream_Down ) * 0.5

            self.Log_Pvalue_Peak_dict[i]   = Score_for_peak
            self.Log_Pvalue_Valley_dict[i] = Score_for_V
            sys.stdout.flush()
        print('...... Statistic tests are finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
            

    def write_nucleosome_file(self):
        print('Generating nucleosome collection file ......')
        nucleosome_beginning_ending_file=open(self.outputfile_nucleosome,'w')
        nucleosome_beginning_ending_file.write('chromosome='+self.chromosome+'  total_nucleosome:'+str(len(self.final_nucleosome))+'  coordinate:hg18\n')
        nucleosome_beginning_ending_file.write('1:Chromosome\t2:Inflection_Pair_Beginning\t3:Inflection_Pair_Ending\t4:Length_between_Inflection\t5:Height_of_nucleosome\t6:Area_under_CuStream_Downe\t7:Nucleosome_Index\t8:Beginning_Peak_Region\t9:Ending_Peak_Region\t10:Peak_Region_Length\t11:Physical_Property\n')
        line=0
        for k in range(len(self.final_nucleosome)):
            line=line+1
            nucleosome_beginning_ending_file.write(self.chromosome+'\t')  
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][0])+'\t')    
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][1])+'\t')    
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][1]-self.final_nucleosome[k][0]+10)+'\t')   
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][4])+'\t')    
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][5])+'\t')    
            nucleosome_beginning_ending_file.write('Nucleosome:'+str(line)+'\t')    
            nucleosome_beginning_ending_file.write('BeginningPeak:'+str(self.final_nucleosome[k][6][8]))    
            nucleosome_beginning_ending_file.write('('+str(self.final_nucleosome[k][6][0])+'--'+str(self.final_nucleosome[k][6][1])+'--'+str(self.final_nucleosome[k][6][2])+'--'+str(self.final_nucleosome[k][6][3])+')'+'\t')
            nucleosome_beginning_ending_file.write('EndingPeak:'+str(self.final_nucleosome[k][7][8]))    
            nucleosome_beginning_ending_file.write('('+str(self.final_nucleosome[k][7][0])+'--'+str(self.final_nucleosome[k][7][1])+'--'+str(self.final_nucleosome[k][7][2])+'--'+str(self.final_nucleosome[k][7][3])+')'+'\t')
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][7][3]-self.final_nucleosome[k][6][0]+10)+'\t')   
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][8])+'\n')    
        nucleosome_beginning_ending_file.close()
        print('...... Nucleosome collection file is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Total line of the collection file:  ',line)



class Positioning_improvement(NucDetect):
    def __init__(self,FilesParameters,ConvolutionParameters,threshold):
        self.inputfile = FilesParameters['NucPosition_in']
        self.chromosome = FilesParameters['Chosen_Chromosome_Abbreviation']
        self.chrlength = FilesParameters['Chosen_Chromosome_Length']
        self.chrRecordingListLength = 0
        self.outputfile_wiggle = FilesParameters['outputfile_like_wig']
        self.outputfile_nucleosome = FilesParameters['outputfile_like_bed']
        print('Output wiggle file of nucleosome distribution: ',self.outputfile_wiggle)
        print('Output file of nucleosome collection:          ',self.outputfile_nucleosome)
        self.ConvolutionParameters = ConvolutionParameters
        self.threshold = threshold
        self.score_list = []
        self.tag_list = [ ]
        self.Gaussian_list = []
        self.FDoG_list = []
        self.LoG_list = []
        self.TDoG_list = []
        self.score_table = []   
      
        if FilesParameters[ 'single_or_paired' ] == 's':
            NucDetect.score(self)
        elif FilesParameters[ 'single_or_paired' ] == 'p':
            NucDetect.score_for_Paired_end(self)
        
        NucDetect.Convolution_smoothing(self)
        NucDetect.extremum_detection(self)
        NucDetect.inflection_pairs_detection(self)
        NucDetect.inflection_pairs_midpoint_detection(self)
        NucDetect.secondary_inflection_pairs_detection(self)
        NucDetect.preliminary_nucleosome_position(self)
        NucDetect.collect_and_sort_shoulders(self)
        NucDetect.precise_nucleosome_position(self)
        NucDetect.adjust_border(self)
        NucDetect.filter_results(self)
        NucDetect.record_results(self)
        NucDetect.Test_significance( self )
       
       
        wiggle=open(self.outputfile_wiggle,'w')
        wiggle.write('track type=like_wiggle\nvariableStep chromosome='+self.chromosome+' span=10\n')
        wiggle.write('1:Coordinate\t2:Original_nucleosome_profile\t3:Convolution_smoothing\t4:LoG\t5:Minor_LoG\t6:Tag_accumulation\t7:Nucloesomes\n')
        line=0
        for k in range(len(self.score_list)):
            wiggle.write(str(k*10+1)+'\t'+str(round(self.score_list[k],3))+'\t'+str(round(self.Gaussian_list[k],3))+'\t'+str(round(self.LoG_list[k],3))+'\t'+str(round(self.secondary_LoG_list[k],3))+'\t'+str(self.tag_list[k])+'\t'+str(round(self.score_table[k][4],3))+'\n')
            line=line+1
        wiggle.close()

        nucleosome_beginning_ending_file=open(self.outputfile_nucleosome,'w')
        nucleosome_beginning_ending_file.write('chromosome='+self.chromosome+'  total_nucleosome:'+str(len(self.final_nucleosome))+'\n')
        nucleosome_beginning_ending_file.write( ('\t').join( ['1:Chromosome','2:Start_inflection','3:End_inflection','4:Nucleosome_index','5:Width_between_inflection','6:Peak_height','7:Area_under_cuStream_Downe','8:Physical_property','9:"-log10(Pvalue_of_peak)"','10:"-log10(Pvalue_of_valley)"'] ) + '\n' )
        line=0
        for k in range(len(self.final_nucleosome)):
            line=line+1
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][0])+'\t')    
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][1])+'\t')   
            nucleosome_beginning_ending_file.write('Nucleosome:'+str(line)+'\t')    
            nucleosome_beginning_ending_file.write(str(self.final_nucleosome[k][1]-self.final_nucleosome[k][0]+10)+'\t')    
            nucleosome_beginning_ending_file.write(str(round(self.final_nucleosome[k][4],3))+'\t')    
            nucleosome_beginning_ending_file.write(str(round(self.final_nucleosome[k][5],3))+'\t')    
            if self.final_nucleosome[k][8] == 'Main_nucleosome:isolated':
                Physical_Property = 'MainPeak'
            elif self.final_nucleosome[k][8] == 'Main_nucleosome:integrated':
                Physical_Property = 'MainPeak+Shoulder'
            elif self.final_nucleosome[k][8] == 'Main_nucleosome:integrated:doublet':
                Physical_Property = 'MainPeak:doublet'
            elif ( 'shoulder' in self.final_nucleosome[k][8] ) or ( 'Shoulder' in self.final_nucleosome[k][8] ):
                Physical_Property = 'Shoulder'
            else:
                Physical_Property = 'Other'
            nucleosome_beginning_ending_file.write( Physical_Property +'\t')    ### 8:Physical_Property
            nucleosome_beginning_ending_file.write( str( self.Log_Pvalue_Peak_dict[k]   ) + '\t' )    ### 9:"-log10(Pvalue_of_peak)"
            nucleosome_beginning_ending_file.write( str( self.Log_Pvalue_Valley_dict[k] ) + '\n' )    ### 10:"-log10(Pvalue_of_valley)"
        nucleosome_beginning_ending_file.close()
        print('...... Nucleosome collection file is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))


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



