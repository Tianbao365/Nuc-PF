#!/usr/bin/env python

import math
import re
import pandas as pd
import os
import sys

def collect_and_sort_shoulders(self):
    peak_nuc_grp_shoulder_simple_collection={}
    for one_shoulder in self.shoulder_list:
        inflection_pair_index=one_shoulder[4]    
        peak_index=one_shoulder[-1][-1]    
        if peak_index in self.nuc_grp_dict.keys():
            if peak_index not in peak_nuc_grp_shoulder_simple_collection.keys():
                peak_nuc_grp_shoulder_simple_collection[peak_index]=[]
                peak_nuc_grp_shoulder_simple_collection[peak_index].append(inflection_pair_index)
            elif peak_index in peak_nuc_grp_shoulder_simple_collection.keys():
                peak_nuc_grp_shoulder_simple_collection[peak_index].append(inflection_pair_index)
        elif peak_index not in self.nuc_grp_dict.keys():
            pass

        for peak_index in sorted( peak_nuc_grp_shoulder_simple_collection.keys() ):
            self.shoulders_sorted_with_main_nuc_grps[peak_index]={}   
            main_nuc_grp_in_the_peak=self.nuc_grp_dict[peak_index]
            left_part=[]
            right_part=[]
            for inflection_pair_index in peak_nuc_grp_shoulder_simple_collection[peak_index]:
                if self.inflection_pairs_dict_with_midpoint[inflection_pair_index][0]>self.nuc_grp_dict[peak_index][1]:
                    right_part.append(inflection_pair_index)
                elif self.inflection_pairs_dict_with_midpoint[inflection_pair_index][1]<self.nuc_grp_dict[peak_index][0]:
                    left_part.append(inflection_pair_index)
            if len(left_part)==1:
                self.shoulders_sorted_with_main_nuc_grps[peak_index][-1]=left_part[0]
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
                    self.shoulders_sorted_with_main_nuc_grps[peak_index][(k+1)*(-1)]=left_part[k]
            if len(right_part)==1:
                self.shoulders_sorted_with_main_nuc_grps[peak_index][1]=right_part[0]
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
                    self.shoulders_sorted_with_main_nuc_grps[peak_index][(k+1)]=right_part[k]

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
    else:
        relationship='independent'
    return relationship



class Positioning_improvement(NucGroup):
    def __init__(self,FilesParameters,ConvolutionParameters,threshold):
        self.inputfile = FilesParameters['NucPosition_in']
        self.chromosome = FilesParameters['Chosen_Chromosome_Abbreviation']
        self.chrlength = FilesParameters['Chosen_Chromosome_Length']
        self.chrRecordingListLength = 0
        self.outputfile_wiggle = FilesParameters['outputfile_like_wig']
        self.outputfile_nuc_grp = FilesParameters['outputfile_like_bed']
        print('Output wiggle file of nuc_grp distribution: ',self.outputfile_wiggle)
        print('Output file of nuc_grp collection:          ',self.outputfile_nuc_grp)
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
            NucGroup.score(self)
        elif FilesParameters[ 'single_or_paired' ] == 'p':
            NucGroup.score_for_Paired_end(self)
        
        NucGroup.Convolution_smoothing(self)
        NucGroup.extremum_detection(self)
        NucGroup.inflection_pairs_detection(self)
        NucGroup.inflection_pairs_midpoint_detection(self)
        NucGroup.secondary_inflection_pairs_detection(self)
        NucGroup.preliminary_nuc_grp_position(self)
        NucGroup.collect_and_sort_shoulders(self)
        NucGroup.precise_nuc_grp_position(self)
        NucGroup.adjust_border(self)
        NucGroup.filter_results(self)
        NucGroup.record_results(self)
        NucGroup.Test_significance( self )
       
       
        wiggle=open(self.outputfile_wiggle,'w')
        wiggle.write('track type=like_wiggle\nvariableStep chromosome='+self.chromosome+' span=10\n')
        wiggle.write('1:Coordinate\t2:Original_nuc_grp_profile\t3:Convolution_smoothing\t4:LoG\t5:Minor_LoG\t6:Tag_accumulation\t7:Nucloesomes\n')
        line=0
        for k in range(len(self.score_list)):
            wiggle.write(str(k*10+1)+'\t'+str(round(self.score_list[k],3))+'\t'+str(round(self.Gaussian_list[k],3))+'\t'+str(round(self.LoG_list[k],3))+'\t'+str(round(self.secondary_LoG_list[k],3))+'\t'+str(self.tag_list[k])+'\t'+str(round(self.score_table[k][4],3))+'\n')
            line=line+1
        wiggle.close()

        nuc_grp_beginning_ending_file=open(self.outputfile_nuc_grp,'w')
        nuc_grp_beginning_ending_file.write('chromosome='+self.chromosome+'  total_nuc_grp:'+str(len(self.final_nuc_grp))+'\n')
        nuc_grp_beginning_ending_file.write( ('\t').join( ['1:Chromosome','2:Start_inflection','3:End_inflection','4:nuc_grp_index','5:Width_between_inflection','6:Peak_height','7:Area_under_cuStream_Downe','8:Physical_property','9:"-log10(Pvalue_of_peak)"','10:"-log10(Pvalue_of_valley)"'] ) + '\n' )
        line=0
        for k in range(len(self.final_nuc_grp)):
            line=line+1
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][0])+'\t')    
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][1])+'\t')   
            nuc_grp_beginning_ending_file.write('nuc_grp:'+str(line)+'\t')    
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][1]-self.final_nuc_grp[k][0]+10)+'\t')    
            nuc_grp_beginning_ending_file.write(str(round(self.final_nuc_grp[k][4],3))+'\t')    
            nuc_grp_beginning_ending_file.write(str(round(self.final_nuc_grp[k][5],3))+'\t')    
            if self.final_nuc_grp[k][8] == 'Main_nuc_grp:isolated':
                Physical_Property = 'MainPeak'
            elif self.final_nuc_grp[k][8] == 'Main_nuc_grp:integrated':
                Physical_Property = 'MainPeak+Shoulder'
            elif self.final_nuc_grp[k][8] == 'Main_nuc_grp:integrated:doublet':
                Physical_Property = 'MainPeak:doublet'
            elif ( 'shoulder' in self.final_nuc_grp[k][8] ) or ( 'Shoulder' in self.final_nuc_grp[k][8] ):
                Physical_Property = 'Shoulder'
            else:
                Physical_Property = 'Other'
            nuc_grp_beginning_ending_file.write( Physical_Property +'\t')    ### 8:Physical_Property
            nuc_grp_beginning_ending_file.write( str( self.Log_Pvalue_Peak_dict[k]   ) + '\t' )    ### 9:"-log10(Pvalue_of_peak)"
            nuc_grp_beginning_ending_file.write( str( self.Log_Pvalue_Valley_dict[k] ) + '\n' )    ### 10:"-log10(Pvalue_of_valley)"
        nuc_grp_beginning_ending_file.close()
        print('...... nuc_grp collection file is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
            
        def relation_determination(self,peak_index,shoulders_in_a_peak):
            shoulders_in_a_peak_relation={}
            for position in sorted(shoulders_in_a_peak.keys()):
                inflection_pair_index=shoulders_in_a_peak[position]
                if len(self.inflection_pairs_dict_with_midpoint[inflection_pair_index])>5:
                    midpoints_list_for_pair=self.inflection_pairs_dict_with_midpoint[inflection_pair_index][5:]
                else:
                    midpoints_list_for_pair=[self.inflection_pairs_dict_with_midpoint[inflection_pair_index][:]]
                midpoints_index_list_for_pair=[p[2] for p in midpoints_list_for_pair]+[p[3] for p in midpoints_list_for_pair]
                if position<=-1:
                    fragment_beginning_index = midpoints_index_list_for_pair
                    fragment_ending_index = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position+1]][2]-1
                    gap_length = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position+1]][2]-self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]
                    LeftONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]]
                    RightONE = self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position+1]]
                    independent_or_shifting = preparation(self,fragment_beginning_index,fragment_ending_index,gap_length,'left',peak_index,LeftONE,RightONE)
                elif position>=1:
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
        main_nuc_grp_integrated_with_shoulder_counting=0
        for peak_index in sorted(self.shoulders_sorted_with_main_nuc_grps.keys()):
            shoulders_in_a_peak=self.shoulders_sorted_with_main_nuc_grps[peak_index]
            shoulders_in_a_peak_relation=relation_determination(self,peak_index,shoulders_in_a_peak)
            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
            main_nuc_grp_has_been_treated='no'  
            for position in sorted(shoulders_in_a_peak_relation.keys()):
                shoulder_candidates_counting=shoulder_candidates_counting+1
                if position<-1:
                    if shoulders_in_a_peak_relation[position]=='independent':
                        independent_counting=independent_counting+1   
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
                            one_independent_shoulder['right_index']=self.nuc_grp_dict[peak_index][3]    
                            self.nuc_grp_dict[peak_index][2]=one_independent_shoulder['left_index']
                            self.nuc_grp_dict[peak_index][3]=one_independent_shoulder['right_index']
                            self.nuc_grp_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                            self.nuc_grp_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                            main_nuc_grp_integrated_with_shoulder_counting=main_nuc_grp_integrated_with_shoulder_counting+1
                            self.nuc_grp_dict[peak_index].append('Main_nuc_grp:integrated')
                            main_nuc_grp_has_been_treated='yes'    
                            one_independent_shoulders={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif max(shoulders_in_a_peak_relation.keys())>=1:
                            pass
                elif position==1:
                    if shoulders_in_a_peak_relation[position]=='independent':
                        independent_counting=independent_counting+1  
                        if one_independent_shoulder['left_index']!='none':   
                            one_independent_shoulder['right_index']=self.nuc_grp_dict[peak_index][3]  
                            self.nuc_grp_dict={peak_index:[pairBC, pairEC, pairBI, pairEI, pair_index, 'height', 'area', [peakBC, maxBC, maxEC, peakEC, peakBI, maxBI, maxEI, peakEI, peak_index]], ...}
                            self.nuc_grp_dict[peak_index][2]=one_independent_shoulder['left_index']
                            self.nuc_grp_dict[peak_index][3]=one_independent_shoulder['right_index']
                            self.nuc_grp_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                            self.nuc_grp_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                            main_nuc_grp_integrated_with_shoulder_counting=main_nuc_grp_integrated_with_shoulder_counting+1
                            self.nuc_grp_dict[peak_index].append('Main_nuc_grp:integrated')
                            main_nuc_grp_has_been_treated='yes'    
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif one_independent_shoulder['left_index']=='none':   
                        
                            main_nuc_grp_has_been_treated='yes'
                          
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
                            one_independent_shoulder['left_index']=self.nuc_grp_dict[peak_index][2]    
                        if max(shoulders_in_a_peak_relation.keys())==position:    
                            one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position]][3]    
                            self.nuc_grp_dict[peak_index][2]=one_independent_shoulder['left_index']
                            self.nuc_grp_dict[peak_index][3]=one_independent_shoulder['right_index']
                            self.nuc_grp_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                            self.nuc_grp_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                            main_nuc_grp_integrated_with_shoulder_counting=main_nuc_grp_integrated_with_shoulder_counting+1
                            self.nuc_grp_dict[peak_index].append('Main_nuc_grp:integrated')
                            main_nuc_grp_has_been_treated='yes'    
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'}    
                        elif max(shoulders_in_a_peak_relation.keys())>position:
                            pass
                elif position>1:
                    if shoulders_in_a_peak_relation[position]=='independent':
                        independent_counting=independent_counting+1    
                        if one_independent_shoulder['left_index']!='none':    
                            one_independent_shoulder['right_index']=self.inflection_pairs_dict_with_midpoint[shoulders_in_a_peak[position-1]][3]    
                            if main_nuc_grp_has_been_treated=='yes':    
                                self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                            elif main_nuc_grp_has_been_treated=='no':    
                                self.nuc_grp_dict[peak_index][2]=one_independent_shoulder['left_index']
                                self.nuc_grp_dict[peak_index][3]=one_independent_shoulder['right_index']
                                self.nuc_grp_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                                self.nuc_grp_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                                main_nuc_grp_integrated_with_shoulder_counting=main_nuc_grp_integrated_with_shoulder_counting+1
                                self.nuc_grp_dict[peak_index].append('Main_nuc_grp:integrated')
                                main_nuc_grp_has_been_treated='yes'    
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
                            if main_nuc_grp_has_been_treated=='yes':   
                                self.more_independent_shoulders.append([one_independent_shoulder['left_index']*10+1,one_independent_shoulder['right_index']*10+1,one_independent_shoulder['left_index'],one_independent_shoulder['right_index'],one_independent_shoulder['peak_index']])
                            elif main_nuc_grp_has_been_treated=='no':   
                                self.nuc_grp_dict[peak_index][2]=one_independent_shoulder['left_index']
                                self.nuc_grp_dict[peak_index][3]=one_independent_shoulder['right_index']
                                self.nuc_grp_dict[peak_index][0]=one_independent_shoulder['left_index']*10+1
                                self.nuc_grp_dict[peak_index][1]=one_independent_shoulder['right_index']*10+1
                                main_nuc_grp_integrated_with_shoulder_counting=main_nuc_grp_integrated_with_shoulder_counting+1
                                self.nuc_grp_dict[peak_index].append('Main_nuc_grp:integrated')
                                main_nuc_grp_has_been_treated='yes'   
                            one_independent_shoulder={'peak_index':peak_index, 'left_index':'none', 'right_index':'none'} 
                        elif max(shoulders_in_a_peak_relation.keys())>position:
                            pass

       
        self.nuc_grp_integrated=[]
        shoulder_counting=0
        shoulder_counting_suplimit=len(self.more_independent_shoulders)-1
        peak_index=1
        if len( self.nuc_grp_dict.keys() ) > 0:
            peak_index_suplimit = sorted( list( self.nuc_grp_dict.keys() ) )[-1]  
        else:
            peak_index_suplimit = 0
        total_shoulders=0
        total_main_nuc_grps=0
        isolated_main_nuc_grp_counting=0
        if len( self.more_independent_shoulders ) > 0:
            for bp in range(len(self.score_table)):
                if self.more_independent_shoulders[shoulder_counting][0]==self.score_table[bp][0]:
                    self.nuc_grp_integrated.append([self.more_independent_shoulders[shoulder_counting][0],self.more_independent_shoulders[shoulder_counting][1],self.more_independent_shoulders[shoulder_counting][2],self.more_independent_shoulders[shoulder_counting][3],'height','area',[self.peak_dict[peak_index][0],self.peak_dict[peak_index][1],self.peak_dict[peak_index][2],self.peak_dict[peak_index][3],self.peak_dict[peak_index][4],self.peak_dict[peak_index][5],self.peak_dict[peak_index][6],self.peak_dict[peak_index][7],peak_index]])
                    self.nuc_grp_integrated[-1].append('Independent_shoulder')
                    total_shoulders=total_shoulders+1   
                    if shoulder_counting<shoulder_counting_suplimit:
                        shoulder_counting=shoulder_counting+1
                    else:
                        pass
                else:
                    while peak_index not in self.nuc_grp_dict.keys():
                        peak_index=peak_index+1
                    if self.nuc_grp_dict[peak_index][0]==self.score_table[bp][0]:
                        self.nuc_grp_integrated.append(self.nuc_grp_dict[peak_index])
                        self.nuc_grp_integrated[-1].pop(4)
                        if self.nuc_grp_integrated[-1][-1]=='Main_nuc_grp:integrated':
                            pass
                        elif self.nuc_grp_integrated[-1][-1]!='Main_nuc_grp:integrated':
                            self.nuc_grp_integrated[-1].append('Main_nuc_grp:isolated')
                            isolated_main_nuc_grp_counting=isolated_main_nuc_grp_counting+1
                        total_main_nuc_grps=total_main_nuc_grps+1
                        if peak_index<peak_index_suplimit:
                            peak_index=peak_index+1
                        else:
                            pass
        else:
            for bp in range(len(self.score_table)):
                while peak_index not in self.nuc_grp_dict.keys():
                    peak_index=peak_index+1
                if self.nuc_grp_dict[peak_index][0]==self.score_table[bp][0]:
                    self.nuc_grp_integrated.append(self.nuc_grp_dict[peak_index])
                    self.nuc_grp_integrated[-1].pop(4)
                    if self.nuc_grp_integrated[-1][-1]=='Main_nuc_grp:integrated':
                        pass
                    elif self.nuc_grp_integrated[-1][-1]!='Main_nuc_grp:integrated':
                        self.nuc_grp_integrated[-1].append('Main_nuc_grp:isolated')
                        isolated_main_nuc_grp_counting=isolated_main_nuc_grp_counting+1
                    total_main_nuc_grps=total_main_nuc_grps+1
                    if peak_index<peak_index_suplimit:
                        peak_index=peak_index+1
                    else:
                        pass
        print('\tChecking ==> ',total_main_nuc_grps,'+',total_shoulders,'=',total_main_nuc_grps+total_shoulders)
        self.final_nuc_grp=[]
        for i in range(len(self.nuc_grp_integrated)):    
            self.final_nuc_grp.append(self.nuc_grp_integrated[i][0:6]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][7]])
        print('...... Precise nuc_grp positioning is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))


            
def adjust_border(self):
print('Adjust the border of every inflection pair with secondary inflection detection ......')
def influence_coefficient(self,ii,influence_from_left_or_right):
    if influence_from_left_or_right=='left':
        influence_index=ii-1
        distance=self.nuc_grp_integrated[ii][2]-self.nuc_grp_integrated[influence_index][3]    
        if self.nuc_grp_integrated[ii][-1]=='Main_nuc_grp:integrated' or self.nuc_grp_integrated[ii][-1]=='Main_nuc_grp:isolated':
            center_1=self.nuc_grp_integrated[ii][6][5]
        elif self.nuc_grp_integrated[ii][-1]=='Independent_shoulder':
            center_1=(self.nuc_grp_integrated[ii][2]+self.nuc_grp_integrated[ii][3])/2
        if self.nuc_grp_integrated[influence_index][-1]=='Main_nuc_grp:integrated' or self.nuc_grp_integrated[influence_index][-1]=='Main_nuc_grp:isolated':
            center_2=self.nuc_grp_integrated[influence_index][6][6]
        elif self.nuc_grp_integrated[influence_index][-1]=='Independent_shoulder':
            center_2=(self.nuc_grp_integrated[influence_index][2]+self.nuc_grp_integrated[influence_index][3])/2
        distance_center=center_1-center_2    
    elif influence_from_left_or_right=='right':
        influence_index=ii+1
        distance=self.nuc_grp_integrated[influence_index][2]-self.nuc_grp_integrated[ii][3]   
        if self.nuc_grp_integrated[influence_index][-1]=='Main_nuc_grp:integrated' or self.nuc_grp_integrated[influence_index][-1]=='Main_nuc_grp:isolated':
            center_1=self.nuc_grp_integrated[influence_index][6][5]
        elif self.nuc_grp_integrated[influence_index][-1]=='Independent_shoulder':
            center_1=(self.nuc_grp_integrated[influence_index][2]+self.nuc_grp_integrated[influence_index][3])/2
        if self.nuc_grp_integrated[ii][-1]=='Main_nuc_grp:integrated' or self.nuc_grp_integrated[ii][-1]=='Main_nuc_grp:isolated':
            center_2=self.nuc_grp_integrated[ii][6][6]
        elif self.nuc_grp_integrated[ii][-1]=='Independent_shoulder':
            center_2=(self.nuc_grp_integrated[ii][2]+self.nuc_grp_integrated[ii][3])/2
        distance_center=center_1-center_2    
    temp_height_self=[]
    temp_smoothedheight_self=[]
    for nnp in range(self.nuc_grp_integrated[ii][2],self.nuc_grp_integrated[ii][3]+1):
        temp_height_self.append(self.score_table[nnp][1])
        temp_smoothedheight_self.append(self.score_table[nnp][2])
    height_self=max(temp_height_self)
    area_self=sum(temp_height_self)
    smoothedheight_self=max(temp_smoothedheight_self)
    smoothedarea_self=sum(temp_smoothedheight_self)
    temp_height_influence=[]
    temp_smoothedheight_influence=[]
    for nnp in range(self.nuc_grp_integrated[influence_index][2],self.nuc_grp_integrated[influence_index][3]+1):
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
for i in range(len(self.nuc_grp_integrated)):
    influence_list.append([])
    peak_index=self.nuc_grp_integrated[i][-2]
    if i==0:
        right_effect=influence_coefficient(self,i,'right')
        if right_effect[0]=='yes':
            influence_list[-1].append('right')
        else:
            pass
    elif i>0 and i<len(self.nuc_grp_integrated)-1:
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
    elif i==len(self.nuc_grp_integrated)-1:
        left_effect=influence_coefficient(self,i,'left')
        if left_effect[0]=='yes':
            influence_list[-1].append('left')
        else:
            pass

    if len( self.nuc_grp_dict.keys() ) > 0:
        peak_index_suplimit = sorted( list(self.nuc_grp_dict.keys()) )[-1] 
    else:
        peak_index_suplimit = 0
    for j in range(len(self.nuc_grp_integrated)):
        peak_index=self.nuc_grp_integrated[j][-2][-1]
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
        else:
            pass
    temp_left_index=0
    temp_left=secondary_canditates[0][2]-self.nuc_grp_integrated[j][2]
    temp_right_index=len(secondary_canditates)-1
    temp_right=secondary_canditates[len(secondary_canditates)-1][3]-self.nuc_grp_integrated[j][3]
    for k in range(len(secondary_canditates)):
        if abs(secondary_canditates[k][2]-self.nuc_grp_integrated[j][2])<abs(temp_left):
            temp_left_index=k
            temp_left=secondary_canditates[k][2]-self.nuc_grp_integrated[j][2]
            if 'left' in adjacent_influence:
                if secondary_canditates[k][2]<self.nuc_grp_integrated[j][2] and secondary_canditates[k][2]>(self.nuc_grp_integrated[j-1][3]+self.nuc_grp_integrated[j][2])/2:
                    self.nuc_grp_integrated[j][0]=secondary_canditates[temp_left_index][0]
                    self.nuc_grp_integrated[j][2]=secondary_canditates[temp_left_index][2]
                else:
                    pass
            else:
                pass
        else:
            pass
        if abs(secondary_canditates[len(secondary_canditates)-1-k][3]-self.nuc_grp_integrated[j][3])<abs(temp_right):
            temp_right_index=len(secondary_canditates)-1-k
            temp_right=secondary_canditates[len(secondary_canditates)-1-k][3]-self.nuc_grp_integrated[j][3]
            if 'right' in adjacent_influence:
                if secondary_canditates[len(secondary_canditates)-1-k][3]>self.nuc_grp_integrated[j][3] and secondary_canditates[len(secondary_canditates)-1-k][3]<(self.nuc_grp_integrated[j][3]+self.nuc_grp_integrated[j+1][2])/2:
                    self.nuc_grp_integrated[j][1]=secondary_canditates[temp_right_index][1]
                    self.nuc_grp_integrated[j][3]=secondary_canditates[temp_right_index][3]
                else:
                    pass
            else:
                pass
        else:
            pass

for i in range(len(self.nuc_grp_integrated)):    
    self.final_nuc_grp.append(self.nuc_grp_integrated[i][0:6]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][7]])
print('...... Adjustment is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))



    def filter_results(self):            

        print('Filter and annotate the nuc_grps ......')
        inteStream_Downal_property_list=[]
        for n in range(len(self.nuc_grp_integrated)-1):
            one_inteStream_Downal={'Distance_between_center':'none',  'InteStream_Downal_length':'none',  'Original_InteStream_Downal_Depth':'none',  'Smoothed_InteStream_Downal_Depth':'none',
                          'Left_OriginalHeight':'none',   'Left_OriginalAUC':'none',   'Left_SmoothedHeight':'none',   'Left_SmoothedAUC':'none',
                          'Right_OriginalHeight':'none',  'Right_OriginalAUC':'none',  'Right_SmoothedHeight':'none',  'Right_SmoothedAUC':'none',
                          'Smoothed_InteStream_Downal_Depth_percentage':'none',
                          'HeightRatio':'none',
                          'Merge':'no',
                          'Noise_on_left_or_right':'none'}
            one_inteStream_Downal['Distance_between_center']=(self.nuc_grp_integrated[n+1][0]+self.nuc_grp_integrated[n+1][1])/2-(self.nuc_grp_integrated[n][0]+self.nuc_grp_integrated[n][1])/2
            one_inteStream_Downal['InteStream_Downal_length']=self.nuc_grp_integrated[n+1][0]-self.nuc_grp_integrated[n][1]  
            InteStream_Downal_original_score_fragment=[]
            InteStream_Downal_smoothed_score_fragment=[]
            for nnq in range(self.nuc_grp_integrated[n][3]+1,self.nuc_grp_integrated[n+1][2]):
                InteStream_Downal_original_score_fragment.append(self.score_list[nnq])
                InteStream_Downal_smoothed_score_fragment.append(self.Gaussian_list[nnq])
            one_inteStream_Downal['Original_InteStream_Downal_Depth']=min(InteStream_Downal_original_score_fragment)
            one_inteStream_Downal['Smoothed_InteStream_Downal_Depth']=min(InteStream_Downal_smoothed_score_fragment)

            Left_original_score_fragment=[]
            Left_smoothed_score_fragment=[]
            for nnp in range(self.nuc_grp_integrated[n][2],self.nuc_grp_integrated[n][3]+1):
                Left_original_score_fragment.append(self.score_list[nnp])
                Left_smoothed_score_fragment.append(self.Gaussian_list[nnp])
            one_inteStream_Downal['Left_OriginalHeight']=max(Left_original_score_fragment)
            one_inteStream_Downal['Left_OriginalAUC']=sum(Left_original_score_fragment)
            one_inteStream_Downal['Left_SmoothedHeight']=max(Left_smoothed_score_fragment)
            one_inteStream_Downal['Left_SmoothedAUC']=sum(Left_smoothed_score_fragment)

            Right_original_score_fragment=[]
            Right_smoothed_score_fragment=[]
            for nnp in range(self.nuc_grp_integrated[n+1][2],self.nuc_grp_integrated[n+1][3]+1):
                Right_original_score_fragment.append(self.score_list[nnp])
                Right_smoothed_score_fragment.append(self.Gaussian_list[nnp])
            one_inteStream_Downal['Right_OriginalHeight']=max(Right_original_score_fragment)
            one_inteStream_Downal['Right_OriginalAUC']=sum(Right_original_score_fragment)
            one_inteStream_Downal['Right_SmoothedHeight']=max(Right_smoothed_score_fragment)
            one_inteStream_Downal['Right_SmoothedAUC']=sum(Right_smoothed_score_fragment)
            one_inteStream_Downal['Smoothed_InteStream_Downal_Depth_percentage']=2*one_inteStream_Downal['Smoothed_InteStream_Downal_Depth']/(one_inteStream_Downal['Left_SmoothedHeight']+one_inteStream_Downal['Right_SmoothedHeight'])    
            one_inteStream_Downal['HeightRatio']=min(one_inteStream_Downal['Left_OriginalHeight'],one_inteStream_Downal['Right_OriginalHeight'])/max(one_inteStream_Downal['Left_OriginalHeight'],one_inteStream_Downal['Right_OriginalHeight'])   
 
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
        for n in range(len(self.nuc_grp_integrated)-2):
            if inteStream_Downal_property_list[n]['Right_OriginalHeight']==inteStream_Downal_property_list[n+1]['Left_OriginalHeight'] and inteStream_Downal_property_list[n]['Right_OriginalAUC']==inteStream_Downal_property_list[n+1]['Left_OriginalAUC'] and inteStream_Downal_property_list[n]['Right_SmoothedHeight']==inteStream_Downal_property_list[n+1]['Left_SmoothedHeight'] and inteStream_Downal_property_list[n]['Right_SmoothedAUC']==inteStream_Downal_property_list[n+1]['Left_SmoothedAUC']:
                pass
            else:
                print('ERROR in inteStream_Downal mapping:  ',n,' and ',n+1)


        def merge_and_annotate(self,inteStream_Downal_property_list):
            merging_choice_list=['none']*len(self.nuc_grp_integrated)
            for i in range(len(self.nuc_grp_integrated)):
                if i==0:
                    if inteStream_Downal_property_list[i]['Merge']=='yes':
                        merging_choice_list[i]='right'
                    else:
                        pass
                elif i>0 and i<len(self.nuc_grp_integrated)-1:
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
                elif i==len(self.nuc_grp_integrated)-1:
                    if inteStream_Downal_property_list[i-1]['Merge']=='yes':
                        merging_choice_list[i]='left'
                    else:
                        pass
 
            merging_counting=0
            unmerging_counting=0
            merging_switch='off'
            final_nuc_grp_temp_recorder=[]    
            one_nuc_grp_for_final_nuc_grp_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
            for i in range(len(self.nuc_grp_integrated)):
                if merging_choice_list[i]=='right':
                    if merging_switch=='off':    
                        one_nuc_grp_for_final_nuc_grp_temp_recorder[0]=self.nuc_grp_integrated[i][0]
                        one_nuc_grp_for_final_nuc_grp_temp_recorder[2]=self.nuc_grp_integrated[i][2]
                        merging_switch='on'
                    elif merging_switch=='on':   
                        final_nuc_grp_temp_recorder.append(self.nuc_grp_integrated[i-1][0:6]+[self.nuc_grp_integrated[i-1][6]]+[self.nuc_grp_integrated[i-1][6]]+[self.nuc_grp_integrated[i-1][7]])    ### 记录前一个
                        unmerging_counting=unmerging_counting+1
                        one_nuc_grp_for_final_nuc_grp_temp_recorder[0]=self.nuc_grp_integrated[i][0]
                        one_nuc_grp_for_final_nuc_grp_temp_recorder[2]=self.nuc_grp_integrated[i][2]

                elif merging_choice_list[i]=='none':
                    if merging_switch=='off':
                        final_nuc_grp_temp_recorder.append(self.nuc_grp_integrated[i][0:6]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][7]])    ### 记录这个
                        unmerging_counting=unmerging_counting+1
                    elif merging_switch=='on':
                        final_nuc_grp_temp_recorder.append(self.nuc_grp_integrated[i-1][0:6]+[self.nuc_grp_integrated[i-1][6]]+[self.nuc_grp_integrated[i-1][6]]+[self.nuc_grp_integrated[i-1][7]])    ### 记录前一个
                        final_nuc_grp_temp_recorder.append(self.nuc_grp_integrated[i][0:6]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][7]])    ### 记录这个
                        unmerging_counting=unmerging_counting+2
                        merging_switch='off'
                        one_nuc_grp_for_final_nuc_grp_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
                elif merging_choice_list[i]=='left':
                    if merging_switch=='off':
                        final_nuc_grp_temp_recorder.append(self.nuc_grp_integrated[i][0:6]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][7]])    ### 记录这个
                        unmerging_counting=unmerging_counting+1
                    elif merging_switch=='on':
                        one_nuc_grp_for_final_nuc_grp_temp_recorder[1]=self.nuc_grp_integrated[i][1]
                        one_nuc_grp_for_final_nuc_grp_temp_recorder[3]=self.nuc_grp_integrated[i][3]
                        one_nuc_grp_for_final_nuc_grp_temp_recorder[7]=self.nuc_grp_integrated[i][6]
                        if type(one_nuc_grp_for_final_nuc_grp_temp_recorder[0])==int and type(one_nuc_grp_for_final_nuc_grp_temp_recorder[1])==int and type(one_nuc_grp_for_final_nuc_grp_temp_recorder[2])==int and type(one_nuc_grp_for_final_nuc_grp_temp_recorder[3])==int and len(one_nuc_grp_for_final_nuc_grp_temp_recorder[6])==9 and len(one_nuc_grp_for_final_nuc_grp_temp_recorder[7])==9:
                            if ('nuc_grp' in one_nuc_grp_for_final_nuc_grp_temp_recorder[8]) or ('nuc_grp' in self.nuc_grp_integrated[i][7]):
                                one_nuc_grp_for_final_nuc_grp_temp_recorder[8]='Main_nuc_grp:integrated:doublet'
                            else:
                                one_nuc_grp_for_final_nuc_grp_temp_recorder[8]='Shoulder:doublet'
                            final_nuc_grp_temp_recorder.append(one_nuc_grp_for_final_nuc_grp_temp_recorder)   
                            merging_counting=merging_counting+1
                            merging_switch='off'
                            one_nuc_grp_for_final_nuc_grp_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
                        else:
                            merging_switch='off'
                            one_nuc_grp_for_final_nuc_grp_temp_recorder=['beginning pairBC', 'ending pairEC', 'beginning pairBI', 'ending pairEI', 'height', 'area', '[beginning peak region]', '[ending peak region]', 'main/shoulder']
                            print('Error in merging!')
            print('\tChecking ==> ',merging_counting,'+',unmerging_counting,'=',merging_counting+unmerging_counting)
            final_nuc_grp_temp_recorder_noise_discarded=[]    
            discarded_noise=0
            for i in range(len(final_nuc_grp_temp_recorder)):
                original_score_fragment=[]
                LoG_fragment=[]
                negative_counting=0
                negative_counting_switch='off'
                negative_counting_list=[]
                LoG_Sigma3_length=0
                for nnp in range(final_nuc_grp_temp_recorder[i][2],final_nuc_grp_temp_recorder[i][3]+1):
                    original_score_fragment.append(self.score_table[nnp][1])
                    LoG_fragment.append(self.score_table[nnp][3])
                    if self.score_table[nnp][3]<0:
                        if nnp<final_nuc_grp_temp_recorder[i][3]:
                            negative_counting_switch='on'
                            negative_counting=negative_counting+1
                        elif nnp==final_nuc_grp_temp_recorder[i][3]:
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
                    final_nuc_grp_temp_recorder_noise_discarded.append(final_nuc_grp_temp_recorder[i])
            print('\tChecking ==> Discarded nuc_grps as noise:  ',discarded_noise)
            print('\tChecking ==> Final remaining nuc_grps:  ',len(final_nuc_grp_temp_recorder_noise_discarded))
            return final_nuc_grp_temp_recorder_noise_discarded

        if threshold['filter_switch']=='off':
            self.final_nuc_grp=[]
            for i in range(len(self.nuc_grp_integrated)):    
                self.final_nuc_grp.append(self.nuc_grp_integrated[i][0:6]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][6]]+[self.nuc_grp_integrated[i][7]])
        elif threshold['filter_switch']=='on':
            self.final_nuc_grp=[]
            self.final_nuc_grp=merge_and_annotate(self,inteStream_Downal_property_list)
        print('\tChecking ==> nuc_grps after filtered:  ',len(self.final_nuc_grp))
        print('...... Results are filtered and annotated.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
            
    def Test_significance( self ):
        self.Log_Pvalue_Peak_dict = { }
        self.Log_Pvalue_Valley_dict = { }
        ########
        print('Test significance ......')
        next_done = 'no'
        beginning_coordinate = 1
        ending_coordinate = len( self.tag_list ) * 10 - 9
        for i in range( len( self.final_nuc_grp ) ):
            sys.stdout.write( '\r    nuc_grp: '+str(i+1) )
            LeftInflection  = self.final_nuc_grp[i][0]  
            RightInflection = self.final_nuc_grp[i][1]   
            nuc_grp_mid = ( LeftInflection + RightInflection ) // 2 // 10 * 10 + 1
            Nuc_BG_left_1000  = max( beginning_coordinate  , nuc_grp_mid - 1000 )
            Nuc_BG_right_1000 = min( nuc_grp_mid + 1000 , ending_coordinate     )
            Nuc_BG_left_5000  = max( beginning_coordinate  , nuc_grp_mid - 5000 )
            Nuc_BG_right_5000 = min( nuc_grp_mid + 5000 , ending_coordinate     )
            Nuc_BG_left_10000  = max( beginning_coordinate   , nuc_grp_mid - 10000 )
            Nuc_BG_right_10000 = min( nuc_grp_mid + 10000 , ending_coordinate      )
            BG_1000_for_peak  = sum( [ self.tag_list[index] for index in range( Nuc_BG_left_1000_index , Nuc_BG_right_1000_index+1 ) ] ) / (Nuc_BG_right_1000_index-Nuc_BG_left_1000_index+1) * (RightInflection_index-LeftInflection_index+1)
            BG_5000_for_peak  = sum( [ self.tag_list[index] for index in range( Nuc_BG_left_5000_index , Nuc_BG_right_5000_index+1 ) ] ) / (Nuc_BG_right_5000_index-Nuc_BG_left_5000_index+1) * (RightInflection_index-LeftInflection_index+1)
            BG_10000_for_peak = sum( [ self.tag_list[index] for index in range( Nuc_BG_left_10000_index , Nuc_BG_right_10000_index+1 ) ] ) / (Nuc_BG_right_10000_index-Nuc_BG_left_10000_index+1) * (RightInflection_index-LeftInflection_index+1)
            BG_for_peak = max( BG_1000_for_peak , BG_5000_for_peak , BG_10000_for_peak )
            foreground_for_peak = sum( [ self.tag_list[index] for index in range( LeftInflection_index , RightInflection_index+1 ) ] )
            Pvalue_for_peak , Score_for_peak = Poisson_test.greater_fast( BG_for_peak , foreground_for_peak )
            if ( i >= 0 ) and ( i < len( self.final_nuc_grp ) - 1 ):
                if next_done == 'yes':
                    Pvalue_for_Up , Score_for_Up = Pvalue_for_Stream_Down , Score_for_Stream_Down
                    next_done = 'no'
                elif next_done == 'no':
                    Upstream_end = self.final_nuc_grp[i][0] - 10    
                    Upstream_start = max( Upstream_end - 1000 , beginning_coordinate )    
                    Upstream_mid = ( Upstream_start + Upstream_end ) // 2 // 10 * 10 + 1
                    Stream_Up_left_1000  = max( beginning_coordinate  , Upstream_mid - 1000 )
                    Stream_Up_right_1000 = min( Upstream_mid + 1000 , ending_coordinate     )
                    Stream_Up_left_5000  = max( beginning_coordinate  , Upstream_mid - 5000 )
                    Stream_Up_right_5000 = min( Upstream_mid + 5000 , ending_coordinate     )
                    Stream_Up_left_10000  = max( beginning_coordinate  , Upstream_mid - 10000 )
                    Stream_Up_right_10000 = min( Upstream_mid + 10000 , ending_coordinate)
                    BG_1000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_1000_index , Stream_Up_right_1000_index+1 ) ] ) / (Stream_Up_right_1000_index-Stream_Up_left_1000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_5000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_5000_index , Stream_Up_right_5000_index+1 ) ] ) / (Stream_Up_right_5000_index-Stream_Up_left_5000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_10000_for_Up = sum( [ self.tag_list[index] for index in range( Stream_Up_left_10000_index , Stream_Up_right_10000_index+1 ) ] ) / (Stream_Up_right_10000_index-Stream_Up_left_10000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_for_Up = min( BG_1000_for_Up , BG_5000_for_Up , BG_10000_for_Up )
                    if BG_for_Up < 0:
                        BG_for_Up = 0
                    foreground_for_Up = sum( [ self.tag_list[index] for index in range( Upstream_start_index , Upstream_end_index+1 ) ] )
                    Pvalue_for_Up , Score_for_Up = Poisson_test.less_fast( BG_for_Up , foreground_for_Up )
                Stream_Down_start = self.final_nuc_grp[i][1] + 10   
                if self.final_nuc_grp[i+1][0] - Stream_Down_start <= 1000:
                    Stream_Down_end = self.final_nuc_grp[i+1][0] - 10    
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
            elif i == len( self.final_nuc_grp ) - 1:
                if next_done == 'yes':
                    Pvalue_for_Up , Score_for_Up = Pvalue_for_Stream_Down , Score_for_Stream_Down
                    next_done = 'no'
                elif next_done == 'no':
                    Upstream_end = self.final_nuc_grp[i][0] - 10   
                    Upstream_start = max( Upstream_end - 1000 , beginning_coordinate )    
                    Upstream_mid = ( Upstream_start + Upstream_end ) // 2 // 10 * 10 + 1
                    Stream_Up_left_1000  = max( beginning_coordinate  , Upstream_mid - 1000 )
                    Stream_Up_right_1000 = min( Upstream_mid + 1000 , ending_coordinate     )
                    Stream_Up_left_5000  = max( beginning_coordinate  , Upstream_mid - 5000 )
                    Stream_Up_right_5000 = min( Upstream_mid + 5000 , ending_coordinate     )
                    Stream_Up_left_10000  = max( beginning_coordinate  , Upstream_mid - 10000 )
                    Stream_Up_right_10000 = min( Upstream_mid + 10000 , ending_coordinate     )
                    BG_1000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_1000_index , Stream_Up_right_1000_index+1 ) ] ) / (Stream_Up_right_1000_index-Stream_Up_left_1000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_5000_for_Up  = sum( [ self.tag_list[index] for index in range( Stream_Up_left_5000_index , Stream_Up_right_5000_index+1 ) ] ) / (Stream_Up_right_5000_index-Stream_Up_left_5000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_10000_for_Up = sum( [ self.tag_list[index] for index in range( Stream_Up_left_10000_index , Stream_Up_right_10000_index+1 ) ] ) / (Stream_Up_right_10000_index-Stream_Up_left_10000_index+1) * (Upstream_end_index-Upstream_start_index+1)
                    BG_for_Up = min( BG_1000_for_Up , BG_5000_for_Up , BG_10000_for_Up )
                    if BG_for_Up < 0:
                        BG_for_Up = 0
                    foreground_for_Up = sum( [ self.tag_list[index] for index in range( Upstream_start_index , Upstream_end_index+1 ) ] )
                    Pvalue_for_Up , Score_for_Up = Poisson_test.less_fast( BG_for_Up , foreground_for_Up )
                Stream_Down_start = self.final_nuc_grp[i][1] + 10    
                Stream_Down_end = min( Stream_Down_start + 300 , ending_coordinate )    
                Stream_Down_mid = ( Stream_Down_start + Stream_Down_end ) // 2 // 10 * 10 + 1
                Stream_Down_left_1000  = max( beginning_coordinate   , Stream_Down_mid - 1000 )
                Stream_Down_right_1000 = min( Stream_Down_mid + 1000 , ending_coordinate      )
                Stream_Down_left_5000  = max( beginning_coordinate   , Stream_Down_mid - 5000 )
                Stream_Down_right_5000 = min( Stream_Down_mid + 5000 , ending_coordinate      )
                Stream_Down_left_10000  = max( beginning_coordinate   , Stream_Down_mid - 10000 )
                Stream_Down_right_10000 = min( Stream_Down_mid + 10000 , ending_coordinate      )

                BG_1000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_1000_index , Stream_Down_right_1000_index+1 ) ] ) / (Stream_Down_right_1000_index-Stream_Down_left_1000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_5000_for_Stream_Down  = sum( [ self.tag_list[index] for index in range( Stream_Down_left_5000_index , Stream_Down_right_5000_index+1 ) ] ) / (Stream_Down_right_5000_index-Stream_Down_left_5000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_10000_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_left_10000_index , Stream_Down_right_10000_index+1 ) ] ) / (Stream_Down_right_10000_index-Stream_Down_left_10000_index+1) * (Stream_Down_end_index-Stream_Down_start_index+1)
                BG_for_Stream_Down = min( BG_1000_for_Stream_Down , BG_5000_for_Stream_Down , BG_10000_for_Stream_Down )
                if BG_for_Stream_Down < 0:
                    BG_for_Stream_Down = 0
                foreground_for_Stream_Down = sum( [ self.tag_list[index] for index in range( Stream_Down_start_index , Stream_Down_end_index+1 ) ] )
                Pvalue_for_Stream_Down , Score_for_Stream_Down = Poisson_test.less_fast( BG_for_Stream_Down , foreground_for_Stream_Down )
            else:
                print('\tError ==> nuc_grp: ' , i )
            Score_for_V = ( Score_for_Up + Score_for_Stream_Down ) * 0.5

            self.Log_Pvalue_Peak_dict[i]   = Score_for_peak
            self.Log_Pvalue_Valley_dict[i] = Score_for_V
            sys.stdout.flush()
        print('...... Statistic tests are finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
            

    def write_nuc_grp_file(self):
        print('Generating nuc_grp collection file ......')
        nuc_grp_beginning_ending_file=open(self.outputfile_nuc_grp,'w')
        nuc_grp_beginning_ending_file.write('chromosome='+self.chromosome+'  total_nuc_grp:'+str(len(self.final_nuc_grp))+'  coordinate:hg18\n')
        nuc_grp_beginning_ending_file.write('1:Chromosome\t2:Inflection_Pair_Beginning\t3:Inflection_Pair_Ending\t4:Length_between_Inflection\t5:Height_of_nuc_grp\t6:Area_under_CuStream_Downe\t7:nuc_grp_Index\t8:Beginning_Peak_Region\t9:Ending_Peak_Region\t10:Peak_Region_Length\t11:Physical_Property\n')
        line=0
        for k in range(len(self.final_nuc_grp)):
            line=line+1
            nuc_grp_beginning_ending_file.write(self.chromosome+'\t')  
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][0])+'\t')    
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][1])+'\t')    
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][1]-self.final_nuc_grp[k][0]+10)+'\t')   
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][4])+'\t')    
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][5])+'\t')    
            nuc_grp_beginning_ending_file.write('nuc_grp:'+str(line)+'\t')    
            nuc_grp_beginning_ending_file.write('BeginningPeak:'+str(self.final_nuc_grp[k][6][8]))    
            nuc_grp_beginning_ending_file.write('('+str(self.final_nuc_grp[k][6][0])+'--'+str(self.final_nuc_grp[k][6][1])+'--'+str(self.final_nuc_grp[k][6][2])+'--'+str(self.final_nuc_grp[k][6][3])+')'+'\t')
            nuc_grp_beginning_ending_file.write('EndingPeak:'+str(self.final_nuc_grp[k][7][8]))    
            nuc_grp_beginning_ending_file.write('('+str(self.final_nuc_grp[k][7][0])+'--'+str(self.final_nuc_grp[k][7][1])+'--'+str(self.final_nuc_grp[k][7][2])+'--'+str(self.final_nuc_grp[k][7][3])+')'+'\t')
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][7][3]-self.final_nuc_grp[k][6][0]+10)+'\t')   
            nuc_grp_beginning_ending_file.write(str(self.final_nuc_grp[k][8])+'\n')    
        nuc_grp_beginning_ending_file.close()
        print('...... nuc_grp collection is finished.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime()))
        print('\tChecking ==> Total line of the collection file:  ',line)





