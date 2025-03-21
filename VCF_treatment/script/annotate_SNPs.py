#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 08:34:27 2020

@author: julien_joseph
modified by Melodie Bastian to read once the genome and write the sequence of interest in a dict and count the stop codon
"""

from math import*
import numpy as np
from random import*
import os
import sys
#import matplotlib.pyplot as plt
import argparse
import time
from Bio import SeqIO

parser = argparse.ArgumentParser(
        description='This program allows your to annotate your SNPs in function if they are'
                    'located in non coding sequence if they are synonymous or non-synonymous'
                    'For more information read the userguide',
        usage='python annotate_SNPs.py -g [REFERENCE GENOME] -gff [ANNOTATION FILE] -snps [SNPS FILE] -o [OUTPUT NAME]',
        formatter_class=argparse.RawTextHelpFormatter, add_help=False)
parser.add_argument("-g", "--genome", help="genome file")
parser.add_argument("-gff", "--annotation_file", help="annotation file")
parser.add_argument("-vcf", "--vcf_file", help="vcf file")
parser.add_argument("-o", "--output", help="the name for the filtered and annotated GFF")
parser.add_argument("-h", "--help", help="help", action="store_true")
args = parser.parse_args()





def annotate_SNPs(fasta,gff,vcf,output):
    if output in os.getcwd():
        os.remove(output)
    t1=time.perf_counter()
    """ function which takes in argument a GFF file, a reference genome, the list of chromosome 
    in the order of the reference genome, the chromosome of interest, a file containing a list of SNPs 
    and the name of the desired output file, and returns a file with the list of SNPs with their annotation
    NC for non coding S for synonymous, and NS for non synonymous"""
    #define the correspondance between codon and amion acid
    codontable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'}
    #define the complementary base of each base pair
    complementary = {'A':'T','T':'A','C':'G','G':'C','a':'t', 't':'a', 'c':'g', 'g':'c','n':'n', 'N':'N' ,'W':'W', 'R':'R', 'Y':'Y', 'K':'K', 'M':'M', 'S':'S','V':'V', 'H':'H', 'D':'D'}
    #create a dictionnary indexes which will attribute a number to each transcript
    fichier=open(gff,'r')
    chro=[]
    #make a list txt with all the elements of the first line
    txt = fichier.readline().split()
    sortie=open(output, 'w')
    while txt!=[]:
        if not txt[0].startswith('#'):
                  chro.append(txt[0])
        del txt
        txt = fichier.readline().split()
    chro=list(set(chro))
    chro.sort()
    if len(chro)==0:
       print ("list of scaffold empty")
       quit()
    print()
    print('List of chromosomes in annotation file',chro)
    print()
    fichier.close()
    scaf=0
    countS=0
    countNS=0
    countNC=0
    countNA=0
    stop=0

    # load seq of interest in a dict
    print("read the genome from the chromosome list")
    chro2seq = dict()
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in chro:
                sequence = 'N' + str(record.seq)#0-1 base
                chro2seq[record.id] = sequence

    for chromosome in chro:
        print(chromosome)
        scaf=scaf+1
        print('scaffold '+str(scaf)+' over '+str(len(chro)))
        fichier=open(gff,'r')
        #make a list txt with all the elements of the first line
        txt = fichier.readline().split('\t')
        indexes={}
        #create a dictionnary transcripts, which will contain for each transcripts, the nucleotide sequence, the phase sequence, and the corresponding genomic position of each nucleotide of the transcript
        transcripts={}
        #create a list genome where in each positions is stored the nucleotide of the reference genome
        #genome=read_genome(fasta,chromosome)
        genome=list(chro2seq[chromosome].upper())
        #create the list transcription position where is stored at each position of the genome the index of the corresponding transcript or - if there is no transcript
        transcript_position=len(genome)*['-']
        #k is the index of the transcripts (index of first transcript=1)
        k=0
        #while the next line is not empty
        while txt!=['']:
            #print(txt)
            if txt[0]==chromosome:
                if txt[1]=='CDS':
                    #read the start column
                    start=int(txt[2])
                    #read the end column
                    end=int(txt[3])
                    #read the phase column
                    phase=txt[5]
                    #read the strand column
                    strand=txt[4]
                    #read the transcript_ID column
                    transcript_ID=txt[6]
                    #if the line has a gene ID and a transcript ID
                    #test wether the transcript is in the dictionnary
                    if transcript_ID not in transcripts:
                        #change the transcript index
                        k=k+1
                        #assign the new index to the transcript
                        indexes[k]=transcript_ID
                        indexes[transcript_ID]=k
                        # add a new transcripts in the dictionnary, first item is the sequence of nucleotide, second the phase sequence, third the corresponding positions, and fourth the strand
                        transcripts[transcript_ID]=[[],[],[],strand]
                    #if the gene is on the + strand
                    if strand=='+':
                        #add the sequence of the exon to the sequence of the transcript
                        transcripts[transcript_ID][0]=transcripts[transcript_ID][0]+genome[start:end+1]
                        #create the list of corresponding phase of the exon
                        if phase=='0':
                            phase_list=floor((end-start+1)/3+1)*[1,2,3]
                            while len(phase_list)!=len(genome[start:end+1]):
                                phase_list.pop(-1)
                        elif phase=='1':
                            phase_list=floor((end-start+1)/3+1)*[3,1,2]
                            while len(phase_list)!=len(genome[start:end+1]):
                                phase_list.pop(-1)
                        elif phase=='2':
                            phase_list=floor((end-start+1)/3+1)*[2,3,1]
                            while len(phase_list)!=len(genome[start:end+1]):
                                phase_list.pop(-1)
                    #if the gene is on the - strand
                    elif strand=='-':
                        #add the complementary sequence of the exon to the sequence of the transcript
                        transcripts[transcript_ID][0]=transcripts[transcript_ID][0]+[complementary[genome[i]] for i in range(start,end+1)]
                        #create the list of corresponding phase of the exon
                        if phase=='0':
                            phase_list=floor((end-start+1)/3+1)*[3,2,1]
                            while len(phase_list)!=len(genome[start:end+1]):
                                phase_list.pop(0)
                        elif phase=='1':
                            phase_list=floor((end-start+1)/3+1)*[2,1,3]
                            while len(phase_list)!=len(genome[start:end+1]):
                                phase_list.pop(0)
                        elif phase=='2':
                            phase_list=floor((end-start+1)/3+1)*[1,3,2]
                            while len(phase_list)!=len(genome[start:end+1]):
                                phase_list.pop(0)
                    for i in range(start,end+1):
                        transcript_position[i]=indexes[transcript_ID]
                    #add the list of corresponding phase of the exon to the transcript
                    transcripts[transcript_ID][1]=transcripts[transcript_ID][1]+phase_list
                    #add a list of the corresponding positions
                    transcripts[transcript_ID][2]=transcripts[transcript_ID][2]+list(range(start,end+1))
            del txt
            #read the next line
            txt = fichier.readline().split('\t')
        #close the GFF
        fichier.close()
        for transc in transcripts:
            #if the strand is -
            if transcripts[transc][3]=='-':
                #reverse the nucleotide sequence so the START codon is at first position
                transcripts[transc][0].reverse()
                #reverse the phase sequence
                transcripts[transc][1].reverse()
                #reverse the position sequence
                transcripts[transc][2].reverse()
        #open the vcf file
        vcf_file=open(vcf,'r')
        #read the first line
        txt = vcf_file.readline().split()
        #open the output file
        #while the line is not empty
        while txt!=[]:
            if txt[0]==chromosome:
                #read the position of the SNP
                position=int(txt[1])
                #print(position)
                #read the reference allele
                ref=txt[3]
                #read the alternative allele
                alt=txt[4]
                #if the SNPs is inside a region in which there is no CDS
                if transcript_position[position]=='-':
                    #annotate the SNPs as non coding (NC)
                    sortie.write(txt[0]+'\t'+txt[1]+'\t'+txt[2]+'\t'+ref+'\t'+alt+'\t'+'NC'+'\t'+'.'+'\t'+'.'+'\t')
                    countNC=countNC+1
                    for i in range(3,len(txt)-1):
                        sortie.write(txt[i]+'\t')
                    sortie.write(txt[-1]+'\n')
                    #print('reference codon . ','alternative codon . ')
                else:
                    #read the transcript in which the SNPs occur
                    transcript=transcripts[indexes[transcript_position[position]]]
                    #read the position in the CDS sequence where the SNPs occur (different from position in the genome)
                    pos_in_transcript=transcript[2].index(position)
                    print(position)
                    #if the SNPs is the first nucleotide of the codon
                    if transcript[1][pos_in_transcript]==1:
                        #if the codon is not full
                        if pos_in_transcript>=len(transcript[1])-2:
                            #assign . to both codons because they are incomplete
                            codonalt='.'
                            codonref='.'
                        #if the codon is full
                        else:
                            #if the phase is correct (1,2,3)
                            if str(transcript[1][pos_in_transcript])+str(transcript[1][pos_in_transcript+1])+str(transcript[1][pos_in_transcript+2])=='123':
                                #read the codon of the reference allele
                                if transcript[0][pos_in_transcript]==genome[position]:
                                    codonref=ref+transcript[0][pos_in_transcript+1]+transcript[0][pos_in_transcript+2]
                                #read the codon of the alternative allele
                                    codonalt=alt+transcript[0][pos_in_transcript+1]+transcript[0][pos_in_transcript+2]
                                else:
                                    codonref=complementary[ref]+transcript[0][pos_in_transcript+1]+transcript[0][pos_in_transcript+2]
                                #read the codon of the alternative allele
                                    codonalt=complementary[alt]+transcript[0][pos_in_transcript+1]+transcript[0][pos_in_transcript+2]
                            #if the phase is incorrect
                            else:
                                #assign . to both codons because the phase is wrong
                                codonalt='.'
                                codonref='.'
                    #same operation as above but for a SNPs on the second nucleotide of the codon
                    elif transcript[1][pos_in_transcript]==2:
                        if pos_in_transcript==len(transcript[1])-1 or pos_in_transcript==0:
                            codonalt='.'
                            codonref='.'
                        else:
                            if str(transcript[1][pos_in_transcript-1])+str(transcript[1][pos_in_transcript])+str(transcript[1][pos_in_transcript+1])=='123':
                                if transcript[0][pos_in_transcript]==genome[position]:
                                    codonref=transcript[0][pos_in_transcript-1]+ref+transcript[0][pos_in_transcript+1]
                                #read the codon of the alternative allele
                                    codonalt=transcript[0][pos_in_transcript-1]+alt+transcript[0][pos_in_transcript+1]
                                else:
                                    codonref=transcript[0][pos_in_transcript-1]+complementary[ref]+transcript[0][pos_in_transcript+1]
                                #read the codon of the alternative allele
                                    codonalt=transcript[0][pos_in_transcript-1]+complementary[alt]+transcript[0][pos_in_transcript+1]
                            else:
                                codonalt='.'
                                codonref='.'
                    #same operation as above but for a SNPs on the second nucleotide of the codon
                    elif transcript[1][pos_in_transcript]==3:
                        if pos_in_transcript<=1:
                            codonalt='.'
                            codonref='.'
                        else:
                            if str(transcript[1][pos_in_transcript]-2)+str(transcript[1][pos_in_transcript-1])+str(transcript[1][pos_in_transcript])=='123':
                                if transcript[0][pos_in_transcript]==genome[position]:
                                    codonref=transcript[0][pos_in_transcript-2]+transcript[0][pos_in_transcript-1]+ref
                                #read the codon of the alternative allele
                                    codonalt=transcript[0][pos_in_transcript-2]+transcript[0][pos_in_transcript-1]+alt
                                else:
                                    codonref=transcript[0][pos_in_transcript-2]+transcript[0][pos_in_transcript-1]+complementary[ref]
                                #read the codon of the alternative allele
                                    codonalt=transcript[0][pos_in_transcript-2]+transcript[0][pos_in_transcript-1]+complementary[alt]
                            else:
                                codonalt='.'
                                codonref='.'
                    #if the codon does not match a codon in the codontable
                    if codonref not in codontable or codonalt not in codontable:
                        #annotate NA to the SNPs because we can't identify the codon in whitch it appears
                        sortie.write(txt[0]+'\t'+txt[1]+'\t'+txt[2]+'\t'+ref+'\t'+alt+'\t'+'NA'+'\t'+codonref+'\t'+codonalt+'\t')
                        for i in range(3,len(txt)-1):
                            sortie.write(txt[i]+'\t')
                        sortie.write(txt[-1]+'\n')
                        #print('reference codon',codonref,'alternative codon',codonalt,'NA')
                        countNA=countNA+1
                    else:
                        if codontable[codonref]=="X":
                            print("ref stop")
                            stop+=1
                            #sys.exit(1)
                        #if the alternative and reference codon are the same
                        if codontable[codonref]==codontable[codonalt]:
                            countS=countS+1
                            #annotate the SNP as synonymous
                            sortie.write(txt[0]+'\t'+txt[1]+'\t'+txt[2]+'\t'+ref+'\t'+alt+'\t'+'S'+'\t'+codonref+'\t'+codonalt+'\t')
                            print('written S')
                            for i in range(3,len(txt)-1):
                                sortie.write(txt[i]+'\t')
                            sortie.write(txt[-1]+'\n')
                            #print('reference codon',codonref,'alternative codon',codonalt,'S')
                        else:
                            countNS=countNS+1
                            #annotate the SNP as non-synonymous
                            sortie.write(txt[0]+'\t'+txt[1]+'\t'+txt[2]+'\t'+ref+'\t'+alt+'\t'+'NS'+'\t'+codonref+'\t'+codonalt+'\t')
                            print('written NS')
                            for i in range(3,len(txt)-1):
                                sortie.write(txt[i]+'\t')
                            sortie.write(txt[-1]+'\n')
                            #print('reference codon',codonref,'alternative codon',codonalt,'NS')
            #delete the line
            del txt
            #read the next line
            txt=vcf_file.readline().split()
        #close the vcf and output files
    tot=countNC+countNA+countS+countNS
    sortie.close()
    vcf_file.close()
    if tot!=0:
        print()
        print()
        print('Unknown: ',countNA/tot*100,'%')
        print('Non-coding: ',countNC/tot*100,'%')
        print('Synonymous: ',countS/tot*100,'%')
        print('Non-synonymous: ',countNS/tot*100,'%')
        print('Stop: ', str(stop))
    t2=time.perf_counter()
    print('time in seconds ',t2-t1)
    #print('Estimated dN/dS ',(countNS/196)/(countS/67))
    return()
        
if args.help:
    print('This program allows your to filter you GFF annotations so that only non overlapping\n'
                    'transcripts remain with an index of confidence regarding its functionality.\n'
                    'For more information read the README.\n'
                    'usage=python3 filter_GFF.py -dirpath [WORKING DIRECTORY PATHWAY] -g [REFERENCE GENOME] -gff [ANNOTATION FILE] -o [OUTPUT NAME] ' )
if not args.help:
    if not args.genome:
        raise Exception('Genome file missing')
    if not args.annotation_file:
        raise Exception('GFF file missing')
    if not args.vcf_file:
        raise Exception('SNPs file missing')
    if not args.output:
        raise Exception('Output name missing')
    annotate_SNPs(args.genome,args.annotation_file,args.vcf_file,args.output)


    




