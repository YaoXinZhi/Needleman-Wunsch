
# file NEddle
# /usr/bin/python
# coding:utf-8

from numpy import *
import numpy as np
import sys
import getopt

def compare(str1,str2):
    
    if (str1 == str2 ):
        return match_score
    else:
        return mismatch_score
    
def loss_function(matrix,i,j,Label_Dir):
    
    #i == j == 1
    a = matrix[i-1,j] + compare(Seq1[i-1],'-')
    b = matrix[i,j-1] + compare(Seq2[j-1],'-')
    c = matrix[i-1,j-1] + compare(Seq1[i-1],Seq2[j-1])
    save_CellLabel(i,j,a,b,c,Label_Dir)
    matrix[i,j] = max(a,b,c)
    
    return matrix

def generate_matrix(file1,file2):

    x=len(file1)+1
    y=len(file2)+1
    distance_matrix = mat(zeros((x,y)))
    Cell_Label={}

    for i in range(1,x):
        distance_matrix[i,0] = distance_matrix[i-1,0]-1    
    for j in range(1,y):
        distance_matrix[0,j] = distance_matrix[0,j-1]-1          
    for i in range(1,x):
        for j in range(1,y):
            distance_matrix = loss_function(distance_matrix,i,j,Cell_Label)

    return distance_matrix ,Cell_Label

def insert(original,new,pos):
    return original[:pos] + new + original[pos:]

def print_matrix(result):
    print ('socore matrix: ')
    row , col = shape(result)
    print()
    print("        ", end = '')
    for i in Seq2:
        print(i, end='\t')
    print()
    for i in range(0, row):
        if (i > 0):
            print(Seq1[i-1], end = '\t')
        else:
            print("  ", end = '')
    
        for j in range(0, col):
            print('%d' % result[i,j], end = '\t')
        print()
        
def string_scoring(str1,str2):
    score=0
    gap=0
    for i in range(len(str1)):  
        if str1[i]==str2[i]:
            score=score+1
        elif str1[i]!=str2[i]:
            score=score-1
            gap=gap+1
            
    return score,gap
        
def save_CellLabel(m,n,a,b,c,matrix_to_return):
    #change global label matrix
    position = '%s%s' % (str(m),str(n))
    matrix_to_return[position]=[]
    #c--Oblique a--Up b--Left
    if max(a,b,c) == a:
        matrix_to_return[position].append('U')
    if max(a,b,c) == b:
        matrix_to_return[position].append('L')
    if max(a,b,c) == c:
        matrix_to_return[position].append('O')
    
    
        
def get_OptimalPath(Scoring_matrix,Label_dir):
    #get OptimalPath base Label Matrix
    SeqA = [Seq1]
    SeqB = [Seq2]
    row ,col = shape(Scoring_matrix)

    while row>1:
        while col>1:
            #There can not be two mismatch
            position = '%s%s' % (str(row-1),str(col-1))
            if 'L' in Label_dir[position]:
                new_seq = SeqA
                new_seq = insert(new_seq[0],'-',row-1)
                if new_seq not in SeqB:
                    SeqA.append(new_seq)
                if len(Label_dir[position]) == 1:
                    col = col-1
            elif 'U' in Label_dir[position]:
                new_seq=SeqB
                new_seq=insert(new_seq[0],'-',col-1)
                if new_seq not in SeqB:
                    SeqB.append(new_seq)
                if len(Label_dir[position]) == 1:
                    row = row-1
            elif 'O' in Label_dir[position]:
                row , col = row-1 ,col-1
            if len(Label_dir[position]) > 1:
                row , col = row-1 ,col-1            
                
    return SeqA,SeqB

def print_result(SeqA,SeqB):
    print ('Results:'+'\n')
    print ('---------'+'    '+'--------------------')

    print ('SequenceA:  '+str([i for i in SeqA]))
    print ('SequenceB:  '+str([i for i in SeqB])+'\n')
    score , gap =string_scoring(SeqA[0],SeqB[0])
    print_Aligments(SeqA[0],SeqB[0],score,gap)
    
    for a in range(1,len(SeqA)):
        for b in range(1,len(SeqB)):
            score , gap = string_scoring(SeqA[a],SeqB[b])
            print_Aligments(SeqA[a],SeqB[b],score,gap)

def print_Aligments(seqA,seqB,score,gap):
    
    print ('Alignment:  '+seqA)
    print ('            '+seqB)
    print ('score:  '+str(score))
    print ('gap:    '+str(gap)+'\n')

    
    
def get_string(file):
    
    with open(file,'r') as f:
        for line in f:
            l = line.replace('\n','')
            if '>' in l:
                pass
            else:
                string = l
                
    return string

def usage():
    print ('please read the file of Readme.txt and to be good lucky !!')

def main():
    Distance_matrix , Cell_Label = generate_matrix(Seq1,Seq2)
    print_matrix(Distance_matrix)   
    SequenceA,SequenceB = get_OptimalPath(Distance_matrix , Cell_Label)
    print_result(SequenceA,SequenceB)

if __name__ == '__main__':
    opts,args = getopt.getopt(sys.argv[1:],'h',['match=','mismatch=','file1=','file2='])
    InputFile1 = '/home/yaoxinzhi/桌面/Bioinformatics/Sequence_Alignment/Needleman-Wunsch/Sequence/Seq1.fasta'
    InputFile2 = '/home/yaoxinzhi/桌面/Bioinformatics/Sequence_Alignment/Needleman-Wunsch/Sequence/Seq2.fasta'
    match_score = 1
    mismatch_score = -1
    
    for op,value in opts:
        if op == '--file1':
            print (value+'111')
            InputFile1 = value
            
        elif op == '--file2':
            InputFile2 = value
            print (value)
        elif op == '-o':
            OutputFile = value
        elif op == '-h':
            usage()
            sys.exit()
        elif op == '--match':
            match_score = float(value)
        elif op == '--mismatch':
            mismatch_score = float(value)
        
    Seq1=get_string(InputFile1)
    Seq2=get_string(InputFile2)
    print ('matrix shape: '+str(len(Seq1))+'*'+str(len(Seq2)))
    main()
    

    
    
        
    
