#Readme.txt

There are some examples in the Squence folder

You can run the algorithm with the following command line:
> python Needleman-Wunsch.py
The default input file is:
'Sequence/Seq1.fasta'	and	'Sequence/Seq2.fasta'

If you want to change the default intput file,You can use the command line:
> python Needleman-Wunsch.py --file1=filename 
For example, you can change the default file1 to 'Squence/Seq3.fasta'
> python Needleman-Wunsch.py --file1=~/Sequence/Seq3.fasta

similarly,you can change the defaule file2:
> python Needleman-Wumsch.py --file2=filename

However, it should be noted that the file format you import must be the same as the instance file, that is, the standard FASTA format file.

In this algorithm, the default scoring system, match score is 1, and mismatch the score is -1, if you want to change the score arrangement, you can use the following command line:
> python Needleman-Wumsch.py --match = num1 --mismatch = num2
For example, you can change the default score system, change the matching score to 2, and mismatch the corresponding score to -2.
>python Needleman-Wumsch.py --match = 2 --mismatch = -2

If you want to ask for help, you can use the following command line, although it does not make much difference:
> python Needleman-Wunsch.py -h
