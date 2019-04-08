
import subprocess

import argparse

desc = '''
-----------------------------------------------------------------------------------------


                        Hello, welcome to pre-mirRNA-plotting!

   We've developed this Python script to facilitate the high-throughput generation of
secondary structure prediction images of miRNA precursors using RNAfold and RNAplot.
One of our main intentions as well is to highlight the position of the miRNAs sequences
within the precursor, as it is an important criteria for miRNA selection and filtering
after the prediction was realized.

   The script will generate a folder called 'prediction_script', where your data will be.
To any valid input files given, a folder will be created with it's name and inside each
one there are going to be another two folders: 'colored_structures' and 'raw_folding'.
One has the secondary structures images with the miRNAs sequences highlighted and the 
other has the TXT files containing the prediction by RNAfold and the uncolored image of
the pre-miRNA structure. The miRNAs are highlighted in green (5p) and red (3p).

   Please read the help info below to have more details about the parameters that best fit
your data and intentions!

-----------------------------------------------------------------------------------------'''

epi = '''


---------------------------------------------------------------


Thanks for using our script! 

This was the first script created in Python by our group so we'd really appreciate if you have any comments, complaints, 
doubts or suggestions. Please don't hesitate to contact our main responsible for this project at igorpaim8@gmail.com.

Have a nice work and let's keep making science evolve!'''

input_help = '''
This script accepts tab-separated text files containing
possibly 4 columns:

1) Sequence ID (optional): Some sort of annotation or ID
information about the sequence (e.g. 'ath-miR-171' or 'seq1');

2) Precursor sequence: The pre-miRNA sequence.

3) miRNA1 sequence: One of the miRNAs sequences.

4) miRNA2 sequence: The other miRNA sequence.


File format examples:


| pre-miRNA |   miRNA1   |  miRNA2  |

>>>>>>>>>>>>>>>>>>>>>> OR <<<<<<<<<<<<<<<<<<<<<<<

|     ID    |  pre-miRNA |  miRNA1  |   miRNA2  |


There is no problem if you don't have both miRNAs sequences; 
you can inform just one and the script will work just fine.
The highlighted miRNA sequence will be in red in that case. 

Check out the parameters descriptions to have more information
about the arguments that best fit your data and intentions.

Also, you can inform more than one file. We suggest that the
filenames should be should somehow informative because they will
be used for naming the generated folder and written with the
file format in the end, such as 'homo_sapiens.txt' or 'hsa.txt'.'''

input_opt_help = '''
Inform the names of the files that you wish to use.\n\n'''

annot_help = '''    T or F (default is False).

Informs if you have some sort of sequence ID, such as a miRNA
family annotation (e.g.'ath-miRNA-171', 'seq1'), necessarily 
on the first column, so that the generated image files can be 
named according to that ID. The default is FALSE and the script
will generate names for the files like 'miRNA-precursor_0' onward.\n\n'''

extra_help = '''    T or F (default is True).

Provides an additional file with basic data about the informed 
pre-miRNAs: their individual length, max and min length, individual
MFE, max and min MFE, and a boxplot of all that data.\n\n'''

usage_text = ' python3 premiRNA-plotting [OPTIONS] [INPUTS]'

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=desc, epilog=epi, add_help=False, usage=usage_text)

parser._positionals.title = 'Positionals arguments'

parser._optionals.title = 'Optional arguments'

parser.add_argument('input', nargs='*', help=input_help, metavar='(INPUT)')

parser.add_argument('-h', '--help', default=argparse.SUPPRESS, action='help', help='Show this help message and exit.')

parser.add_argument('-i', '--input', help=input_opt_help, nargs='*', metavar='\b', dest='inputopt')

parser.add_argument('-a', '--annotation', default='F', help=annot_help, choices=['T', 'F'], metavar='\b')

parser.add_argument('-e', '--extra_info', default=True, help=extra_help, choices=['T', 'F'], metavar='\b')

args = parser.parse_args()

annot = True if args.annotation == 'T' else False

extra = args.extra_info

inputs = args.inputopt if args.input == [] else args.input

subprocess.run('mkdir prediction_script/', shell=True)

pwd = subprocess.run('pwd', shell=True, stdout=subprocess.PIPE, universal_newlines=True).stdout[:-1]

def isin(premirna, mirna, i, filename):
    if mirna in premirna:
        pass
    else:
        print()
        print('''Error: Unfortunately sequence {} at line {} in file '{}' could not be found in {}..., please correct this!'''.format(mirna, str(i+1), filename, premirna[:20]))
        exit()


def repetition(premirna, mirna, i):
    if premirna.count(mirna) > 1:
        print()
        print('Warning! miRNA sequence {} at line {} was found more than once in the precursor sequence, so its last occurrence will be used for the plotting, be aware!'.format(mirna, str(1+i)))


def pos(premirna, mirna1, mirna2):
    if mirna2 == '':
        return [premirna.find(mirna1) + 1, premirna.find(mirna1) + len(mirna1)]
    elif premirna.count(mirna1) == 1:
        return [premirna.find(mirna1) + 1, premirna.find(mirna1) + len(mirna1), premirna.rfind(mirna2) + 1, premirna.rfind(mirna2) + len(mirna2)]
    elif premirna.count(mirna2) == 1:
        return [premirna.find(mirna2) + 1, premirna.find(mirna2) + len(mirna2), premirna.rfind(mirna1) + 1, premirna.rfind(mirna1) + len(mirna1)]


def folding(filename):
    name = filename.split('/')[-1][:-4]
    subprocess.run('mkdir prediction_script/{}'.format(name), shell=True)
    subprocess.run('mkdir prediction_script/{}/foldings'.format(name), shell=True)
    subprocess.run('mkdir prediction_script/{}/colored_structures'.format(name), shell=True)
    with open(filename) as arc:
        lines = arc.readlines()
        for i in range(len(lines)):
                line = lines[i][:-1].split('\t')
                line = [f.upper().replace(' ', '') for f in line]

                # First, defining the existing variables:

                if annot == True:
                    annotation = line[0]
                    precursor = line[1]
                    if len(line) == 4:
                        mirna1 = line[2]
                        mirna2 = line[3]
                    else:
                        mirna1 = line[2]
                else:
                    precursor = line[0]
                    mirna1 = line[1]
                    if len(line) > 2:
                        mirna2 = line[2]

                # Now, checking if there are repeated sequences or if the mirnas are not found in the pre-miRNA

                if annot == True and ({'C', 'T', 'A', 'G'} == set(annotation) or {'C', 'U', 'A', 'G'} == set(annotation)):
                    print('''The informed parameters don't match the file formatting''')
                    exit()

                isin(precursor, mirna1, i, name)

                repetition(precursor, mirna1, i)

                if 'mirna2' in locals():
                    isin(precursor, mirna2, i, name)
                    repetition(precursor, mirna2, i)


                print('\n','######## Data check complete for {}'.format(filename.split('/')[-1]), '\n')

                mirname = annotation if annot == True else 'precursor_{}'.format(i)

                # Run RNAfold to obtain the secondary structure prediction files:

                rnafoldout = open('prediction_script/{}/foldings/{}_fold.txt'.format(name, mirname), 'w')
                rnafold = subprocess.Popen('RNAfold ', stdin=subprocess.PIPE, stdout=rnafoldout, shell=True, universal_newlines=True, cwd='prediction_script/{}/foldings'.format(name))
                rnafold.communicate('>{}\n{}'.format(mirname, precursor))

                # Run RNAplot and color the miRNA sequence within the precursor

                if 'mirna2' not in locals():
                    [init1, fin1] = pos(precursor, mirna1, '')
                else:
                    [init1, fin1, init2, fin2] = pos(precursor, mirna1, mirna2)

                rnaplot = subprocess.Popen(['RNAplot --pre "{} {} 10 {} omark {}"'.format(init1, fin1, 'RED' if (fin1 < (len(precursor)/2.0)) else 'GREEN', '' if 'mirna2' not in locals() else '/ {} {} 10 {} omark'.format(init2, fin2, 'RED' if (fin2 < (len(precursor)/2.0)) else 'GREEN'))], stdin=subprocess.PIPE, cwd='prediction_script/{}/colored_structures'.format(name), shell=True, universal_newlines=True)
                rnaplot.communicate(open('prediction_script/{}/foldings/{}_fold.txt'.format(name, mirname)).read().split(' ')[0])
                subprocess.run('gs -q -dSAFER -dBATCH -dNOPAUSE -sPAPERSIZE=a4 -r200 -sDEVICE=pngalpha -sOutputFile={}.png {}_ss.ps'.format(mirname, mirname), shell=True, cwd='prediction_script/{}/colored_structures'.format(name))
                print(''' # '{}' file generation succesfull'''.format(mirname))

                if 'mirna2' in locals():
                    del mirna2

    subprocess.run('rm *.ps', shell=True, cwd='prediction_script/{}/colored_structures'.format(name))


for file in inputs:
    folding(pwd + '/' + file)



