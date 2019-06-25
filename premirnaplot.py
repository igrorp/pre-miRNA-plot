
import subprocess

import argparse

import matplotlib.pyplot as plt

import os

import numpy as np

from sklearn.linear_model import LinearRegression

desc = '''
-----------------------------------------------------------------------------------------


                        Hello, welcome to pre-mirRNA-plotting!

   We've developed this Python program to facilitate the high-throughput generation of
secondary structure prediction images of miRNA precursors using RNAfold and RNAplot.
One of our main intentions as well is to highlight the position of the miRNAs sequences
within the precursor, as it is an important criteria for miRNA selection and filtering
after the prediction was done.

   The program will generate a folder called 'premirnaplot', where your data will be.
To any valid input file given, a folder will be created with it's name and inside each
one there are going to be another two folders: 'colored_structures' and 'raw_folding'.
One has the secondary structures images with the miRNAs sequences highlighted and the 
other has the TXT files containing the prediction by RNAfold and the uncolored image of
the pre-miRNA structure. The miRNAs are highlighted in green (5p) and red (3p).

   Please read the help info below to have more details about the parameters that best fit
your data and intentions.

-----------------------------------------------------------------------------------------'''

epi = '''


---------------------------------------------------------------


Thanks for using our program! 

This was the first script created in Python by our group so we'd really appreciate if you have any comments, complaints, 
doubts or suggestions. Please don't hesitate to contact our main responsible for this project at igorpaim8@gmail.com.

Have a nice work and let's keep making science evolve!'''

input_help = '''
This program accepts tab-separated text files containing
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
you can inform just one and the program will work just fine. 

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
named according to that ID. The default is FALSE and the program
will generate names for the files like 'miRNA-precursor_0' onward.\n\n'''

extra_help = '''    T or F (default is True).

Provides additional files with basic data about the informed 
pre-miRNAs: one boxplot with their minimum free energy values,
as calculated by RNAfold, and another boxplot of their length.\n\n'''

color_help = '''
COLOR COLOR or RGBCODE RGBCODE (default is GREEN and RED)

You can choose what colors to paint the miRNA sequence within
the precursor. RNAplot provides the predefined colors BLACK,
RED, GREEN, WHITE and BLUE, but you can also inform an RGB code.
Always provide the 5p and 3p, respectively, if you have two miRNA
sequences.

Ex: '-c BLUE GREEN' for blue 5p and green 3p
Ex: '-c 255 255 0 153 0 204' for yellow 5p and purple 3p

'''

usage_text = ' python3+ premiRNA-plotting.py [OPTIONS] ... [INPUTS] ...'

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=desc, epilog=epi, add_help=False, usage=usage_text)

parser._positionals.title = 'Positionals arguments'

parser._optionals.title = 'Optional arguments'

parser.add_argument('input', nargs='*', help=input_help, metavar='(INPUT)')

parser.add_argument('-h', '--help', default=argparse.SUPPRESS, action='help', help='Shows this help message and exit.\n')

parser.add_argument('-i', '--input', help=input_opt_help, nargs='+', metavar='', dest='inputopt')

parser.add_argument('-a', '--annotation', default='F', help=annot_help, choices=['T', 'F'], metavar='\b')

parser.add_argument('-e', '--extra_info', default='T', help=extra_help, choices=['T', 'F'], metavar='\b')

parser.add_argument('-c', '--colors', nargs='+', metavar='', help=color_help, default=['RED', 'GREEN'])

args = parser.parse_args()

annot = True if args.annotation == 'T' else False

extra = True if args.extra_info == 'T' else False

inputs = args.inputopt if args.input == [] else args.input

colors = args.colors

if not inputs:
    print('\n\nError! No input files were given\n\n')
    quit()

if len(colors) == 1:
    colors.append(colors[0])
if len(colors) == 2:
    fc, tc = colors[0].upper(), colors[1].upper()
if len(colors) == 3:
    fc, tc = '{} {} {}'.format(colors[0], colors[1], colors[2]), '{} {} {}'.format(colors[0], colors[1], colors[2])
if len(colors) == 6:
    fc, tc = '{} {} {}'.format(colors[0], colors[1], colors[2]), '{} {} {}'.format(colors[3], colors[4], colors[5])


subprocess.run('mkdir premirnaplot/', shell=True)

def initial_check(filename):
    name = filename.split('/')[-1][:-4]
    with open(filename) as arc:
        lines = arc.readlines()
        for i in range(len(lines)):
                line = lines[i][:-1].split('\t')
                line = [f.upper().replace(' ', '') for f in line]

                if annot:
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

                isin(precursor, mirna1, i, name)

                repetition(precursor, mirna1, i)

                if 'mirna2' in locals():
                    isin(precursor, mirna2, i, name)
                    repetition(precursor, mirna2, i)
                    del mirna2

        print('\n','######## Data check complete for {}'.format(filename.split('/')[-1]), '\n')


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
    mfelst = []
    sizelst = []
    name = filename.split('/')[-1][:-4]
    subprocess.run('mkdir premirnaplot/{} premirnaplot/{}/foldings premirnaplot/{}/colored_structures'.format(name, name, name), shell=True)
    with open(filename) as arc:
        os.chdir('premirnaplot/{}/'.format(name))

        with open('precursor_data.txt', 'a') as data:
            data.writelines(['miRNA name', '\t', 'Sequence', '\t', 'MFE', '\t', 'Lenght', '\n'])

        lines = arc.readlines()
        for i in range(len(lines)):
                line = lines[i][:-1].split('\t')
                line = [f.upper().replace(' ', '') for f in line]

                if annot:
                    annotation = line[0].lower()
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

                mirname = annotation if annot == True else 'precursor_{}'.format(i)

                rnafoldout = open('foldings/{}_fold.txt'.format(mirname), 'w')
                rnafold = subprocess.Popen('RNAfold ',  stdin=subprocess.PIPE,
                                                        stdout=rnafoldout,
                                                        shell=True,
                                                        universal_newlines=True,
                                                        cwd='foldings/')
                rnafold.communicate('>{}\n{}'.format(mirname, precursor))

                if 'mirna2' not in locals():
                    [init1, fin1] = pos(precursor, mirna1, '')
                else:
                    [init1, fin1, init2, fin2] = pos(precursor, mirna1, mirna2)

                p1 = fc if (fin1 < (len(precursor) / 2.0)) else tc
                p2 = '' if 'mirna2' not in locals() else '/ {} {} 10 {} omark'.format(init2, fin2, fc if (fin2 < (len(precursor) / 2.0)) else tc)

                rnaplot = subprocess.Popen(['RNAplot --pre "{} {} 10 {} omark {}"'.format(init1, fin1, p1, p2)],
                                                                                                                stdin=subprocess.PIPE,
                                                                                                                cwd='colored_structures/',
                                                                                                                shell=True,
                                                                                                                universal_newlines=True)
                rnaplot.communicate(open('foldings/{}_fold.txt'.format(mirname)).read().split(' ')[0])


                if extra:
                    mfe = open('foldings/{}_fold.txt'.format(mirname)).read().split(' ')[-1][1:-2]
                    mfelst.append(float(mfe))
                    sizelst.append(float(len(precursor)))
                    with open('precursor_data.txt', 'a') as data:
                        data.writelines([mirname, '\t', precursor, '\t', mfe, '\t', str(len(precursor)), '\n'])


                subprocess.run('gs -q -dSAFER -dBATCH -dNOPAUSE -sPAPERSIZE=a4 -r200 -sDEVICE=pngalpha -sOutputFile={}.png {}_ss.ps'.format(mirname, mirname), shell=True, cwd='colored_structures/')
                print(''' # '{}' file generation succesfull'''.format(mirname))

                if 'mirna2' in locals():
                    del mirna2

    subprocess.run('rm *.ps', shell=True, cwd='colored_structures/')

    if extra:
        plt.clf()
        plt.boxplot(sizelst)
        plt.title('Precursor length')
        plt.ylabel('Sequence length (nt)')
        plt.savefig('length.png', dpi=300)

        plt.clf()
        plt.boxplot(mfelst)
        plt.title('Predicted minimum free energy')
        plt.ylabel('Minimum free energy (kJ/mol)')
        plt.savefig('mfe.png', dpi=300)
        
        plt.clf()
        plt.scatter(sizelst, mfelst, edgecolors='black', color='green')
        x = np.array(sizelst).reshape((-1, 1))
        y = np.array(mfelst)
        model = LinearRegression().fit(x, y)
        plt.plot(x, model.predict(x), color='black')
        plt.rcParams.update({'font.size':8})
        plt.xlabel('Precursor length')
        plt.ylabel('Minimum free energy (kJ/mol)')
        plt.savefig('mfexlength.pdf')


    os.chdir('../../')


for file in inputs:
    initial_check(file)
    folding(file)

