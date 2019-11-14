
# Importing necessary packages

import subprocess

import argparse

import matplotlib.pyplot as plt

import os

from src.help import *

import concurrent.futures as cf

from src.imgparser import SVGconstructor as constructor 

# Adding parameters and a general help message

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=desctxt, add_help=False)

parser.add_argument('input', nargs='*')

parser.add_argument('-h', '--help', default=argparse.SUPPRESS)

parser.add_argument('-i', '--input', nargs='+', metavar='', dest='inputopt')

parser.add_argument('-a', '--annotation', default='F', choices=['T', 'F'])

parser.add_argument('-e', '--extra_info', default='T', choices=['T', 'F'])

parser.add_argument('-c', '--colors', nargs='+', default=['red', 'green'], type=str)

parser.add_argument('-t', '--threads', nargs=1, default=1, type=int)

parser.add_argument('-s', '--style', nargs=1, default=1, type=int)

parser.add_argument('-f', '--outfmt', nargs=1, default='svg', choices=['pdf', 'svg'], type=str)

parser.add_argument('-o', '--outdir', nargs=1, default='premirnaplot', type=str)

args = parser.parse_args()


# Parsing the arguments

inputs = args.inputopt if args.input == [] else args.input

if not inputs:
    print(desctxt)
    quit()

annot = True if args.annotation == 'T' else False

extra = True if args.extra_info == 'T' else False

nthreads = args.threads[0]

pdf = True if args.outfmt[0] == 'pdf' else None

print(pdf)
print(args.outfmt)

outdir = args.outdir[0]

if len(args.colors) == 2:

    color1, color2 = args.colors

    if color1 not in defcolors or color2 not in defcolors:
        
        raise Exception("ERROR! One of the colors you informed is incorrect, please check your spelling or review the predefined colors.")

    else:
        color1, color2 = defcolors[color1], defcolors[color2]

elif len(args.colors) == 6:

    for color in args.colors:
        if color < 0 or color > 255:
            raise Exception("ERROR! Please use RGB code values between 0 and 255")
    
    color1, color2 = ' '.join([colors[0], colors[1], colors[2]]), ' '.join([colors[3], colors[4], colors[5]])

else:
    print("\nThere was an error checking the colors you provided, please review them")


# Functions used in the program


def initial_check(filename):
    
    name = filename.split('/')[-1][:-4]
    
    with open(filename) as arc:
        for index, line in enumerate(arc.readlines()):
                
                line = line[:-1].upper().replace(' ', '').split('\t')

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

                isin(precursor, mirna1, index, name)

                repetition(precursor, mirna1, index)

                if 'mirna2' in locals():
                    isin(precursor, mirna2, index, name)
                    repetition(precursor, mirna2, index)
                    del mirna2

        print('\n','########  Data check complete for {}'.format(filename.split('/')[-1]), '\n')


def isin(premirna, mirna, i, filename):
    if mirna in premirna:
        pass
    else:
        print('''   ERROR: Sequence {} at line {} in file '{}' could not be found in {}..., please correct this!'''.format(mirna, str(i+1), filename, premirna[:20]))
        valid = 0


def repetition(premirna, mirna, i):
    if premirna.count(mirna) > 1:
        print('   Warning! miRNA sequence {} at line {} was found more than once in the precursor sequence, so its last occurrence will be used for the plotting, be aware!'.format(mirna, str(1+i)))


def pos(premirna, mirna):

    if mirna in premirna:
        return (premirna.find(mirna) + 1, premirna.find(mirna) + len(mirna))
    else:
        return None


valid = 1

subprocess.run(f'mkdir {outdir}/', shell=True)

for file in inputs:
    
    print(f"#######  Checking if {file} is ok..\n")
    initial_check(file)

if not valid:
    print("Plase correct the errors and run the program again.")
    quit()


for file in inputs:
    
    mfelst = []
    sizelst = []
    mirdict = {}
    name = file.split('/')[-1][:-4]
    
    subprocess.run(f'mkdir {outdir}/{name} {outdir}/{name}/foldings {outdir}/{name}/colored_structures', shell=True)
    
    with open(file) as arc:

        os.chdir(f'{outdir}/{name}/')

        with open('precursor_data.txt', 'a') as data:
            data.write('miRNA name\tSequence\tMFE\tLenght\n')

        for index, line in enumerate(arc.readlines()):
                
                line = line[:-1].upper().replace(' ', '').split('\t')
                
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

                mirname = annotation if annot else 'precursor_{}'.format(i)

                # Determining the positions of the miRNAs inside the precursor

                pos1 = pos(precursor, mirna1)

                if 'mirna2' in locals():
                    pos2 = pos(precursor, mirna2)

                # Running RNAfold to predict the data in the foldings/ folder

                rnafoldout = open(f'foldings/{mirname}_fold.txt', 'w')
                
                rnafold = subprocess.Popen('RNAfold ',  stdin=subprocess.PIPE,
                                                        stdout=rnafoldout,
                                                        shell=True,
                                                        universal_newlines=True,
                                                        cwd='foldings/')

                rnafold.communicate('>{}\n{}'.format(mirname, precursor))

                # Running RNAplot to generate the initial SVG in the colored_structures/ folder

                rnaplot = subprocess.Popen(['RNAplot -o svg --filename-full'], stdin=subprocess.PIPE, cwd='colored_structures/', shell=True, universal_newlines=True)
                rnaplot.communicate(open('foldings/{}_fold.txt'.format(mirname)).read().split(' ')[0])

                mirdict[mirname] = (pos1, pos2)
                
                if extra:
                    mfe = open('foldings/{}_fold.txt'.format(mirname)).read().split(' ')[-1][1:-2]
                    mfelst.append(float(mfe))
                    sizelst.append(float(len(precursor)))
                    with open('precursor_data.txt', 'a') as data:
                        data.writelines([mirname, '\t', precursor, '\t', mfe, '\t', str(len(precursor)), '\n'])


                if 'mirna2' in locals():
                    del mirna2

    
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
        plt.rcParams.update({'font.size':8})
        plt.xlabel('Precursor length')
        plt.ylabel('Minimum free energy (kJ/mol)')
        plt.savefig('mfexlength.pdf')

    
    with cf.ThreadPoolExecutor(max_workers=nthreads) as executor:
        
        for key in mirdict:
            pos1, pos2 = mirdict[key]
            executor.submit(constructor('colored_structures/' + key + '_ss.svg', args.style[0], pos1, pos2, color1, color2, 'grey', pdf=pdf))
            # constructor('colored_structures/' + key + '_ss.svg', args.style[0], pos1, pos2, color1, color2, 'grey', pdf=pdf)
            # print('Image for', key, 'created')


    os.chdir('../../')

