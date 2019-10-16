
# Importing necessary packages

import subprocess

import argparse

import matplotlib.pyplot as plt

import os

from help import *

import concurrent.futures as cf


# Adding parameters and a general help message


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=desctxt, add_help=False)

parser.add_argument('input', nargs='*', metavar='(INPUT)')

parser.add_argument('-h', '--help', default=argparse.SUPPRESS)

parser.add_argument('-i', '--input', nargs='+', metavar='', dest='inputopt')

parser.add_argument('-a', '--annotation', default='F', choices=['T', 'F'], metavar='\b')

parser.add_argument('-e', '--extra_info', default='T', choices=['T', 'F'], metavar='\b')

parser.add_argument('-c', '--colors', nargs='+', default=[defcolors['red'] + defcolors['green']], type=str)

parser.add_argument('-t', '--threads', nargs=1, default=1, type=int)

parser.add_argument('-q', '--quality', nargs=1, default=200, type=int)

args = parser.parse_args()

annot = True if args.annotation == 'T' else False

extra = True if args.extra_info == 'T' else False

inputs = args.inputopt if args.input == [] else args.input

colors = args.colors

nthreads = args.threads[0]

qty = args.quality[0]

if not inputs:
    print(desctxt)
    quit()

if len(colors) == 2:

    fc, tc = args.colors[0], args.colors[1]

    if fc not in defcolors or tc not in defcolors:
        print("\nOne of the colors you informed is not predefined, please check your spelling or use the RGB code")
        quit()
    else:
        fc, tc = ' '.join(defcolors[fc]), ' '.join(defcolors[tc])

elif len(colors) == 6:
    
    for i in range(len(colors)):
        colors[i] = str(int(colors[i]) / 255)

    fc, tc = ' '.join([colors[0], colors[1], colors[2]]), ' '.join([colors[3], colors[4], colors[5]])
else:
    print("\nThere was an error checking the colors you provided, please review them")


# Functions used in the program


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
        print('''   ERROR: Sequence {} at line {} in file '{}' could not be found in {}..., please correct this!'''.format(mirna, str(i+1), filename, premirna[:20]))
        valid = 0


def repetition(premirna, mirna, i):
    if premirna.count(mirna) > 1:
        print('   Warning! miRNA sequence {} at line {} was found more than once in the precursor sequence, so its last occurrence will be used for the plotting, be aware!'.format(mirna, str(1+i)))


def pos(premirna, mirna1, mirna2):
    if mirna2 == '':
        return [premirna.find(mirna1) + 1, premirna.find(mirna1) + len(mirna1)]
    elif premirna.count(mirna1) == 1:
        return [premirna.find(mirna1) + 1, premirna.find(mirna1) + len(mirna1), premirna.rfind(mirna2) + 1, premirna.rfind(mirna2) + len(mirna2)]
    elif premirna.count(mirna2) == 1:
        return [premirna.find(mirna2) + 1, premirna.find(mirna2) + len(mirna2), premirna.rfind(mirna1) + 1, premirna.rfind(mirna1) + len(mirna1)]


def converter(mirname):
    subprocess.run(f'gs -q -dSAFER -dBATCH -dNOPAUSE -sPAPERSIZE=a4 -r{qty} -sDEVICE=pngalpha -sOutputFile={mirname}.png {mirname}_ss.ps', shell=True, cwd="colored_structures/")
    subprocess.run('rm {}'.format(mirname), shell=True, cwd="colored_structures/")
    
    return f"Image for {mirname} createad succesfully"


# Starting the program

valid = 1

subprocess.run('mkdir premirnaplot/', shell=True)

for file in inputs:
    
    print(f"#######  Checking if {file} is ok..\n")
    initial_check(file)

if not valid:
    print("Plase correct the errors and run the program again.")
    quit()

for file in inputs:
    
    mfelst = []
    sizelst = []
    mirlist = []
    name = file.split('/')[-1][:-4]
    
    subprocess.run(f'mkdir premirnaplot/{name} premirnaplot/{name}/foldings premirnaplot/{name}/colored_structures', shell=True)
    
    with open(file) as arc:

        os.chdir(f'premirnaplot/{name}/')

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

                mirname = annotation if annot else 'precursor_{}'.format(i)

                rnafoldout = open(f'foldings/{mirname}_fold.txt', 'w')
                
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

                mirlist.append(mirname)

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
        results = executor.map(converter, mirlist)

        for result in results:
            print(result)


    os.chdir('../../')

