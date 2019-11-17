
# Importing necessary packages

import subprocess

import argparse

import matplotlib.pyplot as plt

import os

from src.help import *

import concurrent.futures as cf

from src.imgparser import SVGconstructor as constructor

from src.help import Precursor

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

outdir = args.outdir

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

print(args.style)

filedata = {}

# Functions used in the program


def initial_check(filename):
    
    name = filename.split('/')[-1][:-4]
    
    with open(filename) as arc:

        prelist = []

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
                    mirna2 = None
            else:
                precursor = line[0]
                mirna1 = line[1]
                if len(line) > 2:
                    mirna2 = line[2]

            mirname = annotation if annot else 'precursor_{}'.format(index)

            prelist.append(Precursor(mirname, mirna1, mirna2, precursor))

    filedata[filename] = prelist


def folding(prec):

    rnafold = subprocess.Popen(f'RNAfold > {prec.name}_fold.txt',  stdin=subprocess.PIPE,
                                            shell=True,
                                            universal_newlines=True,
                                            cwd='foldings/')

    rnafold.communicate('>{}\n{}'.format(prec.name, prec.premirna))

    # Running RNAplot to generate the initial SVG in the colored_structures/ folder

    with open('foldings/{}_fold.txt'.format(prec.name)) as dot:
        stuff, mfe = dot.read().split(' ')
        prec.premfe = float(mfe[1:-2])

    rnaplot = subprocess.Popen(['RNAplot -o svg --filename-full'], stdin=subprocess.PIPE, cwd='colored_structures/', shell=True, universal_newlines=True)
    rnaplot.communicate(stuff)

    constructor('colored_structures/' + prec.name + '_ss.svg', args.style[0], prec.pos1, prec.pos2, color1, color2, 'grey', pdf=pdf)

    return f'# Created {prec.name} image'


subprocess.run(f'mkdir {outdir}/', shell=True)


for file in inputs:
    
    filedata[file] = []

    print(f"#######  Checking if {file} is ok..\n")
    
    initial_check(file)

    print('\n#######  Data check complete for {}'.format(file))


for file in filedata:
    
    mfelst = []
    sizelst = []
    mirdict = {}
    name = (file.split('/')[-1][:-4] if '.' in file else name)
    
    subprocess.run(f'mkdir {outdir}/{name} {outdir}/{name}/foldings {outdir}/{name}/colored_structures', shell=True)
    
    os.chdir(f'{outdir}/{name}/')

    results = []

    with cf.ProcessPoolExecutor(max_workers=nthreads) as executor:
        
        results = executor.map(folding, filedata[file])

        for result in results:
            print(result)
    

    os.chdir('../../')

