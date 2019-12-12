
# Importing necessary packages

import pandas as pd

import subprocess

import argparse

import matplotlib.pyplot as plt

import os

from src.help import *

import concurrent.futures as cf

from sklearn.linear_model import LinearRegression

from src.imgparser import SVGconstructor as constructor

from src.help import Precursor

# Adding parameters and a general help message

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=desctxt, add_help=False)

parser.add_argument('input', nargs='*')

parser.add_argument('-h', '--help', default=argparse.SUPPRESS)

parser.add_argument('-i', '--input', nargs='+', metavar='', dest='inputopt')

parser.add_argument('-a', '--annotation', default='F', choices=['T', 'F'])

parser.add_argument('-c', '--colors', nargs='+', default=['red', 'green'], type=str)

parser.add_argument('-t', '--threads', default=1, type=int)

parser.add_argument('-s', '--style', default=3, type=int)

parser.add_argument('-f', '--outfmt', default='svg', choices=['pdf', 'svg'], type=str)

parser.add_argument('-o', '--outdir', default='premirnaplot', type=str)

args = parser.parse_args()


# Parsing the arguments

inputs = args.inputopt if args.input == [] else args.input

if not inputs:
    print(desctxt)
    quit()

annot = True if args.annotation == 'T' else False

nthreads = args.threads

pdf = True if args.outfmt == 'pdf' else None

outdir = args.outdir

if len(args.colors) == 2:

    color1, color2 = args.colors

    if color1 not in defcolors or color2 not in defcolors:
        
        raise Exception("ERROR! One of the colors you informed is incorrect, please check your spelling or review the predefined colors.")

    else:
        color1, color2 = defcolors[color1], defcolors[color2]

elif len(args.colors) == 6:

    for color in args.colors:
        if int(color) < 0 or int(color) > 255:
            raise Exception("ERROR! Please use RGB code values between 0 and 255")
    
    color1 = "#{:02x}{:02x}{:02x}".format(int(args.colors[0]), int(args.colors[1]), int(args.colors[2]))
    color2 = "#{:02x}{:02x}{:02x}".format(int(args.colors[3]), int(args.colors[4]), int(args.colors[5]))

else:
    raise Exception("\nThere was an error checking the colors you provided, please review them")


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
                mirna1 = line[2]
                if len(line) == 4:
                    mirna2 = line[3]
                else:
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

    rnafold = subprocess.Popen('RNAfold ',  stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE,
                                            shell=True,
                                            universal_newlines=True,
                                            cwd='foldings/')

    alldata, _ = rnafold.communicate('>{}\n{}'.format(prec.name, prec.premirna))

    # Running RNAplot to generate the initial SVG in the colored_structures/ folder

    with open('foldings/{}_fold.txt'.format(prec.name), 'w') as dot:
        # alldata = dot.read().split(' ')
        # stuff = alldata[0]
        dot.write(alldata)
        mfe = float(alldata.split(' ')[-1][1:-2])


    # print(prec.secpred)
    # prec.secpred = alldata
    # print(prec.secpred)

    secpred = alldata.split('\n')[2].split(' ')[0]

    posa, posb = prec.pos1

    posc, posd = prec.pos2

    mm = secpred[posa:posb].count('.') + secpred[posc:posd].count('.')

    rnaplot = subprocess.Popen(['RNAplot -o svg --filename-full'], stdin=subprocess.PIPE, cwd='colored_structures/', shell=True, universal_newlines=True)
    rnaplot.communicate(alldata)

    constructor('colored_structures/' + prec.name + '_ss.svg', args.style, prec.pos1, prec.pos2, color1, color2, pdf=pdf)

    return prec.name, mfe, mm, secpred


subprocess.run('mkdir {}/'.format(outdir), shell=True)


for file in inputs:
    
    filedata[file] = []

    print(f"#######  Checking if {file} is ok..\n")
    
    initial_check(file)

    print('\n#######  Data check complete for {}\n'.format(file))


for file in filedata:
    
    mfelst = []
    sizelst = []
    mirdict = {'Names':[], 
    'MFEs':[0 for _ in range(len(filedata[file]))], 
    'Lenghts':[], 
    'Sequences':[], 
    'Mismatches':[0 for _ in range(len(filedata[file]))],
    'Secondary structure':[0 for _ in range(len(filedata[file]))]}

    name = (file.split('/')[-1][:-4] if '.' in file else name)
    
    subprocess.run('mkdir {}/{} {}/{}/foldings {}/{}/colored_structures'.format(outdir, name, outdir, name, outdir, name,), shell=True)
    
    os.chdir('{}/{}/'.format(outdir, name))

    for prec in filedata[file]:
        mirdict['Names'].append(prec.name)
        mirdict['Lenghts'].append(prec.prelen)
        mirdict['Sequences'].append(prec.premirna)

    data = pd.DataFrame(mirdict)
    data = data.set_index('Names')

    with cf.ProcessPoolExecutor(max_workers=nthreads) as executor:
        
        for result in executor.map(folding, filedata[file]):
            name, mfe, mm, secpred = result
            print('# Created {} image'.format(name))
            data.loc[name, 'MFEs'] = mfe
            data.loc[name, 'Mismatches'] = mm
            data.loc[name, 'Secondary structure'] = secpred

    print(data)

    data.to_csv(path_or_buf='precursor_data.txt', sep='\t')
    
    plt.clf()
    plt.boxplot(data['Lenghts'])
    plt.title('Precursor length')
    plt.ylabel('Sequence length (nt)')
    plt.savefig('length.png', dpi=500)

    plt.clf()
    plt.boxplot(data['MFEs'])
    plt.title('Predicted minimum free energy')
    plt.ylabel('Minimum free energy (kJ/mol)')
    plt.savefig('mfe.png', dpi=500)
    
    plt.clf()
    plt.scatter(data['Lenghts'], data['MFEs'], edgecolors='black', color=color1, zorder=2)
    x = data['Lenghts'].values.reshape((-1, 1))
    y = data['MFEs']
    model = LinearRegression().fit(x, y)
    plt.plot(x, model.predict(x), color='black', zorder=1)
    plt.rcParams.update({'font.size':8})
    plt.xlabel('Precursor length')
    plt.ylabel('Minimum free energy (kJ/mol)')
    plt.savefig('mfexlength.png', dpi=500)

    os.chdir('../../')

