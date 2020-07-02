#!/usr/bin/env python

# Importing necessary packages

import pandas as pd

import numpy as np

import subprocess

import argparse

import matplotlib.pyplot as plt

import os

from help import *

import concurrent.futures as cf

from sklearn.linear_model import LinearRegression

from imgparser import SVGconstructor as constructor

from precursor import Precursor

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





#! Criar arquivo _fold.txt ???


def folding(prec):

	# with open('foldings/{}_fold.txt'.format(prec.name), 'w') as dot:

	#constructor('colored_structures/' + prec.name + '_ss.svg', args.style, prec.pos1, prec.pos2, color1, color2, pdf=pdf)

	return prec


# Creating the output directory

subprocess.run('rm -r {}'.format(outdir), shell=True)

sp = subprocess.run('mkdir {}/'.format(outdir), shell=True, capture_output=True, universal_newlines=True)

if sp.stderr:

	raise Exception('Could not create the output directory')


# Creating the folder for each file
# Instantiating each one of the precursors from the given files

precs = []

for file in inputs:

	prec += [Precursor.from_file(file)]



for file in filedata:

	mfelst = []
	sizelst = []

	name = (file.split('/')[-1][:-4] if '.' in file else name)

	subprocess.run('mkdir {}/{}'.format(outdir, name), shell=True)

	os.chdir('{}/{}/'.format(outdir, name))

	#data = pd.DataFrame({'Names':[]})

	data = []

	with cf.ThreadPoolExecutor(max_workers=nthreads) as executor:

		for idx, precursor in enumerate(executor.map(folding, filedata[file])):

			#print('# Created {} image'.format(precursor.name))

			fields = vars(precursor)

			#params = {'name':'Names', 'premirna':'Precursor sequence', 'predsec':'Secondary structure', 'mirna1':'miRNA5p', 'mirna2':'miRNA3p', 'prelen':'Precursor length', 'premfe':'MFE', 'mfeden':'MFEden',
            #			'duplexmm':'Duplex MM', 'mirna1mm':'miRNA5p mm', 'mirna2mm':'miRNA3p mm', 'gccontent':'''%GC''',  'mirna1gc':'''%GC miRNA5p''', 'mirna2gc':'''%GC miRNA3p''',}
			print(len(precursor.features()))
			#print(precursor.features())

			data.append(list(precursor.features().values()))

	#data = data.set_index('Names')

	data = pd.DataFrame(data, columns=list(precursor.features().keys()))
	
	data.to_csv(path_or_buf='precursor_data.txt', sep='\t')

	plt.clf()
	plt.boxplot(data['seqlen'])
	plt.title('Precursor length')
	plt.ylabel('Sequence length (nt)')
	plt.savefig('length.png', dpi=500)

	plt.clf()
	plt.boxplot(data['mfe'])
	plt.title('Predicted minimum free energy')
	plt.ylabel('Minimum free energy (kJ/mol)')
	plt.savefig('mfe.png', dpi=500)

	plt.clf()
	plt.scatter(data['seqlen'], data['mfe'], edgecolors='black', color=color1, zorder=2)
	x = data['seqlen'].values.reshape((-1, 1))
	y = data['mfe']
	model = LinearRegression().fit(x, y)
	plt.plot(x, model.predict(x), color='black', zorder=1)
	plt.rcParams.update({'font.size':8})
	plt.xlabel('Precursor length')
	plt.ylabel('Minimum free energy (kJ/mol)')
	plt.savefig('mfexlength.png', dpi=500)

	os.chdir('../../')

