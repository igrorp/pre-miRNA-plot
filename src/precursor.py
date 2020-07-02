
from itertools import product

import subprocess

from help import refmfe, SHIFT_CONST

from svglib.svglib import svg2rlg

from reportlab.graphics import renderPDF

import xml.etree.ElementTree as et

import svgwrite as sw

import statistics as st

import math

from imgparser import SVGconstructor


class Precursor():

	def __init__(self, name, precursor, mirna1, mirna2):

		self.name = name

		# The precursor nucleotide sequence

		#todo Maybe add some validation to the initial given parameters, i.e. make sure that it's nucleotides and it's not empty

		self.sequence = precursor.replace('T', 'U')

		if mirna1 and mirna2:

			mirna1 = mirna1.replace('T', 'U')
			mirna2 = mirna2.replace('T', 'U')

			# self.mirnas = (mirna1 if mirna1 else '', mirna2 if mirna2 else '')

			posa, posb = self.__pos(mirna1)

			posc, posd = self.__pos(mirna2)

			if posa < posc:
				self.mirna1 = mirna1
				self.mirna2 = mirna2
				self.pos1 = (posa, posb)
				self.pos2 = (posc, posd)
			else:
				self.mirna1 = mirna2
				self.mirna2 = mirna1
				self.pos1 = (posc, posd)
				self.pos2 = (posa, posb)

		elif mirna1 and not mirna2:

			self.mirna1 = mirna1
			self.mirna2 = ''

		elif mirna2 and not mirna1:

			self.mirna1 = ''
			self.mirna2 = mirna2

		# The length of the pre-miRNA

		self.seqlen = len(self.sequence)

		self.__rnafold()

		# The GC content of the precursor

		self.gccontent = self.gc(self.sequence)

		self.mirna1gc = self.gc(mirna1) if self.mirna1 else 'N.A.'

		self.mirna2gc = self.gc(mirna2) if self.mirna2 else 'N.A.'

		# The MFEdensity of the precursor

		self.mfeden

		# The number of mismatches in the region of the miRNA duplex

		self.mismatches()
  
		#todo Calculate number of stems and number of loops


		####* Calculating the features
  
		# Normalized minimum free energy of folding (dG)
			
		self.dg = self.mfe / self.seqlen
  
		# Minimum free energy index 1 (MFEI1)

		self.mfei1 = self.dg / self.gccontent
  
		# Minimum free energy index 2 (MFEI2)

		#! self.mfei2 = self.dg / self.n_stems
  
		# Normalized base pair propensity (dP)
  
		#! self.dp = self.tot_bases / self.seqlen
  
		# Normalized Shannon entropy (dQ)

		#! ????
  
		# Normalized base-pair distance (dD)
  
		self.dd = self.diversity / self.seqlen
  
		# The second (Fielder) eigenvalue - degree of compactness (dF)

		#! ????
  
  		# The normalized values of all those before (zG, zP, zQ, zD, zF)

		#! ????

		# Minimum free energy index 3 (MFEI3)
  
		#! self.mfei3 = self.dg / self.n_loops
  
		# Minimum free energy index 4 (MFEI4)
  
		#! self.mfei4 = self.mfe / self.tot_bases
  
		# Normalized ensemble free energy (NEFE)
  
		self.nefe = self.efe / self.seqlen
  
		# Difference (from microPred)
  
		self.diff = abs(self.mfe - self.efe) / self.seqlen


	@classmethod

	def from_file(cls, filename):

		prelist = []

		with open(filename) as arc:

			for index, line in enumerate(arc.readlines()):

				info = line[:-1].upper().replace(' ', '').replace('N', '').split('\t')

				if len(info) < 2 or len(info) > 4:

					raise Exception("There was an error parsing your file {} at line {}".format(filename, index+1))

				annot = set(info[1]) - set('ACTG')

				if annot:

					annotation, precursor, *mirnas = info
					mirna1, mirna2 = mirnas +  (2 - len(mirnas)) * [None]
				
				else:

					precursor, *mirnas = info
					mirna1, mirna2 = mirnas +  (2 - len(mirnas)) * [None]

				mirname = annotation if annot else 'precursor_{}'.format(index)

				prelist.append(cls(mirname, precursor, mirna1, mirna2))

		return prelist

	@property

	def tot_bases(self):

		return self.sequence.count('(') + self.sequence.count(')')
	

	def __rnafold(self, folder='./'):

		''' Runs the RNAfold secondary structure prediction with the partition function and
		pairing probability matrix calculation. Also, parses the output and can direct the
		generated PS image file to a specific path. '''

		rnafold = subprocess.Popen(['RNAfold -p --noPS --noDP'],
										stdin=subprocess.PIPE,
										stdout=subprocess.PIPE,
										stderr=subprocess.PIPE,
										shell=True,
										universal_newlines=True,
										cwd=folder)

		data, errormsg = rnafold.communicate('>{}\n{}'.format(self.name, self.sequence))

		if rnafold.returncode:

			raise Exception('There was an error running RNAfold for precursor {} with the following message:\n{}'.format(self.name, errormsg))

		elif data == '':

			raise Exception('There was an error running RNAfold for precursor {}: the output was empty')

		else:

			self.__rnafold_parser(data)


	def __rnafold_parser(self, data):

		''' Parses the output from RNAfold -p and sets the Precursor properties'''

		try:

			*_, mfedata, pseudo, centroid, ensemble, _ = data.split('\n')

		except ValueError:

			raise Exception('Could not parse the RNAfold output')

		#print(centroid)

		# The predicted secondary strucuture with the lowest energy in dot-bracket notation and its MFE value

		secondary, *mfe = mfedata.split(' ')

		self.secondary = secondary

		self.mfe = float(''.join(mfe)[1:-1])

		# Pseudo bracket notation of pair probabilities and ensemble free energy (EFE)

		pseudo, *efe = pseudo.split(' ')

		self.pseudo = pseudo

		self.efe = float(''.join(efe)[1:-1])

		# Centroid ensemble structure dot bracket notation, its free energy and distance from the ensemble

		*notation, energy, dist = centroid.split(' ')

		self.centroid = notation

		self.ctdenergy = float(energy.replace(' ', '')[1:])

		self.ctddist = float(dist.replace(' ', '')[2:-1])

		# The frequency of the MFE structure and the ensemble strucutral diversity (mean base pair distance)

		separated = ensemble.split(';')

		self.freq = float(separated[0].split(' ')[-1])

		self.diversity = float(separated[1].split(' ')[3])


	def __pos(self, mirna):

		if mirna:

			if mirna in self.sequence:

				if self.sequence.count(mirna) > 1:

					print('WARNING! miRNA {} was found more than once in the precursor sequence {}..., but its last occurrence will be used!'.format(mirna, self.sequence[25:]))

				return (self.sequence.find(mirna), self.sequence.find(mirna) + len(mirna))

			else:

				raise Exception('ERROR! Could not find sequence {} inside {}, please correct this'.format(mirna, self.sequence))

		else:

			return None, None

	@property

	def mfeden(self):

		if self.seqlen >= 40 and self.seqlen <= 600:

			return round(100 * (self.mfe - refmfe[self.seqlen]) / (self.seqlen - SHIFT_CONST), 2)

		else:

			return 'N.A.'

	@staticmethod

	def gc(seq):

		return round((seq.count('G') + seq.count('C')) / len(seq) , 2)


	def mismatches(self):

		if self.mirna1:

			posa, posb = self.pos1

			self.mirna1mm = self.secondary[posa:posb].count('.')

		else:

			self.mirna1mm = 'N.A.'

		if self.mirna2:

			posc, posd = self.pos2

			self.mirna2mm = self.secondary[posc:posd].count('.')

		else:

			self.mirna2mm = 'N.A.'

		if self.mirna1 and self.mirna2:

			self.duplexmm =  self.mirna1mm + self.mirna2mm

		else:

			self.duplexmm = 'N.A.'


	def setpredsec(self, secstruct):

		if self.seqlen != len(secstruct):

			raise Exception('Predicted secondary structure and precursor sequence have different lengths')

		else:

			self.secondary = secstruct


	def loop(self):

		lastopen, firstclose = 0, len(self.sequence)

		for idx, nt in enumerate(self.secondary):

			if nt == '(':

				lastopen = idx

			elif nt == ')':

				firstclose = idx

				break

		return lastopen, firstclose


	def triplets(self):

		triplets = {
			'A.((':0, 'A(..':0, 'A..(':0, 'A((.':0, 'A(((':0, 'A...':0, 'A(.(':0, 'A.(.':0,
			'C.((':0, 'C(..':0, 'C..(':0, 'C((.':0, 'C(((':0, 'C...':0, 'C(.(':0, 'C.(.':0,
			'U.((':0, 'U(..':0, 'U..(':0, 'U((.':0, 'U(((':0, 'U...':0, 'U(.(':0, 'U.(.':0,
			'G.((':0, 'G(..':0, 'G..(':0, 'G((.':0, 'G(((':0, 'G...':0, 'G(.(':0, 'G.(.':0,
		}

		if self.secondary == '':

			raise Exception('Its necessary to have a secondary structure')

		else:

			new = self.secondary.replace(')', '(')


		loopinit, loopend = self.loop()

		# fivepstem = self.secondary[self.secondary.find('(') : loopinit + 1]
		# threepstem = self.secondary[loopend : self.secondary.rfind(')') + 1]

		for n in range(self.secondary.find('(') + 1, loopinit):

			triplets[self.sequence[n] + new[n-1:n+2]]+=1

		for n in range(loopend, self.secondary.rfind(')') - 1):

			triplets[self.sequence[n] + new[n-1:n+2]]+=1

		soma = sum(triplets.values())

		triplets = {triplet : round(freq / soma, 5) for triplet, freq in triplets.items()}

		return triplets


	def porcents(self):

		return {nt1 + nt2 : round(self.sequence.count(nt1 + nt2) / (self.seqlen - 1) * 100, 2) for nt1, nt2 in product(['A', 'T', 'C', 'G'], repeat=2)}


	def features(self):

		features = vars(self)

		features.update(self.porcents())

		features.update(self.triplets())

		return features

	
	#? Image definitions

	def __rna_plot(self, folder='./'):

		''' Creates the RNAplot SVG image '''

		rnaplot = subprocess.Popen(['RNAplot -o svg'],
										stdin=subprocess.PIPE,
										stdout=subprocess.PIPE,
										stderr=subprocess.PIPE,
										shell=True,
										universal_newlines=True,
										cwd=folder)

		data, errormsg = rnaplot.communicate('>{}\n{}\n{}'.format(self.name, self.sequence, self.secondary))

		if rnaplot.returncode:

			raise Exception('There was an error running RNAplot for precursor {} with the following message:\n{}'.format(self.name, errormsg))


	def createSVG(self, style, color1, color2, folder='./', pdf=False):

		self.__rna_plot(folder=folder)

		SVGconstructor(folder + self.name + '_ss.svg', style, self.pos1, self.pos2, color1, color2, pdf=pdf)


prec = Precursor('let7-1', 'UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA', 'UGGGAUGAGGUAGUAGGUUGU', 'AUCUACUGUCUUUCCUA')

