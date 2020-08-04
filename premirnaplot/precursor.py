
from itertools import product

import subprocess

import svgwrite as sw

from svglib.svglib import svg2rlg

from reportlab.graphics import renderPDF

import xml.etree.ElementTree as et

import statistics as st

import math

from copy import copy

from . import extra

#import extra

from . import imgparser
#import imgparser


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

		self.mfei2 = self.dg / (self.n_stems if self.n_stems else 1)

		# Normalized base pair propensity (dP)
  
		self.dp = self.bp_number / self.seqlen
  
		# Normalized Shannon entropy (dQ)

		#! ????
  
		# Normalized base-pair distance (dD)
  
		self.dd = self.diversity / self.seqlen
  
		# The second (Fielder) eigenvalue - degree of compactness (dF)

		#! ????
  
  		# The normalized values of all those before (zG, zP, zQ, zD, zF)

		#! ????

		# Minimum free energy index 3 (MFEI3)
  
		self.mfei3 = self.dg / self.loop_number
  
		# Minimum free energy index 4 (MFEI4)
  
		self.mfei4 = self.mfe / self.bp_number
  
		# Normalized ensemble free energy (NEFE)
  
		self.nefe = self.efe / self.seqlen
  
		# Difference (from microPred)
  
		self.diff = abs(self.mfe - self.efe) / self.seqlen


	@classmethod

	def from_file(cls, filename):

		prelist = []

		with open(filename) as arc:

			for index, line in enumerate(arc.readlines()):

				info = line[:-1].upper().replace(' ', '').replace('N', '').replace('T', 'U').split('\t')

				if len(info) < 2 or len(info) > 4:

					raise Exception("There was an error parsing your file {} at line {}".format(filename, index+1))

				annot = set(info[0]) - set('ACUG')

				if annot:

					annotation, precursor, *mirnas = info
					mirna1, mirna2 = mirnas +  (2 - len(mirnas)) * [None]
				
				else:

					precursor, *mirnas = info
					mirna1, mirna2 = mirnas +  (2 - len(mirnas)) * [None]

				mirname = annotation.lower() if annot else 'precursor_{}'.format(index)

				
				prelist.append(cls(mirname, precursor, mirna1, mirna2))
				

		return prelist

	@property

	def bp_number(self):

		''' The total number of base pairs '''

		return (self.secondary.count('(') + self.secondary.count(')')) / 2.0
	
	@property

	def n_stems(self):

		''' The number of stems in the secondary stucture. A stem is considered a motif with more than 3 consecutive base pairings '''

		return len(self.stem_positions())


	def pairs(self):

		''' Returns a list of tuples containing the indexes (positions) of paired bases '''

		buff, pairs = [], []

		for idx, symb in enumerate(self.secondary):

			if symb == '(':

				buff.append(idx)

			elif symb == ')':

				pairs.append((buff.pop(), idx))

			else:	# The base in unpaired so it's not registered or there's a unexpected symbol

				pass

		if self.bp_number != len(pairs):

			raise Exception('There was an error calculating the number of base pairings')

		return sorted(pairs)


	def bp_counts(self):

		''' Counts the number of each type of base pair '''

		freqs = {'AU':0,'GC':0,'GU':0}

		for init, end in self.pairs():

			try:
				
				freqs[self.sequence[init] + self.sequence[end]]+=1

			except:

				try:

					freqs[self.sequence[end] + self.sequence[init]]+=1

				except:

					raise Exception('Could not account for this base pair: {}'.format(self.sequence[init] + self.sequence[end]))

		return freqs


	def bp_props(self):

		''' The proportion of each type of base pairing divide by the total number of base pairs '''

		freqs = self.bp_counts()

		nfreqs = {}

		for base_pair in freqs:

			nfreqs[base_pair + '\\bp_number'] = freqs[base_pair] / self.bp_number

		return nfreqs

	
	def bp_per_stems(self):

		''' The proportion of each type of base pairing divide by the total number of stems '''

		freqs = self.bp_counts()

		nfreqs = {}

		for base_pair in freqs:

			nfreqs[base_pair + '\\n_stems'] = freqs[base_pair] / self.n_stems

		return nfreqs


	def stem_positions(self):

		''' Returns a list of 4 indexes (positions of the stem): beggining and end of strand 1 and 2 '''

		buff = []
		cons = False
		stems = []
		pos = self.pairs()
		pos.append((-1,-1))
		# Estou chegando o pareamento de só uma das fitas, assumindo que estão perfeitamente pareadas
		
		while len(pos) > 1:
			
			i1, e2 = pos.pop(0)
			i2, e1 = pos[0]
			
			if i1 + 1 == i2 and e2 - 1 == e1:
				
				if not cons:
					buff.append((i1,e2))
				
				cons = True
			
			else:
				
				if cons:
					
					buff.append((i1,e2))
					cons = False
					
					if buff[1][0] - buff[0][0] > 2:
						
						(i1, e2),(e1, i2) = buff
						stems.append((i1, e1, i2, e2))
				
				buff = []

		return stems


	def avg_bp_stems(self):

		''' The average base pair proportion in the stems '''

		n = self.n_stems

		freqs = {'AU':0,'GC':0,'GU':0}

		for init1, end1, init2, end2 in self.stem_positions():

			print(init1, end1, init2, end2)
			print(self.sequence[init1:end1+1], self.sequence[init2:end2+1][::-1])

			for nt1, nt2 in zip(self.sequence[init1:end1+1], self.sequence[init2:end2+1][::-1]):

				print(nt1, nt2)

				try:
					
					freqs[nt1 + nt2]+=1

				except KeyError:

					try:

						freqs[nt2 + nt1]+=1

					except KeyError:

						raise Exception('Could not account for this base pair: {}'.format(nt1 + nt2))


		return {'avg_' + base_pair + r'\n_stems':freqs[base_pair] / n for base_pair in freqs}


	def longest_stem(self):

		''' The size of the longest stem '''

		return max([
			end1 - init1 for init1, end1, *_ in self.stem_positions()
		])


	def nt_props(self):

		''' The proportion of each nucleotide in the precursor sequence '''

		props = {'A':0,'C':0,'U':0,'G':0}

		for nt in self.sequence:

			props[nt]+=1

		return {'%' + nt : round(props[nt] / self.seqlen * 100, 2) for nt in props}


	@property

	def avg_stem_length(self):

		return sum([end - init for init, end, *_ in self.stem_positions()]) / len(self.n_stems)


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

			return round(100 * (self.mfe - extra.refmfe[self.seqlen]) / (self.seqlen - extra.SHIFT_CONST), 2)

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


	def dint_props(self):

		''' The dinucleotide frequencies in the precursor sequence '''

		return {'%' + nt1 + nt2 : round(self.sequence.count(nt1 + nt2) / (self.seqlen - 1) * 100, 2) for nt1, nt2 in product(['A', 'U', 'C', 'G'], repeat=2)}

	@property

	def bulge_number(self):

		''' Determines the number of bulges inside the secondary structure '''

		return len(self.bulges_pos())

	
	def bulges_pos(self):

		bulges = []

		stems = self.stem_positions()

		for n in range(len(stems) - 1):

			# Getting the positions from the first stem
			fi1, fe1, fi2, fe2 = stems[n]
			
			# Getting the positions from the second stem
			si1, se1, si2, se2 = stems[n+1]

			# Checking if the first strand is consecutive
			strand1 = fe1 + 1 == si1

			# Checking if the second strand is consecutive
			strand2 = se2 + 1 == fi2

			# If the first strand is consec. and the second is not, there is a bulge in the second one
			if strand1 and not strand2:

				bulges.append((se2, fi2))

			# If the second strand is consec. and the first one is not, there is a bulge in the second one
			elif not strand1 and strand2:

				bulges.append((fe1, si1))

			# If both conditions are False, that means it's a loop

			elif not strand1 and not strand2:

				pass

			# If both conditions are True, that means the two stems are actually connected and should have been joined

			else:

				raise Exception('There is an incompatibility between the stem joining and stem number')

		return bulges

	@property

	def loop_number(self):

		''' Determines the number of loops inside the secondary structure '''

		return len(self.loops_pos())

	
	def loops_pos(self):

		loops = []

		stems = self.stem_positions()

		for n in range(len(stems) - 1):

			# Getting the positions from the first stem
			fi1, fe1, fi2, fe2 = stems[n]
			
			# Getting the positions from the second stem
			si1, se1, si2, se2 = stems[n+1]

			# Checking if the first strand is consecutive
			strand1 = fe1 + 1 == si1

			# Checking if the second strand is consecutive
			strand2 = se2 + 1 == fi2

			# If both conditions are False, that means it's a loop

			if strand1 and strand2:

				raise Exception('There is an incompatibility between the stem joining and stem number')

			elif not strand1 and not strand2:

				loops.append((fe1, si1, se2, fi2))

			else:

				pass

		return loops


	def features(self):

		features = copy(vars(self))
		
		for nonfeature in ['name', 'sequence', 'mirna1', 'mirna2', 'pos1', 'pos2', 'secondary', 'pseudo', 'centroid']:

			features.pop(nonfeature)

		features.update(self.nt_props())

		features.update(self.dint_props())

		features.update(self.triplets())

		#if self.n_stems:
		
		features.update(self.bp_props())

		features.update(self.bp_per_stems())

		features.update(self.avg_bp_stems())


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


	def __rnaplot_svg_parser(self, filepath):

		sequence = None
		locations = []
		pairs = []
		transform = None
		seqtransform = None
		radius = 0

		tree = et.parse(filepath)
		root = tree.getroot()

		for child in root:
			
			self.transform = child.attrib.get('transform', None)

		for child in root:
			
			if child.tag == '{http://www.w3.org/2000/svg}g':
				
				container = child
				break

		if container is None:
			
			raise Exception("Could not parse SVG: Cannot find container")

		box = (float(root.attrib['width']), float(root.attrib['height']))
		
		locations = __locations__(container)
		
		pairs = __pairs__(container)

		if not pairs:
			
			raise Exception("Did not find any pairs")

		if not locations:
			
			raise Exception("Did not find drawing coordinates (locations)")

		return box, locations, pairs


	def createSVG(self, style=3, color1='red', color2='green', folder='./', pdf=False):

		self.__rna_plot(folder=folder)

		imgparser.SVGconstructor(folder + self.name + '_ss.svg', style, self.pos1, self.pos2, color1, color2, pdf=pdf)


prec = Precursor('let7-1', 'UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA', 'UGGGAUGAGGUAGUAGGUUGU', 'AUCUACUGUCUUUCCUA')

prec.stem_positions()