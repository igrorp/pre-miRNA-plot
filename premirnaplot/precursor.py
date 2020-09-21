
from itertools import product

from collections import defaultdict as dfd

import subprocess

import svgwrite as sw

from svglib.svglib import svg2rlg

from reportlab.graphics import renderPDF

import xml.etree.ElementTree as et

import statistics as st

import math

from copy import copy

import re

from . import extra

#import extra

from . import imgparser
#import imgparser


class Precursor():

	def __init__(self, name, precursor, mirna1, mirna2):

		self.name = name

		#todo Maybe add some validation to the initial given parameters, i.e. make sure that it's nucleotides and it's not empty

		
		######! -----> Sequence related features

		self.sequence = precursor.upper().replace('T', 'U')
		
		self.seqlen = len(self.sequence)

		# self.nt_freqs()

		# self.dint_freqs()

		self.gccontent = self.gc(self.sequence)

		self.gcratio = self.sequence.count('G') / self.sequence.count('C')

		
		######! ----> miRNA related features
		
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
	
		self.mirna1gc = self.gc(mirna1) if self.mirna1 else 'N.A.'

		self.mirna2gc = self.gc(mirna2) if self.mirna2 else 'N.A.'

		self.__rnafold()

		# The number of mismatches in the region of the miRNA duplex

		self.mismatches()

		
		#####! ----> Thermodinamics realated features
		
		# The MFEdensity of the precursor

		self.mfeden


  
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
  
		self.mfei3 = self.dg / (self.n_loops if self.n_loops else 1) 
  
		# Minimum free energy index 4 (MFEI4)
  
		self.mfei4 = self.mfe / self.bp_number
  
		# Normalized ensemble free energy (NEFE)
  
		self.nefe = self.efe / self.seqlen
  
		# Difference (from microPred)
  
		self.diff = abs(self.mfe - self.efe) / self.seqlen

		# The stem triplets

		ignoreidx = []

		stems = self.stem_positions()

		for n in range(len(stems) - 1):

			# Getting the positions from the first stem
			fi1, fe1, fi2, fe2 = stems[n]
			
			# Getting the positions from the second stem
			si1, se1, si2, se2 = stems[n+1]

			# If both conditions are False, that means it's a loop

			if si1 - fe1 > 3:

				ignoreidx.extend([n for n in range(fe1, si1+1)])

			elif fi2 - se2 > 3:

				ignoreidx.extend([n for n in range(se2, fi2+1)])

		print(ignoreidx)
		
		self.stem_triplets = self.triplets(ignore=ignoreidx)


	#! Class methods

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


	#! Sequence related features

	@staticmethod

	def gc(seq):

		return round((seq.count('G') + seq.count('C')) / len(seq) , 2)


	def nt_freqs(self):

		''' The proportion of each nucleotide in the precursor sequence '''

		props = {'A':0,'C':0,'U':0,'G':0}

		for nt in self.sequence:

			props[nt]+=1

		return {'%' + nt : round(props[nt] / self.seqlen * 100, 2) for nt in props}


	def dint_freqs(self):

		''' The dinucleotide frequencies in the precursor sequence (%AA, %AT, %AC...) '''

		return {'%' + nt1 + nt2 : round(self.sequence.count(nt1 + nt2) / (self.seqlen - 1) * 100, 2) for nt1, nt2 in product(['A', 'U', 'C', 'G'], repeat=2)}



	#! Secondary structure related features


	def triplets(self, ignore=None):

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


		init1, loopinit, loopend, end2 = self.stem_positions()[-1]	# Getting the last stem (where the terminal loop starts)


		idxs = [n for n in range(self.secondary.find('(') + 1, loopinit)] + [n for n in range(loopend, self.secondary.rfind(')') - 1)]

		# Removing unwanted index positions
		
		if ignore:

			for idx in ignore:

				try:

					idxs.remove(idx)

				except:

					pass

		for idx in idxs:

			triplets[self.sequence[idx] + new[idx-1:idx+2]]+=1

		soma = sum(triplets.values())

		for triplet, freq in triplets.items():

			triplets[triplet] = round(freq / soma, 5)

		return triplets


	#? Stems

	def stem_positions(self):

		''' Returns a list of 4 indexes (positions of the stem): beggining and end of strand 1 and 2 '''

		buff = []
		cons = False
		stems = []
		pos = self.pairs()
		pos.append((-1,-1))
		
		while len(pos) > 1:
			
			i1, e2 = pos.pop(0)
			i2, e1 = pos[0]
			
			if i1 + 1 == i2 and e2 - 1 == e1:
				
				if not cons:

					buff.append((i1, e2))
				
				cons = True
			
			else:
				
				if cons:
					
					buff.append((i1, e2))
					cons = False
					
					if buff[1][0] - buff[0][0] > 2:
						
						(i1, e2), (e1, i2) = buff
						stems.append((i1, e1, i2, e2))
				
				buff = []

		return stems

	
	@property

	def n_stems(self):

		''' The number of stems in the secondary stucture. A stem is considered a motif with more than 3 consecutive base pairings '''

		return len(self.stem_positions())


	def longest_stem(self):

		''' The size of the longest stem '''

		return max([
			end1 - init1 for init1, end1, *_ in self.stem_positions()
		])


	@property

	def avg_stem_length(self):

		return sum([end - init for init, end, *_ in self.stem_positions()]) / len(self.n_stems)


	#? Base pairs

	@property

	def bp_number(self):

		''' The total number of base pairs '''

		return (self.secondary.count('(') + self.secondary.count(')')) / 2.0
	

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


	def bp_freqs(self):

		''' The proportion of each type of base pairing divide by the total number of base pairs '''

		freqs = self.bp_counts()

		nfreqs = {}

		for base_pair in freqs:

			nfreqs[base_pair + '\\bp_number'] = freqs[base_pair] / self.bp_number

		return nfreqs

	
	def bp_stems(self):

		''' The proportion of each type of base pairing divide by the total number of stems '''

		freqs = self.bp_counts()

		nfreqs = {}

		for base_pair in freqs:

			nfreqs[base_pair + '\\n_stems'] = freqs[base_pair] / self.n_stems

		return nfreqs


	def avg_bp_stems(self):

		''' The average base pair proportion in the stems '''

		n = self.n_stems

		freqs = {'AU':0,'GC':0,'GU':0}

		for init1, end1, init2, end2 in self.stem_positions():

			for nt1, nt2 in zip(self.sequence[init1:end1+1], self.sequence[init2:end2+1][::-1]):

				try:
					
					freqs[nt1 + nt2]+=1

				except KeyError:

					try:

						freqs[nt2 + nt1]+=1

					except KeyError:

						raise Exception('Could not account for this base pair: {}'.format(nt1 + nt2))


		return {'avg_' + base_pair + r'\n_stems':freqs[base_pair] / n for base_pair in freqs}


	#? Bulges

	@property

	def n_bulges(self):

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

	
	def total_nt_bulges(self):

		return sum([end - init for init, end in self.bulges_pos()])


	def avg_bulge_size(self):

		return self.total_nt_bulges() / self.seqlen


	#? Loops

	@property

	def n_loops(self):

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


	def longest_loop(self):

		return max([(end1 - init1) + (end2 - init2) for init1, end1, init2, end2 in self.loops_pos()])


	def asym_loops(self):

		c = 0

		for init1, end1, init2, end2 in self.loops_pos():

			if end1 - init1 != end2 - init2:

				c+=1
		
		return c


	def sym_loops(self):

		c = 0

		for init1, end1, init2, end2 in self.loops_pos():

			if end1 - init1 == end2 - init2:

				c+=1
		
		return c


	def nt_asym_loops(self):

		sizes = []

		for init1, end1, init2, end2 in self.loops_pos():

			if end1 - init1 != end2 - init2:

				sizes.append((end1 - init1) + (end2 - init2))
		
		return sum(sizes) / len(sizes)


	def nt_sym_loops(self):

		sizes = []

		for init1, end1, init2, end2 in self.loops_pos():

			if end1 - init1 == end2 - init2:

				sizes.append((end1 - init1) + (end2 - init2))
		
		return sum(sizes) / len(sizes)


	#! miRNA related features

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



	#! Thermodinamics related features

	@property

	def dq(self):

		soma = 0

		for i in self.bpp:

			for j in self.bpp[i]:

				if i < j:

					pb = self.bpp[i][j]

					soma += pb * math.log(pb, 2)

		return -1 * soma / self.seqlen


	#! Image functions


	def __rnafold(self, folder='./'):

		''' Runs the RNAfold secondary structure prediction with the partition function and
		pairing probability matrix calculation. Also, parses the output and can direct the
		generated PS image file to a specific path. '''

		rnafold = subprocess.Popen(['RNAfold -p --noPS'],
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

		pattern = re.compile(r'[.()}{|,]+(?=\s\W[\s|-])|-?\d+\.\d+')

		matches = pattern.findall(data)

		if len(matches) != 9:

			print(f'ERROR!\n{data}')

			raise Exception('Could not parse the RNAfold output')

		secondary, mfe, ppnotation, efe, ctd, ctdenergy, ctddist, freq, div = matches


		# The predicted secondary strucuture with the lowest energy in dot-bracket notation and its MFE value

		self.secondary = secondary

		self.mfe = float(mfe)

		# Pseudo bracket notation of pair probabilities and ensemble free energy (EFE)

		self.pseudo = ppnotation

		self.efe = float(efe)

		# Centroid ensemble structure dot bracket notation, its free energy and distance from the ensemble

		self.centroid = ctd

		self.ctdenergy = float(ctdenergy)

		self.ctddist = float(ctddist)

		# The frequency of the MFE structure and the ensemble strucutral diversity (mean base pair distance)

		self.freq = float(freq)

		self.diversity = float(div)

		# Parsing the RNAfold name_dp.ps for the base pairing probabilites

		self.bpp = dfd(dict)

		with open(f'{self.name}_dp.ps') as dp:

			for line in dp:

				if line == '%start of base pair probability data\n':

					break

			for line in dp:
									
				try:
					
					i, j, pb, ubox = line[:-1].split(' ')

					if ubox == 'lbox':

						break

				except ValueError:

					print(line)

					raise Exception('Could not parse the base pair probability file from RNAfold')
				
				self.bpp[int(i)][int(j)] = float(pb) ** 2  


	@property

	def mfeden(self):

		if self.seqlen >= 40 and self.seqlen <= 600:

			return round(100 * (self.mfe - extra.refmfe[self.seqlen]) / (self.seqlen - extra.SHIFT_CONST), 2)

		else:

			return 'N.A.' 


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


	def features(self):

		features = copy(vars(self))
		
		for nonfeature in ['name', 'sequence', 'mirna1', 'mirna2', 'pos1', 'pos2', 'secondary', 'pseudo', 'centroid']:

			features.pop(nonfeature)

		features.update(self.nt_freqs())

		features.update(self.dint_freqs())

		features.update(self.triplets())

		#if self.n_stems:
		
		features.update(self.bp_freqs())

		features.update(self.bp_stems())

		features.update(self.avg_bp_stems())


		return features

	
	#! Image functions


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



