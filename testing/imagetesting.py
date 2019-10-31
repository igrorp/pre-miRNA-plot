
import xml.etree.ElementTree as et

import statistics as st

import math

import svgwrite as sw

class NoLocationAnnotations(Exception):
	pass

class NoPairsAnnotations(Exception):
	pass

class NoSequenceAnnotation(Exception):
	pass

class UnimplementedParser(Exception):
	pass

class UnparsableSVG(Exception):
	pass


class SVGParser():

	"""This is a class to parse the SVG files produced by RNAplot. """

	def __init__(self, file):

		sequence = None
		locations = []
		box = ()
		pairs = []
		transform = None
		seqtransform = None
		radius = 0

		self.sequence, self.pairs = self.load_data(file)

		self.radius = self.__radius__()

		if not self.pairs:
			raise NoPairsAnnotations("Did not find any pairs")

		if self.sequence is None:
			raise NoSequenceAnnotation("Did not find the sequence")

		if not self.locations:
			raise NoLocationAnnotations("Did not find drawing coordinates")

	
	def load_data(self, file):
		tree = et.parse(file)
		root = tree.getroot()

		for child in root:
			self.transform = child.attrib.get('transform', None)

		container = None
		for child in root:
			if child.tag == '{http://www.w3.org/2000/svg}g':
				container = child
				break

		if container is None:
			raise UnparsableSVG("Cannot find container")

		self.box = self.__box__(root)
		sequence, self.locations = self.__locations__(container)
		
		return sequence, self.__pairs__(container)


	def __pairs__(self, root):
		
		pair_node = None
		size = None

		for child in root:
			if child.attrib.get('id', None) == 'pairs':
				pair_node = child
			if child.attrib.get('id', None) == 'seq':
				size = len(child)

		if pair_node is None:
			raise NoPairsAnnotations("Couldn't find the pairs")

		if size is None:
			raise NoSequenceAnnotation("Couldn't find the sequence")

		pairs = [None] * size
		
		for pair in pair_node:
			pair  = list(map(lambda i: int(i) - 1, pair.attrib['id'].split(',')))
			pairs[pair[0]] = pair[1]
			pairs[pair[1]] = pair[0]

		return pairs


	def __radius__(self):

		x = 0
		a = []

		while x < len(self.locations) - 1:
			x1, y1 = self.locations[x]
			x2, y2 = self.locations[x+1]
			x+=1
			a.append(math.sqrt(((y2 - y1) ** 2) + ((x2 - x1) ** 2)))

		n = 0

		for pair in self.pairs:
			if pair: n+=1

		n = n / len(self.pairs)

		return (st.mean(a) + (st.stdev(a) * n)) / 2


	def __box__(self, root):

		''' Returns the width and heigth of the SVG box'''

		return (int(root.attrib['width']), int(root.attrib['height']))

	
	def __locations__(self, root):
		
		''' Returns the coordinates of the nucleotides text from the SVG
		and the group transform value '''

		seq_node = None
		width, height = self.box
		transform = self.transform
		
		for child in root:
			if child.attrib.get('id', None) == 'seq':
				seq_node = child
				break

		if seq_node is None:
			raise NoSequenceAnnotation("Could not find sequence")

		xyt = transform[transform.find('translate') + len('translate'):]
		xt = round(float(xyt.split(',')[0][1:]), 5)
		yt = round(float(xyt.split(',')[1][:-1]), 5)

		scale = transform[transform.find('scale') + len('scale'):transform.find(' translate')]
		scalex = round(float(scale.split(',')[0][1:]), 5)
		scaley = round(float(scale.split(',')[1][:-1]), 5)

		sequence = ''
		locations = []
		
		for node in seq_node:
			sequence+=node.text
			locations.append((round((float(node.attrib['x']) + xt) * scalex, 5),
				              round((float(node.attrib['y']) + yt) * scaley, 5)))
		
		return sequence, locations


some = SVGParser(open('rna.svg'))

print(some.box)
print(some.sequence)
print(some.pairs)
print(some.locations)
print(some.transform)

radius = some.radius

print(f'this is the {radius}')

width, height = some.box

method = 'lines'

fs = radius * 2 * 0.725

fdy = radius / 2

fdx = radius / 2 * -1 * 1.13

if method == 'spheres':
	strokew = radius / 3
elif method == 'lines':
	strokew = radius * 2 * 0.725

mirpos = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
mir2pos = [69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91]


dwg = sw.Drawing('new.svg', viewBox=f"0, 0, {width}, {height}", preserveAspectRatio="xMidYMid meet")

# Criando o grupo que vai conter todas as estruturas juntas

precursor = dwg.add(dwg.g(id="precursor", transform='translate(0, 10) scale(0.95, 0.95)'))

# Criando todos os grupos

pairsgroup = precursor.add(dwg.g(id='pairs', stroke='black', stroke_width=radius / 6))

circgroup = precursor.add(dwg.g(id='circles', display='none', fill="#9494b8"))

# Adicionando a linha que conecta todos os nucleotideos

precursor.add(dwg.polyline(points=some.locations, stroke='#9494b8', fill="none", stroke_width=strokew))

# Adicionando as linhas que conectam os nucleotídeos pareados

for i in range(len(some.pairs)):
	if some.pairs[i]:
		x1, y1 = some.locations[i]
		x2, y2 = some.locations[int(some.pairs[i])]
		pairsgroup.add(dwg.line(start=some.locations[i], end=some.locations[int(some.pairs[i])])) 

# Adicionando a linha do estilo de traço que destaca a linha da posição dos miRNAs

linha1 = []
linha2 = []

for index, xytuple in enumerate(some.locations):

	if index in mirpos:
		x, y = some.locations[index]
		if index == mirpos[0]:
			y-=radius
		elif index == mirpos[-1]:
			y+=radius
		linha1.append((x,y))
	elif index in mir2pos:
		x, y = some.locations[index]
		if index == mir2pos[0]:
			y+=radius
		elif index == mir2pos[-1]:
			y-=radius
		linha2.append((x,y))

precursor.add(dwg.polyline(points=linha1, fill='none', stroke='red', stroke_width=strokew))
precursor.add(dwg.polyline(points=linha2, fill='none', stroke='orange', stroke_width=strokew))

textgroup = precursor.add(dwg.g(id='nucleotides', transform=f'translate({fdx},{fdy})',font_size=fs,
								fill='black', font_family='Helvetica', font_weight='bold'))

# Adicionando cada um dos círculos do estilo de esferas

for index, xytuple in enumerate(some.locations):
	
	circ = circgroup.add(dwg.circle(center=xytuple, r=radius))

	textgroup.add(dwg.text(some.sequence[index], insert=some.locations[index]))
	
	if index in mirpos:
		circ.fill('red')
	elif index in mir2pos:
		circ.fill('orange')

# Salvando o arquivo com identação

dwg.save(pretty=True)

