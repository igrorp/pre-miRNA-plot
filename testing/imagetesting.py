
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

		self.sequence, self.pairs = self.load_data(file)

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

		print(transform)

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


some = SVGParser(open('rna2.svg'))

print(some.box)
print(some.sequence)
print(some.pairs)
print(some.locations)
print(some.transform)

width, height = some.box


dist = 0
x = 0
a = []

while x < len(some.locations) - 1:
	x1, y1 = some.locations[x]
	x2, y2 = some.locations[x+1]
	x+=1
	a.append(math.sqrt(((y2 - y1) ** 2) + ((x2 - x1) ** 2)))

n = 0

for pair in some.pairs:
	if pair: n+=1

n = n / len(some.pairs)

radius = (st.mean(a) + (st.stdev(a) * n)) / 2

fs = radius * 0.725

dwg = sw.Drawing('new.svg', viewBox=f"0, 0, {width}, {height}", preserveAspectRatio="xMidYMid meet")

# Criando o grupo que vai conter todas as estruturas juntas

precursor = dwg.add(dwg.g(id="precursor", transform='translate(0, 10) scale(0.95, 0.95)'))

# Criando o grupo de pares entre nucleotideos conectados

pairsgroup = precursor.add(dwg.g(id='pairs', stroke='black', stroke_width='2'))

for i in range(len(some.pairs)):
	if some.pairs[i]:
		x1, y1 = some.locations[i]
		x2, y2 = some.locations[int(some.pairs[i])]
		pairsgroup.add(dwg.line(start=some.locations[i], end=some.locations[int(some.pairs[i])])) 


# Adicionando a linha que conecta todos os nucleotideos do precursor

precursor.add(dwg.polyline(points=some.locations, stroke='black', fill="none", stroke_width="2"))

# Criando o grupo que vai conter todos os circulos

circgroup = precursor.add(dwg.g(id='circles', fill="#9494b8"))

mirpos = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
mir2pos = [69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91]

textgroup = precursor.add(dwg.g(id='nucleotides', transform=f'translate({-1 * radius / 2 * 1.3},{radius / 2})',font_size=fs, fill='black', font_family='Helvetica', font_weight='bold'))

linha1 = []
linha2 = []

for index, xytuple in enumerate(some.locations):
	
	circ = circgroup.add(dwg.circle(center=xytuple, r=radius, stroke='none', stroke_width='1.5'))

	textgroup.add(dwg.text(some.sequence[index], insert=some.locations[index]))
	
	if index in mirpos:
		circ.fill('red')
		linha1.append(some.locations[index])
	if index in mir2pos:
		circ.fill('orange')
		linha2.append(some.locations[index])
	

#precursor.add(dwg.polyline(points=linha1, fill='red', stroke='red', stroke_width="10"))
#precursor.add(dwg.polyline(points=linha2, fill='red', stroke='orange', stroke_width="10"))


dwg.save(pretty=True)

