
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



class SVGconstructor(SVGParser):


	def __init__(self, file, style):
		
		''' Precisa criar o elemento dwg e todos os grupos'''

		dwg = None

		precursor = None

		super().__init__(file)

		width, height = self.box

		self.dwg = sw.Drawing('new.svg', viewBox=f"0, 0, {width}, {height}", preserveAspectRatio="xMidYMid meet")

		self.precursor = self.dwg.add(self.dwg.g(id="precursor", transform='translate(0, 10) scale(0.95, 0.95)'))

		self.drawpairs()

		mirpos = [5, 25]
		mir2pos = [69, 91]

		if style == '1':
			# Desenhar pares
			# Desenhar única polyline fina e preta
			self.polyline([0, len(self.locations) - 1], 'black', self.radius / 4)
			# Desenhar circulos bg e colored
			# Desenhar texto em tudo
			self.dwg.save(pretty=True)
			pass
		elif style == '2':
			# Desenhar pares
			# Desenhar única polyline média e bg
			self.polyline([0, len(self.locations) - 1], 'grey', self.radius / 2)
			# Desenhar circulos colored

			# Desenhar apenas texto dos mirnas
			pass
		elif style == '3':
			# Desenhar pares
			# Desenhar unica polyline fina e preta
			self.polyline([len(self.locations) - 1, 0], 'black', self.radius / 4)
			# Desenhar circulos bg e colored com stroke
			self.circgroup = self.precursor.add(self.dwg.g(id='circles', fill='white', stroke='black'))
			self.textgroup = self.precursor.add(self.dwg.g(id='nucleotides', transform=f'translate({fdx},{fdy})', font_size=fs, fill='black', font_family='Helvetica', font_weight='bold'))
			self.drawcircles([0, int(len(self.locations) / 2)], mirpos, 'red')
			self.drawcircles([int(len(self.locations) / 2), len(self.locations)], mir2pos, 'orange')
			self.dwg.save(pretty=True)
			# Desenhar texto em tudo
			pass
		elif style == '4':
			# Desenhar pares
			# Desenhar polyline unica e media preta
			self.polyline([0,len(self.locations) - 1], 'black', self.radius / 2)
			# Desenhar polyline duplas grossas
			self.polyline(mirpos, 'red', 9)
			self.polyline(mir2pos[::-1], 'orange', 9)
			# Desenhar texto com fonte menorzinha
			self.dwg.save(pretty=True)
			pass
		elif style == '5':
			# Desenhar pares
			# Desenhar polyline unica grossa com bg
			self.polyline([0, len(self.locations) - 1], 'black', 9)
			# Desenhar polyline duplas grossa
			self.polyline(mirpos, 'red', 9)
			self.polyline(mirpos, 'orange', 9)
			# Desenhar texto em tudo com fonte menorzinha
			pass
		else:
			raise Exception("Could not identify the style of the image, choose between 1-5")

	
	def __properties__(self):

		''' Setting up the properties of the image, such as width, height,
		font size, stroke width, etc '''

		pass


	def drawcircles(self, poslst, spclst, color):

		''' Creates circles for the given list of positions '''

		for index in range(poslst[0], poslst[1]):
	
			circ = self.circgroup.add(self.dwg.circle(center=self.locations[index], r=self.radius))

			self.textgroup.add(dwg.text(some.sequence[index], insert=some.locations[index]))
			
			if index >= spclst[0] and index <= spclst[1]:
				circ.fill(color)


	def drawpairs(self):

		pairsgroup = self.precursor.add(self.dwg.g(id='pairs', stroke='black', stroke_width=self.radius / 6))

		for i in range(len(self.pairs)):
			if self.pairs[i]:
				x1, y1 = self.locations[i]
				x2, y2 = self.locations[int(self.pairs[i])]
				pairsgroup.add(self.dwg.line(start=self.locations[i], end=self.locations[int(self.pairs[i])])) 


	def text():
		pass
		# textgroup = self.precursor.add(self.dwg.g(id='nucleotides', transform=f'translate({fdx},{fdy})',font_size=fs, fill='black', font_family='Helvetica', font_weight='bold'))
		

	def polyline(self, poslst, color, strokew):
		
		newpos = []

		for index, xytuple in enumerate(self.locations):

			if index >= min(poslst) and index <= max(poslst):
				x, y = self.locations[index]
				if index == poslst[0]:
					y-=self.radius
				elif index == poslst[1]:
					y+=self.radius
				
				newpos.append((x,y))

		self.precursor.add(self.dwg.polyline(points=newpos, fill='none', stroke=color, stroke_width=strokew))


 


some = SVGconstructor(open('rna2.svg'), '3')

print(some.box)
print(some.sequence.index('AGCAGUUU') - len('AGCAGUUU'))
print(some.pairs)
print(some.locations)
print(some.transform)

# radius = some.radius

# print(f'this is the {radius}')

# width, height = some.box

# method = 'spheres'

# fs = radius * 2 * 0.725

# fdy = radius / 2

# fdx = radius / 2 * -1 * 1.13

# if method == 'spheres':
# 	strokew = radius / 3
# elif method == 'lines':
# 	strokew = radius * 2 * 0.725





# # Criando o grupo que vai conter todas as estruturas juntas


# # Criando todos os grupos

# pairsgroup = precursor.add(dwg.g(id='pairs', stroke='black', stroke_width=radius / 6))



# # Adicionando a linha que conecta todos os nucleotideos

# precursor.add(dwg.polyline(points=some.locations, stroke='#9494b8', fill="none", stroke_width=strokew))

# # Adicionando as linhas que conectam os nucleotídeos pareados


# # Adicionando a linha do estilo de traço que destaca a linha da posição dos miRNAs



# # Adicionando cada um dos círculos do estilo de esferas


		
# # Salvando o arquivo com identação

# dwg.save(pretty=True)

