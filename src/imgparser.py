
from svglib.svglib import svg2rlg

from reportlab.graphics import renderPDF

import xml.etree.ElementTree as et

import statistics as st

import math

import svgwrite as sw

import subprocess


class SVGParser():

	""" Parses the SVG files produced by RNAplot """

	def __init__(self, file):

		sequence = None
		locations = []
		box = ()
		pairs = []
		transform = None
		seqtransform = None
		radius = 0

		self.sequence, self.pairs = self.load_data(file)

		if not self.pairs:
			raise Exception("Did not find any pairs")

		if not self.sequence:
			raise Exception("Did not find the sequence")

		if not self.locations:
			raise Exception("Did not find drawing coordinates (locations)")

	
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
			
			raise Exception("Could not parse SVG: Cannot find container")

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
			raise Exception("Couldn't find the pairs")

		if size is None:
			raise Exception("Couldn't find the sequence")

		pairs = [None] * size
		
		for pair in pair_node:
			pair  = list(map(lambda i: int(i) - 1, pair.attrib['id'].split(',')))
			pairs[pair[0]] = pair[1]
			pairs[pair[1]] = pair[0]

		return pairs


	def __box__(self, root):

		''' Returns the width and heigth of the SVG box'''

		return (float(root.attrib['width']), float(root.attrib['height']))

	
	def __locations__(self, root):
		
		''' Returns the coordinates of the nucleotides text from the SVG
		and the group transform value '''

		seq_node = None
		transform = self.transform
		
		for child in root:
			if child.attrib.get('id', None) == 'seq':
				seq_node = child
				break

		if seq_node is None:
			raise Exception("Could not find sequence")

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


	def __init__(self, filepath, style, mirpos, mir2pos, color1, color2, pdf=False):

		dwg = None

		precursor = None

		if mirpos is None: mirpos = (0,0)
		if mir2pos is None: mir2pos = (0,0)

		super().__init__(open(filepath))

		self.__properties__(style)

		width, height = self.box

		self.dwg = sw.Drawing(filepath, viewBox='0, 0, {}, {}'.format(width, height), preserveAspectRatio='xMidYMid meet', debug=False)

		self.precursor = self.dwg.add(self.dwg.g(id='precursor', transform='translate(0, 10) scale(0.95, 0.95)'))

		self.drawpairs()


		if style == 1:

			self.polyline([0, len(self.locations) - 1], 'black', self.strokew)

			self.circgroup = self.precursor.add(self.dwg.g(id='circles', fill='#9494b8'))
			
			self.textgroup = self.precursor.add(self.dwg.g(id='nucleotides', transform='translate(0,{})'.format(self.fdy),
												font_size=self.fs, fill='black', font_family='Helvetica', font_weight='bold'))

			self.drawcircles((0, len(self.locations)), '#9494b8')
			self.drawcircles(mirpos, color1)
			self.drawcircles(mir2pos, color2)

			self.drawtext([0, len(self.locations)])

		elif style == 2:

			self.polyline([0, len(self.locations) - 1], 'black', self.strokew)

			self.circgroup = self.precursor.add(self.dwg.g(id='circles', fill='#9494b8'))
			
			self.textgroup = self.precursor.add(self.dwg.g(id='nucleotides', transform='translate(0,{})'.format(self.fdy),
												font_size=self.fs, fill='black', font_family='Helvetica', font_weight='bold'))

			self.drawcircles(mirpos, color1)
			self.drawcircles(mir2pos, color2)

			self.drawtext([mirpos[0], mirpos[1]])
			self.drawtext([mir2pos[0], mir2pos[1]])

		elif style == 3:

			self.polyline([len(self.locations) - 1, 0], 'black', self.strokew)
			
			self.circgroup = self.precursor.add(self.dwg.g(id='circles', fill='white', stroke='black'))
			
			self.textgroup = self.precursor.add(self.dwg.g(id='nucleotides', transform='translate(0, {})'.format(self.fdy),
												font_size=str(self.fs), fill='black', font_family='Helvetica', font_weight='bold'))
			
			self.drawcircles((0, len(self.locations)), 'white')
			self.drawcircles(mirpos, color1)
			self.drawcircles(mir2pos, color2)

			self.drawtext([0, len(self.locations)])
			
		elif style == 4:

			self.polyline([0,len(self.locations) - 1], 'black', self.strokew)

			self.polyline(mirpos, color1, self.fs * 1.15)
			self.polyline(mir2pos[::-1], color2, self.fs * 1.15)

			self.textgroup = self.precursor.add(self.dwg.g(id='nucleotides', transform='translate(0, {})'.format(self.fdy),
												font_size=self.fs, fill='black', font_family='Helvetica', font_weight='bold'))

			self.drawtext([mirpos[0], mirpos[1]])
			self.drawtext([mir2pos[0], mir2pos[1]])

		elif style == 5:

			self.polyline([0, len(self.locations)], '#9494b8', self.strokew)

			self.polyline(mirpos, color1, self.strokew)
			self.polyline(mir2pos[::-1], color2, self.strokew)

			self.textgroup = self.precursor.add(self.dwg.g(id='nucleotides', transform='translate(0, {})'.format(self.fdy),
									font_size=self.fs, fill='black', font_family='Helvetica', font_weight='bold'))

			self.drawtext([0, len(self.locations)])

		else:

			raise Exception("Could not identify the style of the image, choose between 1-5")
		
		
		self.dwg.saveas(filepath, pretty=True)

		if pdf:
			drawing = svg2rlg(filepath)
			renderPDF.drawToFile(drawing, filepath[:-3] + 'pdf')
			subprocess.run('rm {}'.format(filepath), shell=True)


	def __properties__(self, style):

		''' Setting up the properties of the image, such as width, height,
		font size, stroke width, etc '''

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

		self.radius = (st.mean(a) + (st.stdev(a) * n)) / 2

		if style == 3:
			self.circsw = self.radius * 0.1613
			self.radius-=self.circsw

		self.fs = self.radius * 2 * 0.725

		if style == 4 or style == 5:
			self.fs = self.radius * 2 * 0.625

		self.fdy = self.fs * 0.345

		if style == 1 or style == 3:
			self.strokew = (self.radius * 2) / 6
		elif style == 2 or style == 4:
			self.strokew = self.fs / 2
		elif style == 5:
			self.strokew = self.fs * 1.15


	def drawcircles(self, postup, color):

		''' Creates circles for the given list of positions '''

		init, fin = postup

		for index in range(init, fin):
		
			self.circgroup.add(self.dwg.circle(center=self.locations[index], r=self.radius, fill=color))


	def drawpairs(self):

		pairsgroup = self.precursor.add(self.dwg.g(id='pairs', stroke='black', stroke_width=self.radius / 6))

		for i in range(len(self.pairs)):
			
			if self.pairs[i]:
				
				pairsgroup.add(self.dwg.line(start=self.locations[i], end=self.locations[int(self.pairs[i])])) 


	def drawtext(self, postup):

		''' Creates the nucletotides letters for the given tuple of positions '''

		init, fin = postup

		for index in range(init, fin):

			self.textgroup.add(self.dwg.text(self.sequence[index], insert=self.locations[index], text_anchor='middle'))


	def as_string(self):

		self.dwg.tostring()
					

	def polyline(self, postup, color, strokew):
		
		newpos = []

		init, fin = postup

		for index, _ in enumerate(self.locations):

			if index >= min(postup) and index <= max(postup):
				x, y = self.locations[index]
				if index == init:
					y-=self.radius
				elif index == fin:
					y+=self.radius
				
				newpos.append((x,y))

		self.precursor.add(self.dwg.polyline(points=newpos, fill='none', stroke=color, stroke_width=strokew))


## Example on how you would parse and create a new SVG image using the constructor

# SVGconstructor('rna2.svg', '3', (5, 25), (69, 91), '#cc33ff', '#ffff00', 'grey')

## Example on how you would parse and create a new SVG image and convert it to PDF

# SVGconstructor('rna2.svg', '3', (5, 25), (69, 91), '#cc33ff', '#ffff00', 'grey', pdf=True)
