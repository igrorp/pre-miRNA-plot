
import xml.etree.ElementTree as et

import svgwrite as sw

class NoLocationAnnotations(Exception):
    """This error is raise when we cannot find the part of the postscript which
    specifies where to put the nucleotides.
    """
    pass

class NoPairsAnnotations(Exception):
    """This error is raised when we cannot find the pairing information in the
    psotscript.
    """
    pass

class NoSequenceAnnotation(Exception):
    """This exception is raised if we cannot find the part of the postscript
    that specifies the sequence.
    """
    pass

class UnimplementedParser(Exception):
    """This is raised when we need a parser that we have not yet implemented to
    parse the results of RNAplot.
    """
    pass

class UnparsableSVG(Exception):
    """This is raised if we can't parse the SVG file for any reason.
    """
    pass

class Parser():
    """This is the basic parser for all parsers produced by RNAplot.
    """
    #__metaclass__ = abc.ABCMeta

    def __init__(self, stream):
        sequence = None
        self.locations = []
        """The locations of coordinates to draw."""
        self.box = ()
        """The bounding box of the drawing."""

        sequence, pairs = self.load_data(stream)

        if not pairs:
            raise NoPairsAnnotations("Did not find any pairs")

        if sequence is None:
            raise NoSequenceAnnotation("Did not find the sequence")

        if not self.locations:
            raise NoLocationAnnotations("Did not find drawing coordinates")

        #super(Parser, self).__init__(pairs, sequence=sequence)

    #@abc.abstractmethod
    def load_data(self, stream):
        """This method should load all data from the stream. It should return
        the pairs and set self.sequence, self.locations, and self.box for this
        object.

        :stream: The data stream to read.
        :returns: The pairs.
        """
        pass

class SVGParser():
    """This is a class to parse the svg files produced by RNAplot.
    """


    def __init__(self, stream):

        sequence = None
        locations = []
        """The locations of coordinates to draw."""
        box = ()
        pairs = []
        """The bounding box of the drawing."""

        self.sequence, self.pairs = self.load_data(stream)

        if not self.pairs:
            raise NoPairsAnnotations("Did not find any pairs")

        if self.sequence is None:
            raise NoSequenceAnnotation("Did not find the sequence")

        if not self.locations:
            raise NoLocationAnnotations("Did not find drawing coordinates")

    def load_data(self, stream):
        tree = et.parse(stream)
        root = tree.getroot()
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
        return (int(root.attrib['width']), int(root.attrib['height']))

    def __locations__(self, root):
        # TODO: Deal with the required transforms
        seq_node = None
        for child in root:
            if child.attrib.get('id', None) == 'seq':
                seq_node = child
                break

        if seq_node is None:
            raise NoSequenceAnnotation("Could not find sequence")

        sequence = []
        locations = []
        for node in seq_node:
            sequence.append(node.text)
            locations.append((float(node.attrib['x']),
                              float(node.attrib['y'])))
        return ''.join(sequence), locations

some = SVGParser(open('rna2.svg'))

print(some.box)
print(some.sequence)
print(some.pairs)
print(some.locations)

width, height = some.box

print(width, height)

dwg = sw.Drawing('new.svg', size=some.box, viewBox= f'0 0 {width} {height}')

# Criando o grupo que vai conter todas as estruturas juntas

precursor = dwg.add(dwg.g(id="precursor", transform="scale(0.7, 0.7) translate(206.680191,99.463547)"))

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

textgroup = precursor.add(dwg.g(id='nucleotides', transform='translate(-3.7,3.7)',font_size=11, fill='black', font_family='Helvetica', font_weight='bold'))

linha1 = []
linha2 = []

for index, xytuple in enumerate(some.locations):
	
	circ = circgroup.add(dwg.circle(center=xytuple, r=8, stroke='none', stroke_width='1.5'))

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

