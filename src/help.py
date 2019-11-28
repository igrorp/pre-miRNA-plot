desctxt = '''
 python3+ premiRNA-plot.py [INPUTS] ... [OPTIONS] ...

-----------------------------------------------------------------------------------------

Pre-miRNA-plot is a program for generating multiple custom images of miRNA precursors based
on RNAfold and RNAplot. It allows you to highlight the miRNA location within the precursor 
and obtain general and practical information about your data, so you can filter it or use 
it in publications.

-----------------------------------------------------------------------------------------

The program accepts tab-separated text files containing possibly 4 columns:

	1) Sequence ID (optional): Some sort of annotation or ID
		information about the sequence (e.g. 'ath-miR-171' or 'seq1');
	2) Precursor sequence: The pre-miRNA sequence.
	3) miRNA1 sequence: One of the miRNAs sequences.
	4) miRNA2 sequence (optional): The other miRNA sequence.
	
	File format examples:
	
	| pre-miRNA |   miRNA1   |  miRNA2  |
	
	>>>>>>>>>>>>>>>>>>>>>> OR <<<<<<<<<<<<<<<<<<<<<<<
	
	|     ID    |  pre-miRNA |  miRNA1  |   miRNA2  |

There is no problem if you don't have both miRNAs sequences; 
you can inform just one and the program will work just fine. 

Checkout our github repository for more details on how to use
the program: https://github.com/igrorp/pre-miRNA-plot

Parameters:

-i, --input       (str ...)
	 
	Inform the names of/paths to the files that you want to use.

-a, --annot       (str) --> (T or F)
	
	Informs if you have some sort of sequence ID, such as a miRNA
	family annotation (e.g.'ath-miRNA-171', 'seq1'), necessarily 
	on the first column, so that the generated image files can be 
	named according to that ID.

	Default = F (False)

-s, --style       (int) --> (between 1 and 5)

	The style of created images. Check the repository to see how they
	look like.

	Default = 3

-c, --color       (str str) OR (int int int int int int)
	
	You can choose which colors to paint the miRNA sequence within
	the precursor. Always provide the 5p and 3p colors, respectively.
	You can choose the predefined color names blue, red, green, purple,
	pink, yellow, cyan, white, black and orange; or you can inform the 
	RGB codes of the colors you want.

	Ex: '-c blue green' for blue 5p and green 3p
	Ex: '-c 255 255 0 153 0 204' for yellow 5p and purple 3p
	Default = green and red

-t, --threads     (int)
	
	Choose the number of allowed processors (CPUs) to be used.
	
	Default = 1

-f, --formats     (str) --> (svg or pdf)

	Choose the output format of the images. Choose between PDF and SVG.

	Default = svg

-o, --outdir      (str)

	The output name of directory created containing all the generated data.

	Default = premirnaplot



---------------------------------------------------------------

If you have any comments, complaints, doubts or suggestions, please contact
our main responsible for this project at igorpaim8@gmail.com or create an issue in
out github repository https://github.com/igrorp/pre-miRNA-plot/issues.

Have a nice work and let's keep making science evolve!

'''


defcolors = {
	'blue':'#0000cc',
	'red':'#ff0000',
	'green':'#33cc33',
	'purple':'#9933ff',
	'pink':'#ff66cc',
	'yellow':'#ffff00',
	'cyan':'#00ffff',
	'white':'#ffffff',
	'black':'#000000',
	'orange':'#ff9933'
}

class Precursor():

    def __init__(self, name, mirna1, mirna2, precursor):

        # The prefix of the filename for the pre-miRNA
        self.name = name

        # The precursor sequence
        self.premirna = precursor
        
        # The tuple containing the positions of the miRNA within the pre-miRNA
        self.pos1 = self.pos(precursor, mirna1)

        # The tuple containing the positions of the other miRNA within the pre-miRNA
        self.pos2 = self.pos(precursor, mirna2)

        # The length of the pre-miRNA
        self.prelen = len(precursor)

        # The Minimum Free Energy for the precursor as predicted by RNAfold
        self.premfe = 0


    def pos(self, premirna, mirna):

        if mirna:

            if mirna in premirna:

                if premirna.count(mirna) > 1:
                    print('WARNING! miRNA {} was found more than once in the precursor sequence {}..., but its last occurrence will be used!'.format(mirna, premirna[25:]))

                return (premirna.find(mirna), premirna.find(mirna) + len(mirna))
            
            else:
                
                raise Exception('ERROR! Could not find sequence {} inside {}, please correct this'.format(mirna, premirna))

        else:

            return None

