desctxt = '''
 python3+ premiRNA-plotting.py [OPTIONS] ... [INPUTS] ...

-----------------------------------------------------------------------------------------


                        Hello, welcome to pre-mirRNA-plotting!

   We've developed this Python program to facilitate the high-throughput generation of
secondary structure prediction images of miRNA precursors using RNAfold and RNAplot.
One of our main intentions as well is to highlight the position of the miRNAs sequences
within the precursor, as it is an important criteria for miRNA selection and filtering
after the prediction was done.

   The program will generate a folder called 'premirnaplot', where your data will be.
To any valid input file given, a folder will be created with it's name and inside each
one there are going to be another two folders: 'colored_structures' and 'raw_folding'.
One has the secondary structures images with the miRNAs sequences highlighted and the 
other has the TXT files containing the prediction by RNAfold and the uncolored image of
the pre-miRNA structure. The miRNAs are highlighted in green (5p) and red (3p).

   Please read the help info below to have more details about the parameters that best fit
your data and intentions.

-----------------------------------------------------------------------------------------

This program accepts tab-separated text files containing
possibly 4 columns:

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

Check out the parameters descriptions to have more information
about the arguments that best fit your data and intentions.

Also, you can inform more than one file. We suggest that the
filenames should be somehow informative because they will
be used for naming the generated folder and written with the
file format in the end, such as 'homo_sapiens.txt' or 'hsa.txt'.


Parameters:


 -i, --input

	 Inform the names of the files that you wish to use.


 -a, --annot    <T or F>

	 Informs if you have some sort of sequence ID, such as a miRNA
	 family annotation (e.g.'ath-miRNA-171', 'seq1'), necessarily 
	 on the first column, so that the generated image files can be 
	 named according to that ID.

	 Default = False


 -e, --extra    <T or F>

	 Provides additional files with basic data about the informed 
	 pre-miRNAs: one boxplot with their minimum free energy values,
	 as calculated by RNAfold, and another boxplot of their length.

	 Default = True


 -c, --color

	 COLOR COLOR or RGBCODE RGBCODE

	 You can choose which colors to paint the miRNA sequence within
	 the precursor. Always provide the 5p and 3p colors, respectively.
	 You can choose the predefined color names blue, red, green, purple,
	 pink, yellow, cyan, white, black and orange; or you can inform the 
	 RGB codes of the colors you want.

	 Ex: '-c blue green' for blue 5p and green 3p
	 Ex: '-c 255 255 0 153 0 204' for yellow 5p and purple 3p

	 Default = green and red


 -q, --quality

	 Choose the quality (in dpi) of the generated images.
	
	 Default = 200 dpi


-t, --threads

	 Choose the number of allowed threads (CPUs) to run.

	 Default = 1


---------------------------------------------------------------


Thanks for using our program! 

If you have any comments, complaints, doubts or suggestions. Please don't
hesitate to contact our main responsible for this project at igorpaim8@gmail.com.

Have a nice work and let's keep making science evolve!


'''


defcolors = {
	'blue':['0', '0.4', '1.'],
	'red':['1.', '0', '0'],
	'green':['0', '0.6', '0'],
	'purple':['0.8', '0.2', '1.'],
	'pink':['1.', '0', '0.4'],
	'yellow':['1.', '0.8', '0'],
	'cyan':['0.2', '0.8', '1.'],
	'white':['1.', '1.', '1.'],
	'black':['0', '0', '0'],
	'orange':['1.', '0.458', '0.102']
}