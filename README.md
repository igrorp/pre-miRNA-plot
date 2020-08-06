---

<h1 align="center"> Welcome to premiRNAplot!</h1>



<p>PremiRNAplot is a Python 3.5+ library with a command line interface for extracting features and generating multiple custom images of miRNA precursors. It provides over a 100 different features that can be used to train and test machine learning algorithms for the correct prediction and classification of pre-miRNAs.</p>
<p>It's based on RNAfold and RNAplot from Vienna RNA and simplifies feature obtention and calculation. It also allows you to produce beautiful, costumizable and publication-ready images of the secondary structure of these precursors, highlighting the position of the miRNAs sequences within them.</p>
<ol>
	<li> <a href="#installation">Installation</a><br>
		1.1 <a href="#using-conda">Using Conda</a><br>
		1.2 <a href="#using-docker">Using Docker</a><br>
		1.3 <a href="#using-pip">Using pip</a><br>
		1.4 <a href="#using-this-repository">Using this repository</a><br>
	</li>
	<li>  <a href="#command-line-interface-usage">Command line interface</a><br>
		2.1 <a href="#input-files">Input files</a><br>
		2.2 <a href="#image-styles">Using Conda</a><br>
		2.3 <a href="#colors">Colors</a><br>
		2.4 <a href="#image-formats">Image formats</a><br>
		2.5 <a href="#multiprocessing">Multiprocessing</a><br>
	</li>
	<li> <a href="#3">Library usage</a><br>
		3.1 <a href="#31">Basic properties and modules</a><br>
		3.2 <a href="#32">Features</a><br>
	</li>
</ol>
<h1 id="installation">1. Installation</h1>
<p>There are multiple ways of downloading premiRNAplot. We highly suggest that you use the Conda package, as it would remove the need to manually install packages (ViennaRNA in specific) or having conflicts. Make sure that you have your platform updated and that you are using Python 3.5 or higher.</p>
<h2 id="using-conda">Using Conda</h2>
<p>It’s always a good pratice to create virtual environments to isolate your programs and avoid version conflicts, for example. To create a Conda virtual environment for premiRNAplot, run:</p>
(this step is optional)
<pre><code>conda create -n premirnaplot python=3.7 -y
conda activate premirnaplot
</code></pre>
<p>To install, simply run:</p>
<pre><code>conda install -c igror premirnaplot
</code></pre>
<p>The program will be available through the Bioconda channel soon.</p>
<h2 id="using-docker">Using Docker</h2>
<p>Future instructions on how to use the Docker container.</p>
<h2 id="using-pip">Using pip</h2>
<p>You can install the program from pip as well, and it will automatically install all the Python dependecies, but you need to have the Vienna RNA package manually installed to run, since it’s not a Python package. You can compile Vienna RNA from source or install it using Conda.</p>
<p>To install premirnaplot from pip (you can create a virtual environment using venv):</p>
<pre><code>pip install premirnaplot
</code></pre>
<p>To install the Vienna RNA package from Conda:</p>
<pre><code>conda install -c bioconda viennarna
</code></pre>
<p>To compile the Vienna RNA pacakge from source:</p>
<pre><code>wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz
tar xvzf ViennaRNA-2.4.13.tar.gz
cd ViennaRNA-2.4.13/
./configure
make
sudo make install
</code></pre>
<h2 id="using-this-repository">Using this repository</h2>
<p>To use this repository for the installation, clone it using:</p>
<pre><code>git clone https://github.com/igrorp/pre-miRNA-plot.git
</code></pre>
<p>Now you need to enter the cloned repository folder and run the <strong><a href="http://setup.py">setup.py</a></strong> configuration file, which will download all the necessary Python packages used in premiRNAplot.</p>
<pre><code>cd pre-miRNA-plot/
python setup.py install
</code></pre>
<p>There is going to be a lot of text being displayed informing the packages downloaded and the operation status. After that, you can run our test data to check if all the requirements are satisfied or if there are any problems during the execution.</p>
<h1 id="command-line-interface-usage">3. Command line interface usage</h1>
<p>PremiRNAplot provides a very easy and straight-forward command line interface (CLI). You can call it using <code>premirnaplot -help</code> from the terminal. Using the CLI, you'll automatically generate a folder containing the images of the secondary structure of the precursors and a file called 'precursor_data.txt' with all the features and their values.
</p><h3 id="input-files">Input files</h3>
<p>

The input files are text files with the sequences separated by <strong>tabs</strong> ‘\t' (TSVs), containing:</p>
<p><strong>1)</strong> A label to the pre-miRNA (<strong>optional</strong>, if not provided they will be labeled “precursor_0” onward);<br>
<strong>2)</strong> The pre-miRNA sequence;<br>
<strong>3)</strong> One of the miRNAs sequence;<br>
<strong>4)</strong> The other miRNA sequence (optional);</p>
<p>

One example input file with labels is:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/premirnaplot/imgs/ex1.png?raw=true" alt="enter image description here"></p>
<p>One example input with no labels:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/premirnaplot/imgs/ex2.png?raw=true" alt="enter image description here"></p>
<p>Then, you can run premirnaplot by typing:</p>
<pre><code>premirnaplot file1.txt
</code></pre>
<p>Or with multiple files (also defining the output directory name with the -o parameter):</p>
<pre><code>premirnaplot file1.txt file2.txt -o analysis
</code></pre>
<h3 id="styles">Image styles</h3>
<p>PremiRNAplot has 5 different styles for creating the image for the precursor. This is how the final images look like and how to pick the style you want:</p>
<pre><code>premirnaplot your_file.txt -s 4
</code></pre>
<img src="https://raw.githubusercontent.com/igrorp/pre-miRNA-plot/6f6ea139fd5d2d29cb5e950293e4feca3e231d70/premirnaplot/imgs/all.svg" class="preimg">
<h3 id="colors">Colors</h3>
<p>You can set which colors will be used to highlight the miRNAs within the precursor. Choose between predefined colors (blue, red, green, purple, pink, yellow, cyan, white, black and orange) or select a particular color tone informing its RGB code.
</p><blockquote> You can get RGB codes from selected colors in this <a href="[Color picker](https://www.w3schools.com/colors/colors_picker.asp)">website</a><p></p></blockquote>If you wanted to set the colors to green and blue, for example, you would have to type:
<pre><code>python3 premirnaplot your_file.txt -c green blue</code></pre>
<p>If you wanted to set the colors to a custom tone of purple and yellow, you could type:
</p><pre><code>python3 premirnaplot your_file.txt -c 204 0 205 255 255 102</code></pre>
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/premirnaplot/imgs/im3.svg" align="left" width="400px">
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/premirnaplot/imgs/im4.svg" width="400px">

<h3 id="labels">Image formats</h3>
<p>Pre-miRNA-plot also allows you to save the final images in SVG or PDF format. Note that SVG is a really great format because it hardly loses quality, although some operational systems/web browsers are not compatible. PDF in the other hand is  pretty much universal but it will take a little longer to generate the images.</p>
<pre><code>premirnaplot your_file.txt -f pdf
</code></pre>
<p>Or the default:</p>
<pre><code>premirnaplot your_file.txt -f svg
</code></pre>
<h3 id="multiprocessing">Multiprocessing</h3>
<p>You can set the number of threads to run and speed up a lot the image generation, as the example below: </p>
<pre><code>premirnaplot your_file.txt -t 8</code></pre>
<h1 id="library-usage">4. Library usage</h1>
<p>After succesfully installing premiRNAplot, you can use the library inside a script or application. The fundamental class you’ll need to use is <strong>Precursor</strong>. You can import it using:</p>
<pre class="  language-python"><code class="prism  language-python"><span class="token keyword">from</span> premirnalot<span class="token punctuation">.</span>precursor <span class="token keyword">import</span> Precursor
</code></pre>
<p>And that’s all there is to it. All the methods available for creating the images, extracting and using features or loading an input file are described in the sections below.</p>
<h2 id="basic-properties-and-methods">4.1 Basic properties and methods</h2>
<h2 id="features">4.2 Features</h2>
<p>These are the current features implemented in premiRNAplot. They are divide in Structure, Thermodynamics, Structural-Thermodynamics and miRNA related features.</p>

### Sequence related features


|               **Feature name**              | **Attribute/method** | **Data type** | **Data size** |                      **Description**                     | **Reference** |
|:---------------------------------------:|:----------------:|:---------:|:---------:|:----------------------------------------------------:|:---------:|
|             Sequence length             |    prec.seqlen   |   float   |     1     |                Length of the sequence                |           |
|          Nucleotide frequencies         |   prec.freqs()   |    dict   |     4     |       Frequency of nucleotides (%A, %T, %C, %G)      |           |
|         Dinucleotide frequencies        |  prec.difreqs()  |    dict   |     16    | Frequency of dinucleotides (%AA, %AT, %AC, %AG, ...) |           |
|               G+C content               |                  |           |           |                                                      |           |
|                G/C ratio                |                  |           |           |                                                      |           |
|     Thermodynamics related features     |                  |           |           |                                                      |           |
|               Feature name              | Attribute/method | Data type | Data size |                      Description                     | Reference |
|                 Triplets                |    prec.seqlen   |   float   |     1     |                Length of the sequence                |           |
|             Huang's elements            |   prec.freqs()   |    dict   |     4     |       Frequency of nucleotides (%A, %T, %C, %G)      |           |
|               Stem number               |  prec.difreqs()  |    dict   |     16    | Frequency of dinucleotides (%AA, %AT, %AC, %AG, ...) |           |
|         Base pair type frequency        |                  |           |           |                                                      |           |
|    Base pair type frequency per stem    |                  |           |           |                                                      |           |
|     Average base pair type per stem     |                  |           |           |                                                      |           |
|           Average stem length           |                  |           |           |                                                      |           |
|        Total nucleotides in stem        |                  |           |           |                                                      |           |
|               Longest stem              |                  |           |           |                                                      |           |
|           Terminal loop length          |                  |           |           |                                                      |           |
|         Terminal loop GC content        |                  |           |           |                                                      |           |
|              Bulges number              |                  |           |           |                                                      |           |
|               Loops number              |                  |           |           |                                                      |           |
|               Longest loop              |                  |           |           |                                                      |           |
|         Asymmetric loops number         |                  |           |           |                                                      |           |
|          Symmetric loops number         |                  |           |           |                                                      |           |
|  Average nucleotides in symmetric loops |                  |           |           |                                                      |           |
| Average nucleotides in asymmetric loops |                  |           |           |                                                      |           |
|      Average nucleotides in bulges      |                  |           |           |                                                      |           |
|            Average bulge size           |                  |           |           |                                                      |           |
|           Number of base pairs          |                  |           |           |                                                      |           |
|       Normalized base pair number       |                  |           |           |                                                      |           |

<!--stackedit_data:
eyJoaXN0b3J5IjpbNDI4MTQzMCwtMTg1NDkwNjc3MV19
-->