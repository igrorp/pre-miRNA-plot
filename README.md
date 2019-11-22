<h1 align="center" id="welcome-to-pre-mirna-plot-manual">Welcome to pre-miRNA-plot manual!</h1>
<p></p>
<p>Pre-miRNA-plot is a program for generating multiple custom images of miRNA precursors based on RNAfold and RNAplot. It allows you to highlight the miRNA location within the precursor and obtain general and practical information about your data, so you can filter it or use it in publications. You can see the information in this tutorial in a more visual way in the <a href="https://github.com/igrorp/pre-miRNA-plot/blob/mastee.pdf">documentation</a> file.</p>
<ol>
<li><a href="#1-configuration">Configuration</a><br>
1.1 <a href="#11-python">Python</a><br>
1.2 <a href="#12-vienna-rna-package">Vienna RNA package</a></li>
<li><a href="#2-installation">Installation</a></li>
<li><a href="#3-how-to-use">How to use</a><br>
3.1 <a href="#31-input-files">Input files</a><br>
3.2 <a href="#32-exploring-the-parameters">Exploring the parameters</a></li>

</ol>
<h1 id="configuration">1. Configuration</h1>
Pre-miRNA-plot has some dependencies and you need to check whether you have to install or update them.

<h2 id="python">1.1 Python</h2>
<p>Pre-miRNA-plot runs in Python3+. You can check which version of Python you have installed in your machine with the command bellow:</p>
<pre><code>python --version</code></pre>
<p>Anything higher than 3.0 should work just fine. In  case you have an older version, you can go to the <a href="https://www.python.org/downloads/"> Python website</a> and follow their tutorial to update the platform to a more recent release.</p>

<h2 id="vienna-rna-package">1.2 Vienna RNA package</h2>
<p>Vienna RNA package contains <strong>RNAfold</strong> and <strong>RNAplot</strong>, used to predict the secondary structure of the pre-miRNA and to generate an SVG image of it, respectively. You can visit their <a href="https://www.tbi.univie.ac.at/RNA/documentation.html">website</a> to have more details and get to know more about this amazing package. If you don't have it installed you can download the package and compile it following the steps below:</p>
<pre><code>wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz
tar xvzf ViennaRNA-2.4.13.tar.gz
cd ViennaRNA-2.4.13/
./configure
make
sudo make install
</code></pre>

<h1 id="installation">2. Installation</h1>
<p>To install pre-miRNA-plot you can download this repository as a zip file in this main page, or clone it in your machine:</p>

<pre><code>git clone https://github.com/igrorp/pre-miRNA-plot.git</code></pre>

<h1 id="how-to-use">3. How to use</h1>
<p>Pre-miRNA-plot usage is very simple. The file input has to be a TSV (tab-separeted values) text file containing first the pre-miRNA and then the miRNA sequences. You can specify the colors (in RGB code) to highlight the miRNAs, quality values and other parameters. If you’re not familiar with the files and how to set parameters, the next sessions will explore this properties more. You can see some the information about the parameters of the program by typing <code>premirnaplot -help</code>.</p>

<h2 id="input-files">3.1 Input files</h2>
<p>The input files are text files containing the pre-miRNA sequence and the miRNA sequences, separated by <strong>tabs</strong>. They should look something like this:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/ex1.png" alt="Example 1"></p>
<blockquote>
<p>Note that you don’t need to have necessarly both miRNA sequences, you can have just one of them.</p></blockquote>
<p>If you have labels or some sort of  to your data, you can include them in the first column, like this:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/ex2.png" alt="Example 2"></p>
<blockquote>
<p>The created image files will be named accordingly to the label, as it is a more efficient way of organizing your data</p></blockquote>

<h2 id="exploring-the-parameters">3.2 Exploring the parameters</h2>

<h3 id="labels">Labels</h3>
<p>If you included labels in your input files, as described above, make sure you set the <code>-a</code> parameter to True, like this:</p>
<pre><code>premirnaplot your_file.txt -a T</code></pre>

<h3 id="labels">Styles</h3>

Pre-miRNA-plot has 5 different styles for creating the image for the precursor. This is how the final images look like and how to pick the style you want. 

    python premirnaplot.py your_file.txt -s 4

.preimg img {
	height: 350px;
	margin: 

<div class="preimg">
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/style1.svg"/>
 <img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/style2.svg"/>
  <img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/style3.svg"/> 
  <img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/style4.svg"0/> 
  <img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/style5.svg"/>
</div>

<h3 id="colors">Colors</h3>

<p>You can set which colors will be used to highlight the miRNAs within the precursor. Choose between predefined colors (blue, red, green, purple, pink, yellow, cyan, white, black and orange) or select a particular color tone informing its RGB code.
</p><blockquote> You can get RGB codes from selected colors in this <a href="[https://www.w3schools.com/colors/colors_picker.asp](https://www.w3schools.com/colors/colors_picker.asp)">website</a><p></p></blockquote>If you wanted to set the colors to green and blue, for example, you would have to type:

<pre><code>premirnaplot your_file.txt -c green blue</code></pre>

<p>If you wanted to set the colors to a custom tone of purple and pink, you could type:
</p><pre><code>premirnaplot your_file.txt -c 204 0 205 255 51 153</code></pre>

<p></p><p align="center">
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/im3.svg" height=600>
</p>

<h3 id="labels">Image formats</h3>

Pre-miRNA-plot also allows you to save the final images in SVG or PDF format. Note that the SVG is a really great format because it hardly loses quality, although some operational systems/web browsers are not compatible. PDF in the other hand is  pretty much universal but it will take a little longer to generate the images.

    python premirnaplot.py your_file.txt -f pdf

Or the default:

    python premirnaplot.py your_file.txt -f svg

<h3 id="multiprocessing">Multiprocessing</h3>
<p>You can set the number of allowed processors to run and speed up a lot the program execution, as the example below: </p>
<pre><code>premirnaplot your_file.txt -t 8</code></pre>


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTg1Mzk4NDA4NiwtMTY4MDU5OTgwNywtOD
A4OTMwNjEzLDEwNDE4MTQ4ODAsNjk0MTE2MjAyLC0xNjgyNTgy
NzQ4LC00NjcxMjUxNjcsLTE2Mzc3MzcxMTIsMTg2MjUyODEzNy
wtMjAzODY0NDI5MSw4NjgxMTIzNzcsLTEzNzcyNDA5OTAsMzIx
MTgxNTcsLTEyNDcxMzI0MjIsMTMyMTAwNDkyNywtMjAxMzEwNz
gyNF19
-->