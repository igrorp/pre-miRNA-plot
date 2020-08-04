<h1 align='center'> Welcome to premiRNAplot! </h1>
<p></p>


PremiRNAplot is a Python 3.5+ library with a command line interface for extracting features  and generating multiple custom images of miRNA precursors. It provides over a 100 different features that can be used to train and test machine learning algorithms for the correct prediction and classification of pre-miRNAs.

It's based on RNAfold and RNAplot from Vienna RNA and simplifies feature obtention and calculation. It also allows you to produce beautiful, costumizable and publication-ready images of the secondary structure of these precursors, highlighting the position of the miRNAs sequences within them.

<br>
<ol>
<li> <a href="#1-configuration">Configuration</a><br>
1.1 <a href="#11-python">Python</a><br>
1.2 <a href="#12-vienna-rna-package">Vienna RNA package</a><br>
</li><li>  <a href="#2-installation">Installation</a><br>
</li><li> <a href="#3-how-to-use">How to use</a><br>
3.1 <a href="#31-input-files">Input files</a><br>
3.2 <a href="#32-exploring-the-parameters">Exploring the parameters</a>
<br>
</li></ol>

<h2 id="vienna-rna-package">1.2 Vienna RNA package</h2>
<p>Vienna RNA package contains <strong>RNAfold</strong> and <strong>RNAplot</strong>, used to predict the secondary structure of the pre-miRNA and to generate an SVG image of it, respectively. You can visit their <a href="https://www.tbi.univie.ac.at/RNA/documentation.html">website</a> to have more details and get to know more about this amazing package. If you don't have it installed, you can download the package and compile it following the steps below:</p>
<pre><code>
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz
tar xvzf ViennaRNA-2.4.13.tar.gz
cd ViennaRNA-2.4.13/
./configure
make
sudo make install
</code></pre>
<h2 id="pip">1.3 Pip</h2>
<p>The program dependes on <strong>pip</strong> to install all the Python imaging packages used (matplotlib, svglib, svgwrite, scikit-learn, etc.). Without it the installation would be vary manual, time-consuming and maybe conflicting. If you don’t have it already on your machine, run:</p>
<pre><code>sudo apt update
sudo apt install python3-pip
</code></pre>
<p>After that, check it:</p>
<pre><code>python3 -m pip --version
</code></pre>
<p>Now, install <strong>setuptools</strong> and <strong>wheel</strong>, that will help you download all the required packages:</p>
<pre><code>sudo python3 -m pip install --upgrade pip setuptools wheel
</code></pre>
<h1 id="installation">2. Installation</h1>

<p>There are multiple ways of downloading premiRNAplot. We highly suggest that you use the Docker container or the Conda package, as they would remove the need to manually install packages (ViennaRNA in specific) or having conflicts. Make sure that you have your platform updated and that you are using Python 3.5 or higher.</p>

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
python3 setup.py install
</code></pre>
<p>There is going to be a lot of text being displayed informing the packages downloaded and the operation status. After that, you can run our test data to check if all the requirements are satisfied or if there are any problems during the execution.</p>

<h1 id="command-line-interface-usage">3. Command line interface usage</h1>
<p>Pre-miRNA-plot usage is very simple. The file input has to be a TSV (tab-separeted values) text file containing first the pre-miRNA and then the miRNA sequences. You can specify the colors (in RGB code) to highlight the miRNAs, quality values and other parameters. You can also see this information about the parameters of the program by typing <code>python3 premirnaplot.py -help</code>.</p>
<p>If you are inside the cloned repository folder, you can run the program like this:</p>
<pre><code>python3 premirnaplot.py path/to/your/input_file.txt ... (parameters)
</code></pre>
<p>If you are outside the program repository, and you want to run the program from wherever you are, use:</p>
<pre><code>python3 path/to/where/the/program/is/premirnaplot.py input_file.txt ... (parameters)
</code></pre>
<h2 id="input-files">3.1 Input files</h2>
<p>The input files are text files containing the pre-miRNA sequence and the miRNA sequences, separated by <strong>tabs</strong>. They should look something like this:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/ex2.png" alt="Example 1"></p>
<blockquote>
<p>Note that you don’t need to have necessarly both miRNA sequences, you can have just one of them.</p></blockquote>
<p>If you have labels or some sort of  to your data, you can include them in the first column, like this:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/ex1.png" alt="Example 2"></p>
<blockquote>
<p>The created image files will be named accordingly to the label, as it is a more efficient way of organizing your data</p></blockquote>
<h2 id="exploring-the-parameters">3.2 Exploring the parameters</h2>
<h3 id="labels">Labels</h3>
<p>If you included labels in your input files, as described above, make sure you set the <code>-a</code> parameter to True, like this:</p>
<pre><code>python premirnaplot.py your_file.txt -a T</code></pre>
<h3 id="styles">Styles</h3>
<p>Pre-miRNA-plot has 5 different styles for creating the image for the precursor. This is how the final images look like and how to pick the style you want.</p>
<pre><code>python3 premirnaplot.py your_file.txt -s 4
</code></pre>
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/all.svg" class="preimg">
<h3 id="colors">Colors</h3>
<p>You can set which colors will be used to highlight the miRNAs within the precursor. Choose between predefined colors (blue, red, green, purple, pink, yellow, cyan, white, black and orange) or select a particular color tone informing its RGB code.
</p><blockquote> You can get RGB codes from selected colors in this <a href="[https://www.w3schools.com/colors/colors_picker.asp](https://www.w3schools.com/colors/colors_picker.asp)">website</a><p></p></blockquote>If you wanted to set the colors to green and blue, for example, you would have to type:
<pre><code>python3 premirnaplot your_file.txt -c green blue</code></pre>
<p>If you wanted to set the colors to a custom tone of purple and yellow, you could type:
</p><pre><code>python3 premirnaplot your_file.txt -c 204 0 205 255 255 102</code></pre>
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/im3.svg" align="left" width="400px">
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/im4.svg" width="400px">
<h3 id="labels">Image formats</h3>
<p>Pre-miRNA-plot also allows you to save the final images in SVG or PDF format. Note that SVG is a really great format because it hardly loses quality, although some operational systems/web browsers are not compatible. PDF in the other hand is  pretty much universal but it will take a little longer to generate the images.</p>
<pre><code>python3 premirnaplot.py your_file.txt -f pdf
</code></pre>
<p>Or the default:</p>
<pre><code>python3 premirnaplot.py your_file.txt -f svg
</code></pre>
<h3 id="multiprocessing">Multiprocessing</h3>
<p>You can set the number of allowed processors to run and speed up a lot the program execution, as the example below: </p>
<pre><code>premirnaplot your_file.txt -t 8</code></pre>
<h1 id="library-usage">4. Library usage</h1>
<p>After succesfully installing premiRNAplot, you can use the library inside a script or application. The fundamental class you’ll need to use is <strong>Precursor</strong>. You can import it using:</p>
<pre class=" language-python"><code class="prism  language-python"><span class="token keyword">from</span> premirnalot<span class="token punctuation">.</span>precursor <span class="token keyword">import</span> Precursor
</code></pre>
<p>And that’s all there is to it. All the methods available for creating the images, extracting and using features or loading an input file are described in the sections below.</p>
<h2 id="basic-properties-and-methods">4.1 Basic properties and methods</h2>
<h2 id="features">4.2 Features</h2>
<p>These are the current features implemented in premiRNAplot. They are divide in Structure, Thermodynamics, Structural-Thermodynamics and miRNA related features.</p>

<table>
<thead>
<tr>
<th align="center">Feature name</th>
<th align="center">Class property</th>
<th align="center">Data type</th>
<th align="center">Description</th>
<th align="center">Reference</th>
</tr>
</thead>
<tbody>
<tr>
<td align="center">Length</td>
<td align="center">prec.seqlen</td>
<td align="center">int</td>
<td align="center">The length of the sequence</td>
<td align="center">-</td>
</tr>
<tr>
<td align="center">Triplets</td>
<td align="center">The triplet elements frequencies</td>
<td align="center">dict (size 32)</td>
<td align="center"></td>
<td align="center"></td>
</tr>
</tbody>
</table>
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTYxOTEwMzgzNCwxMjM5MDY3OTQxXX0=
-->