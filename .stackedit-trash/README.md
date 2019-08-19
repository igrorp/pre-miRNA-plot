---


---

<hr>
<hr>
<h1 id="welcome-to-pre-mirna-plot-manual">Welcome to pre-miRNA-plot manual!</h1>
<p>Pre-miRNA-plot is a program for generating multiple images of miRNA precursors using RNAfold and RNAplot. It allows to highlight the miRNA location within the precursor and obtain general and practical information about your data, so you can filter it or use it in publications. You can see the information in this tutorial in a more visual way in the <a href="https://github.com/igrorp/pre-miRNA-plot/blob/master/documentation.pdf">documentation</a> file.</p>
<ol>
<li><a href="#1-configuration">Configuration</a><br>
1.1 <a href="#11-python">Python</a><br>
1.2 <a href="#12-matplotlib">Matplotlib</a><br>
1.3 <a href="#13-ghostscript">Ghostscript</a><br>
1.4 <a href="#14-vienna-rna-package">Vienna RNA package</a></li>
<li><a href="#2-installation">Installation</a></li>
<li><a href="#3-how-to-use">How to use</a><br>
3.1 <a href="#31-input-files">Input files</a><br>
3.2 <a href="#32-exploring-the-parameters">Exploring the parameters</a></li>
</ol>
<h1 id="configuration">1. Configuration</h1>
<p>Pre-miRNA-plot has some dependencies and you need to check whether you have to install or update them.</p>
<h2 id="python">1.1 Python</h2>
<p>Pre-miRNA-plot runs in Python3+. You can check which version of Python you have installed in your machine with the command bellow:</p>
<pre><code>python --version
</code></pre>
<p>Anything higher than 3.0 should work just fine. In  case you have an older version, you can go to the <a href="https://www.python.org/downloads/">Python website</a> and follow their tutorial to update the platform to a more recent release.</p>
<h2 id="matplotlib">1.2 Matplotlib</h2>
<p>Matplotlib is a graphing package for Python; pre-miRNA-plot uses it to generate diferente plots. You can install it using the <strong>pip</strong> Python package manager.</p>
<pre><code>python –m pip install –U matplotlib
</code></pre>
<p>If you don’t have pip installed, you can use <strong>apt</strong> instead:</p>
<pre><code>sudo apt install python3-matplotlib
</code></pre>
<p>If you are having trouble, please visit their <a href="https://matplotlib.org/3.1.1/users/installing.html">website</a> for more details.</p>
<h2 id="ghostscript">1.3 Ghostscript</h2>
<p>Ghostscript is used to convert the Postscript files (.ps) to Portable Network Graphics images (.png). If you do not have it installed, we have to download the tar ball file containing the program, decompress it and then compile it.</p>
<pre><code>wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs927/ghostscript-9.27.tar.gz
tar xvzf ghostscript-9.27.tar.gz
cd ghostscript-9.27/
./configure
make
sudo make install
</code></pre>
<h2 id="vienna-rna-package">1.4 Vienna RNA package</h2>
<p>Vienna RNA package contains <strong>RNAfold</strong> and <strong>RNAplot</strong>, used to predict the secondary structure of the pre-miRNA and to “costumize” it, respectively. You can visit their <a href="https://www.tbi.univie.ac.at/RNA/documentation.html">website</a> to have more details and get to know more about this amazing package.</p>
<pre><code>wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz
tar xvzf ViennaRNA-2.4.13.tar.gz
cd ViennaRNA-2.4.13/
./configure
make
sudo make install
</code></pre>
<blockquote>
<p>If you are not familiar with command line, each one of the lines above has to be run separetly. Copy the command in each line, run it and wait to see if the process was successfull.</p>
</blockquote>
<h1 id="installation">2. Installation</h1>
<p>To install pre-miRNA-plot you can download this repository as a zip file in the main page, or clone it in your machine:</p>
<pre><code>git clone https://github.com/igrorp/pre-miRNA-plot.git
</code></pre>
<p>After decompression or cloning, you have to enter the folder and run the <a href="http://install.sh">install.sh</a> file to make the program executable and to move it to <em>/usr/local/bin/</em> so you can access it from anywhere. You will need superuser permission for that.</p>
<pre><code>./install.sh
</code></pre>
<p>You can check if the program has been successfully moved and installed by testing it with our test data.</p>
<pre><code>premirnaplot test_data/osativa.txt -a T -c x x x x x x
</code></pre>
<h1 id="how-to-use">3. How to use</h1>
<p>Pre-miRNA-plot usage is very simple. The file input has to be a TSV (tab-separeted values) text file containing first the pre-miRNA and then the miRNA sequences. You can specify the colors (in RGB code) to highlight the miRNAs, quality values and other parameters. If you’re not familiar with the files and how to set parameters, the next sessions will explore this properties more. You can see some the information about the parameters of the program by typing <code>premirnaplot --help</code>.</p>
<h2 id="input-files">3.1 Input files</h2>
<p>The input files are text containing the pre-miRNA sequence and the miRNA sequences, separated by <strong>tabs</strong>. They should look something like this:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/ex1.png" alt="Example 1"></p>
<blockquote>
<p>Note that you don’t need to have necessarly both miRNA sequences, you can have just one of them.</p>
</blockquote>
<p>If you have labels or some sort of annotation to your data, you can include them in the first column, like this:</p>
<p><img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/ex2.png" alt="Example 2"></p>
<blockquote>
<p>The created image files will be named accordingly to the label, so it is a more efficient way of organizing your data</p>
</blockquote>
<h2 id="exploring-the-parameters">3.2 Exploring the parameters</h2>
You can choose which colors will be used to highlight the 5p and
3p miRNAs, respectively. You can use predefined colors’ names (green, black, red, blue, white) or
specify the <strong>RGB codes</strong> corresponding to the colors. There is no
problem if you inform just one color. By default, the sequences will
be colored red (5p) and green (3p).


<!--stackedit_data:
eyJoaXN0b3J5IjpbMTgyNjYyMzk2NF19
-->