
<br>
<h1 align="center" id="welcome-to-pre-mirna-plot-manual">Welcome to pre-miRNA-plot!</h1>
<p></p>
<br>
<p>Pre-miRNA-plot is a program for generating multiple custom images of miRNA precursors based on RNAfold and RNAplot. It allows you to highlight the miRNA location within the precursor and obtain general and practical information about your data, so you can filter it or use it in publications.</p>

<br>

<ol>
<li> <a href="#1-configuration">Configuration</font></a><br>
1.1 <a href="#11-python">Python</a><br>
1.2 <a href="#12-vienna-rna-package">Vienna RNA package</a><br>
<li>  <a href="#2-installation">Installation</a><br>
<li> <a href="#3-how-to-use">How to use</a><br>
3.1 <a href="#31-input-files">Input files</a><br>
3.2 <a href="#32-exploring-the-parameters">Exploring the parameters</a>
<br>
</ol>

<h1 id="configuration">1. Configuration</h1>
Pre-miRNA-plot has some dependencies and you need to check whether you have to install or update them.

<h2 id="python"> 1.1 Python</h2>
<p>Pre-miRNA-plot runs in Python3+. You can check which version of Python you have installed in your machine with the command bellow:</p>
<pre><code>python --version</code></pre>
<p>Anything higher than 3.0 should work just fine. In  case you have an older version, you can go to the <a href="https://www.python.org/downloads/"> Python website</a> and follow their tutorial to update the platform to a more recent release.</p>

<h2 id="vienna-rna-package">1.2 Vienna RNA package</h2>
<p>Vienna RNA package contains <strong>RNAfold</strong> and <strong>RNAplot</strong>, used to predict the secondary structure of the pre-miRNA and to generate an SVG image of it, respectively. You can visit their <a href="https://www.tbi.univie.ac.at/RNA/documentation.html">website</a> to have more details and get to know more about this amazing package. If you don't have it installed, you can download the package and compile it following the steps below:</p>
<pre><code>wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz
tar xvzf ViennaRNA-2.4.13.tar.gz
cd ViennaRNA-2.4.13/
./configure
make
sudo make install
</code></pre>

## 1.3 Pip

The program dependes on **pip** to install all the Python imaging packages used (matplotlib, svglib, svgwrite, scikit-learn, etc.). Without it the installation would be vary manual, time-consuming and maybe conflicting. If you don't have it already on your machine, run:

    sudo apt update
    sudo apt install python3-pip

After that, check it:

    python3 -m pip --version

Now, install **setuptools** and **wheel**, that will help you download all the required packages:

      
    sudo python3 -m pip install --upgrade pip setuptools wheel


<h1 id="installation">2. Installation</h1>

There are multiple ways of downloading premiRNAplot, and we highly suggest that you use the Docker container or the Conda package, as they would remove the need to manually install packages (Vienna RNA) or having conflicts.

## Using Conda

It's always a good pratice to create virtual environments to isolate your programs and avoid version conflicts, for example. To create a virtual environment for premiRNAplot, run (this step is optional): 

    conda create -n premirnaplot python=3.7 -y
    conda activate premirnaplot

To install, simply run:

    conda install -c igror premirnaplot

Soon the program will be available through the Bioconda channel.

## Using Docker

Future instructions on how to use a Docker container.

## Using pip

You can install the program from pip as well, and it will automatically install all the Python packages but you need to have the Vienna RNA package alrea

<pre><code>git clone https://github.com/igrorp/pre-miRNA-plot.git</code></pre>

Now you need to enter the cloned repository folder and run the **setup.py** configuration file, which will download all the necessary Python packages used in pre-miRNA-plot.

    cd pre-miRNA-plot/
    python3 setup.py install
   
   There is going to be a lot of text being displayed informing the packages downloaded and the operation status. After that, you can run our test data to check if all the requirements are satisfied or if there are any problems during the execution.

    python3 premirnaplot.py src/example.txt -a T

<h1 id="how-to-use">3. How to use</h1>
<p>Pre-miRNA-plot usage is very simple. The file input has to be a TSV (tab-separeted values) text file containing first the pre-miRNA and then the miRNA sequences. You can specify the colors (in RGB code) to highlight the miRNAs, quality values and other parameters. You can also see this information about the parameters of the program by typing <code>python3 premirnaplot.py -help</code>.</p>

If you are inside the cloned repository folder, you can run the program like this:

    python3 premirnaplot.py path/to/your/input_file.txt ... (parameters)

If you are outside the program repository, and you want to run the program from wherever you are, use:

    python3 path/to/where/the/program/is/premirnaplot.py input_file.txt ... (parameters)



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

Pre-miRNA-plot has 5 different styles for creating the image for the precursor. This is how the final images look like and how to pick the style you want. 

    python3 premirnaplot.py your_file.txt -s 4


<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/all.svg" class="preimg"/>

<h3 id="colors">Colors</h3>

<p>You can set which colors will be used to highlight the miRNAs within the precursor. Choose between predefined colors (blue, red, green, purple, pink, yellow, cyan, white, black and orange) or select a particular color tone informing its RGB code.
</p><blockquote> You can get RGB codes from selected colors in this <a href="[https://www.w3schools.com/colors/colors_picker.asp](https://www.w3schools.com/colors/colors_picker.asp)">website</a><p></p></blockquote>If you wanted to set the colors to green and blue, for example, you would have to type:

<pre><code>python3 premirnaplot your_file.txt -c green blue</code></pre>

<p>If you wanted to set the colors to a custom tone of purple and yellow, you could type:
</p><pre><code>python3 premirnaplot your_file.txt -c 204 0 205 255 255 102</code></pre>


<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/im3.svg" align="left" width="400px">
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/src/im4.svg" width="400px">


<h3 id="labels">Image formats</h3>

Pre-miRNA-plot also allows you to save the final images in SVG or PDF format. Note that SVG is a really great format because it hardly loses quality, although some operational systems/web browsers are not compatible. PDF in the other hand is  pretty much universal but it will take a little longer to generate the images.

    python3 premirnaplot.py your_file.txt -f pdf

Or the default:

    python3 premirnaplot.py your_file.txt -f svg

<h3 id="multiprocessing">Multiprocessing</h3>
<p>You can set the number of allowed processors to run and speed up a lot the program execution, as the example below: </p>
<pre><code>premirnaplot your_file.txt -t 8</code></pre>

<!--stackedit_data:
eyJoaXN0b3J5IjpbMTkzMDA3NTA3LDE1NDgzMjg0ODEsLTY0MD
I0MzA3OSwxMTgzNjQ1NDYwLDY0OTc4MzQ3Nyw2OTc3NTEwNjQs
MTAwNTM3MDg3MSwzNDY3ODczNjYsLTk2NDAyOTUwLDE4MDM4Mz
MyNDYsMTY5ODExNzIwNiw5NzEwNTA5NTUsMjQ0MTYwNDEzLDEx
OTk3MDg5NTcsMTk3Njc2MTI2NiwtMTY5MDcxOTU4LC0zNTI3Mj
A2NjYsODAwNjE5MjkwLDg4NjM1MTQ0MywxMDE5MjkwN119
-->