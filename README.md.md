\-\-\- \-\-\-

# Welcome to pre-miRNA-plot manual!

Pre-miRNA-plot is a program for generating multiple images of miRNA precursors using RNAfold and RNAplot. It allows to highlight the miRNA location within the precursor and obtain general and practical information about your data, so you can filter it or use it in publications. You can see the information in this tutorial in a more visual way in the [documentation](https://github.com/igrorp/pre-miRNA-plot/blob/master/documentation.pdf) file.

1.  [Configuration](#1-configuration)  
    1.1 [Python](#11-python)  
    1.2 [Matplotlib](#12-matplotlib)  
    1.3 [Ghostscript](#13-ghostscript)  
    1.4 [Vienna RNA package](#14-vienna-rna-package)
2.  [Installation](#2-installation)
3.  [How to use](#3-how-to-use)  
    3.1 [Input files](#31-input-files)  
    3.2 [Exploring the parameters](#32-exploring-the-parameters)

# 1\. Configuration

Pre-miRNA-plot has some dependencies and you need to check whether you have to install or update them.

## 1.1 Python

Pre-miRNA-plot runs in Python3+. You can check which version of Python you have installed in your machine with the command bellow:

```
python --version

```

Anything higher than 3.0 should work just fine. In case you have an older version, you can go to the [Python website](https://www.python.org/downloads/) and follow their tutorial to update the platform to a more recent release.

## 1.2 Matplotlib

Matplotlib is a graphing package for Python; pre-miRNA-plot uses it to generate diferente plots. You can install it using the **pip** Python package manager.

```
python –m pip install –U matplotlib

```

If you don’t have pip installed, you can use **apt** instead:

```
sudo apt install python3-matplotlib

```

If you are having trouble, please visit their [website](https://matplotlib.org/3.1.1/users/installing.html) for more details.

## 1.3 Ghostscript

Ghostscript is used to convert the Postscript files (.ps) to Portable Network Graphics images (.png). If you do not have it installed, we have to download the tar ball file containing the program, decompress it and then compile it.

```
wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs927/ghostscript-9.27.tar.gz
tar xvzf ghostscript-9.27.tar.gz
cd ghostscript-9.27/
./configure
make
sudo make install

```

## 1.4 Vienna RNA package

Vienna RNA package contains **RNAfold** and **RNAplot**, used to predict the secondary structure of the pre-miRNA and to “costumize” it, respectively. You can visit their [website](https://www.tbi.univie.ac.at/RNA/documentation.html) to have more details and get to know more about this amazing package.

```
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz
tar xvzf ViennaRNA-2.4.13.tar.gz
cd ViennaRNA-2.4.13/
./configure
make
sudo make install

```

> If you are not familiar with command line, each one of the lines above has to be run separetly. Copy the command in each line, run it and wait to see if the process was successfull.

# 2\. Installation

To install pre-miRNA-plot you can download this repository as a zip file in the main page, or clone it in your machine:

```
git clone https://github.com/igrorp/pre-miRNA-plot.git

```

After decompression or cloning, you have to enter the folder and run the [install.sh](http://install.sh) file to make the program executable and to move it to _/usr/local/bin/_ so you can access it from anywhere. You will need superuser permission for that.

```
./install.sh

```

You can check if the program has been successfully moved and installed by testing it with our test data.

```
premirnaplot test_data/osativa.txt -a T -c x x x x x x

```

# 3\. How to use

Pre-miRNA-plot usage is very simple. The file input has to be a TSV (tab-separeted values) text file containing first the pre-miRNA and then the miRNA sequences. You can specify the colors (in RGB code) to highlight the miRNAs, quality values and other parameters. If you’re not familiar with the files and how to set parameters, the next sessions will explore this properties more. You can see some the information about the parameters of the program by typing `premirnaplot --help`.

## 3.1 Input files

The input files are text containing the pre-miRNA sequence and the miRNA sequences, separated by **tabs**. They should look something like this:

![Example 1](https://github.com/igrorp/pre-miRNA-plot/blob/master/ex1.png)

> Note that you don’t need to have necessarly both miRNA sequences, you can have just one of them.

If you have labels or some sort of annotation to your data, you can include them in the first column, like this:

![Example 2](https://github.com/igrorp/pre-miRNA-plot/blob/master/ex2.png)

> The created image files will be named accordingly to the label, so it is a more efficient way of organizing your data

## 3.2 Exploring the parameters

### Labels


### Colors
You can choose which **colors** will be used to highlight the 5p and 3p miRNAs, respectively. You can use predefined colors’ names (green, black, red, blue, white) or specify the **RGB codes** corresponding to the colors. There is no problem if you inform just one color. By default, the sequences will
be colored red (5p) and green (3p).
>You can choose colors and get their RGB codes in this [website]([https://www.w3schools.com/colors/colors_picker.asp](https://www.w3schools.com/colors/colors_picker.asp))

If you wanted to set the colors to blue and red you would have to type:

    premirnaplot your_file.txt -c blue red
 If you wanted to set the colors to purple and pink, you would have to type:

    premirnaplot your_file.txt -c 204 0 205 255 51 153
<p align="center">
<img src="https://github.com/igrorp/pre-miRNA-plot/blob/master/colors.png" width=314 height=480/>
</p>

### Resolution

You can also choose the **resolution** of the images, the default is 200 pixels. You can increase this value for better images for publications or to be able to zoom in a particular area of the picture and not losing quality, for example, but be aware that this increases the program execution time.

If you wanted to set the resolution of the images to 1200 pixels, you wold have to type:

    premirnaplot your_file.txt -q 1200 


<!--stackedit_data:
eyJoaXN0b3J5IjpbMTUzMzQwNTMyNCwxNTA4MDgzNDQxLDE1Mj
c2Mjk1MDcsMzYyMzY0NDAzLDE3MjE5ODQyNTgsMTQ4NzcxNDM3
MSwxNTA5NTc2NzQ2XX0=
-->