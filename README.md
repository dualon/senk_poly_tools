SenkPolyTools v0.2 Alpha
========================

SenkPolyTools is an EDF(+) parser and analyzer written in Python 3 for the SE Neurology Clinic. It is in early experimental phase.

Requirements
------------

Python 3.3+, Numpy 1.8+, MatPlotLib 1.4+

1. Install Python 3.3+: The official [Python homepage](https://wiki.python.org/moin/BeginnersGuide/Download) has instructions on installing Python on Windows, linux or Mac boxes.

2. Check from the command line if Python 3 has been installed. Linux/Windows/Mac:
	
		$ python3 --version
	
	See these pages on how to use the command line: [Linux](https://help.ubuntu.com/community/UsingTheTerminal), [Windows](http://windows.microsoft.com/en-us/windows-vista/open-a-command-prompt-window), [Mac](http://www.wikihow.com/Get-to-the-Command-Line-on-a-Mac)

3. Install the required dependencies

	Linux:

		$ sudo apt-get install python3-numpy python3-matplotlib
	
	Windows:

	Download and install numpy and matplotlib from Christoph Gohlke's [excellent repository](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

4. Download the source code of SenkPolyTools from this repository.
Beginners may want to use the "Download ZIP" button in the right side panel of this page.


Usage
-----

./ refers to the main directory of the package, e.g. "some_path/senk_poly_tools/".

1. Put your input EDF files in the ./data folder.
2. Create the ./results folder if it doesn't exist.
3. Open the command line and navigate to ./senk_poly_tools (so now you are in some_path/senk_poly_tools/)
4. Run the tool as
	
	$ python3 senk_poly_tools/senk_poly_tools.py -edf your_input_file.edf

	Where "your_input_file.edf" is the input file ("./data/your_input_file.edf")

5. Find the results in the ./results directory.

Command Line Parameters
-----------------------
-s --sampling_frequency:
	Optional parameter to specify the sampling frequency (in Hz) at which the processed/interpolated data will be exported. The value should be a number larger than zero (float or integer). Default sampling frequency = 0.5
	
	$ python3 senk_poly_tools/senk_poly_tools.py -s 10 -edf your_input_file.edf

-r --raw:
	Optional flag to export the raw data with minimal noise reduction. Processed data will also be exported.
	
	$ python3 senk_poly_tools/senk_poly_tools.py -n -edf your_input_file.edf

Notes
-----

The European Data Format (EDF and EDFPlus) specifications are available at:
	[http://www.edfplus.info/index.html](http://www.edfplus.info/index.html)