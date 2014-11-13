SenkPolyTools v0.1 Alpha
=============

SenkPolyTools is an EDF(+) parser and analyzer for the SE Neurology Clinic. It is in early experimental phase.

Requirements
------------

Python3, Numpy 1.8+, MatPlotLib 1.4+
Put the input EDF files in the ./data folder.
The European Data Format (EDF and EDFPlus) specifications are available at:
	http://www.edfplus.info/index.html

Usage
-----

	$ pwd
	./senk_poly_tools/senk_poly_tools
	$ python3 senk_poly_tools.py -edf filename.edf

SenkPolyTools results are exported to the ./results folder.