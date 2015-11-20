""" Polygraph Tools, Semmelweis Egyetem Neurologiai Klinika

	The code is an experimental parser and analyzer for the EDF files that
	are results of various polygraph recordings.
	
	author: David M. 'dualon' Gyurko gyurko dot david at e-arc dot hu
	
	The EDF IO was inspired by Boris Reuderink's EEGTools:
	https://github.com/breuderink/eegtools/tree/master/eegtools

"""

import sys
import os
import argparse
import csv
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../senk_poly_tools')))
from senk_poly_tools.senk_poly_tools import SenkPolyTools


argp = argparse.ArgumentParser()
argp.add_argument(
	'-edf',
	help="The name of the input EDF file (must be in the ./data folder)")
args = argp.parse_args()

spt = SenkPolyTools()
edfc = spt.loadEdf(args.edf)

results_base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))

print("Annotations... ", end='')
annot_fname = os.path.join(results_base_path, "{}_annotations.txt".format(edfc.file_basename))
with open(annot_fname, 'w', newline='') as annfp:
	csvw = csv.writer(annfp, dialect='excel', delimiter=';')
	
	csvw.writerow(['Onset', 'Duration', 'Annotation(s)'])
	for ann in edfc.annotations:
		if ann['annotation']:
			csvw.writerow([
				ann['onset'],
				ann['duration'],
				ann['annotation']
			])
print("OK")

print(edfc.data_labels[20])

D = edfc.data[20][130000:200000]
nrD = spt.smoothByAvg(D, 30) # <-- THIS NUMBER IS THE IMPORTANT PARAMETER

fig, ax = plt.subplots(1, 1, figsize=(100, 5), dpi=300)
fig.tight_layout()

x = range(len(D))
nrx = range(len(nrD))

ax.plot(x, D, color='#000000', alpha=0.5)
ax.plot(nrx, nrD, color='#AA0000', alpha=0.5)
ax.set_title(edfc.data_labels[20], loc='left')

plt.savefig(os.path.join('..', 'results', '{}_nrtest.png'.format(edfc.file_basename)))
plt.close()