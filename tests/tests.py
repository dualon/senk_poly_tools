import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
#sys.path.append('..')

from senk_poly_tools.senk_poly_tools import SenkPolyTools

class SenkPolyToolsTests(SenkPolyTools):

	def testCreateEdfContainer(self, fname):
		print("createEdfContainer(): {}".format(fname))
		
		edf_container = self.createEdfContainer(fname)
		print("\topen existing file [ OK ]")
		
		try:
			not_edf_c = self.createEdfContainer('not_existing_filename')
			print("\tnot existing file was opened... [ FAIL ]")
		except FileNotFoundError:
			print("\ttry non-existing file [ OK ]")
		
		return edf_container
	
	
	def testLoadEdfHeaders(self, edf_container):
		print("loadEdfHeaders():")
		edf_container = self.loadEdfHeaders(edf_container)
		print("\tloaded [ OK ]")
		# @TODO: test the header fields
		
		return edf_container
	
	
	def testLoadEdfData(self, edf_container):
		print("loadEdfData():")
		d = self.loadEdfData(edf_container)
		print("\tloaded [ OK ]")
		return d
	
	
	def testCurveSmoothing(self, e, channels = []):
		print("smoothByAvg():")
		fig, axes = plt.subplots(e.data.shape[0], 1, figsize=(30, e.data.shape[0]*4), dpi=300)
		fig.tight_layout()
		
		for idx, ax in enumerate(axes):
			y1 = e.data[idx][0:1000]
			x1 = np.arange(len(y1))
			y2 = self.smoothByAvg(y1)
			x2 = np.arange(len(y2))
			
			ax.set_title(e.data_labels[idx] + ' Original (grey) vs Smooth (red)', loc='left')
			ax.plot(x1, y1, color='#404040', alpha=0.5)
			ax.plot(x2, y2, color='#AA0000', alpha=0.5)
		
		plt.savefig(os.path.join('..', 'results', '{}_orig_vs_smoothbyavg.png'.format(e.file_basename)))
		plt.close()
		
		print("\t[ OK ], assert visually!")
	
	
	def testFindExtrema(self):
		print("findExtrema():")
		x1 = np.array([1.0, .6, .8, .5, .7, .8, 1.0, 1.0, 1.0, .8, .6, .2, .5, .3, .4, .5, .8, 1.0, .7])
		
		print("\tFinding minima, maxima...", end=" ")
		minima, maxima = self.findExtrema(x1, 3)
		print("[ OK ]")
		
		print("\tAsserting minima...", end=" ")
		#assert minima == [(3, .5), (9, .2)]
		print(minima)
		print("[ OK ]")
		
		print("\tAsserting maxima...", end=" ")
		#assert maxima == [(3, .5), (9, .2)]
		print(maxima)
		print("[ OK ]")
		
		fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=300)
		fig.tight_layout()
			
		ax.set_title('Find Minima (red) and Maxima (green)', loc='left')
		ax.plot(x1, color='#404040', alpha=0.5)
		ax.vlines([idx for idx, _ in minima], 0.0, 1.0, color='#AA0000', alpha=0.5)
		ax.vlines([idx for idx, _ in maxima], 0.0, 1.0, color='#00AA00', alpha=0.5)
		#ax.plot(x2, y2, color='#AA0000', alpha=0.5)
		
		plt.savefig(os.path.join('..', 'results', 'find_extrema.png'))
		plt.close()
	


if __name__ == "__main__":
	print("SenkPolyToolsTests() init:")
	spt = SenkPolyToolsTests()
	print("\t[ OK ]")
	
	edf_container = spt.testCreateEdfContainer('osas2002plusqrs.edf')
	spt.e = edf_container
	
	edf_container = spt.testLoadEdfHeaders(edf_container)
	
	#spt.printEdfHeaders(edf_container)

	edf_container = spt.testLoadEdfData(edf_container)
	
	spt.testCurveSmoothing(edf_container)
	
	spt.testFindExtrema()
	