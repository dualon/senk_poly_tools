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
		print("\t[ OK ] open existing file")
		
		try:
			not_edf_c = self.createEdfContainer('not_existing_filename')
			print("\t[ FAILED ] not existing file was opened...")
		except FileNotFoundError:
			print("\t[ OK ] try not existing file")
		
		return edf_container
	
	
	def testLoadEdfHeaders(self, edf_container):
		edf_container = self.loadEdfHeaders(edf_container)
		# @TODO: test the header fields
		
		return edf_container
	
	
	def testLoadEdfData(self, edf_container):
		return self.loadEdfData(edf_container)
	
	
	def testCurveSmoothing(self, e, channels = []):

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


if __name__ == "__main__":
	print("SenkPolyToolsTests() init:")
	spt = SenkPolyToolsTests()
	print("\t[ OK ]")
	
	edf_container = spt.testCreateEdfContainer('osas2002plusqrs.edf')
	
	edf_container = spt.testLoadEdfHeaders(edf_container)
	
	spt.printEdfHeaders(edf_container)

	edf_container = spt.testLoadEdfData(edf_container)
	
	spt.testCurveSmoothing(edf_container)
	
	
	