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
		assert minima == [(3, .5), (11, .2)]
		#print(minima)
		print("[ OK ]")
		
		print("\tAsserting maxima...", end=" ")
		assert maxima == [(0, 1.0), (6, 1.0), (17, 1.0)]
		#print(maxima)
		print("[ OK ]")
		
		fig, ax = plt.subplots(1, 1, figsize=(8, 8), dpi=300)
		fig.tight_layout()
			
		ax.set_title('Find Minima (purple) and Maxima (yellow)', loc='left')
		#ax.set_axis_bgcolor('#a0a0a0')
		ax.plot(x1, color='#404040', alpha=0.5)
		ax.vlines([idx for idx, _ in minima], 0.0, 1.0, color='#3B209C', alpha=0.5)
		ax.vlines([idx for idx, _ in maxima], 0.0, 1.0, color='#ED9A00', alpha=0.5)
		#ax.plot(x2, y2, color='#AA0000', alpha=0.5)
		
		plt.savefig(os.path.join('..', 'results', 'find_extrema.png'))
		plt.close()
	
	
	def testSmoothedExtrema(self, e):
		print("visualize smoothed first channel of real data (smoothing window: 3, findExtrema half-window 100)...", end=" ")
		
		chn_0_sm = spt.smoothByAvg(edf_container.data[0][:2000], 8)
		chn_0_mins, chn_0_maxes = spt.findExtrema(chn_0_sm, 40)
		
		fig, ax = plt.subplots(1, 1, figsize=(100, 8), dpi=200)
		fig.tight_layout()
		
		ax.set_title('Orig. (grey) and smoothed (green) data, find Minima (purple) and Maxima (yellow) on smoothed', loc='left')
		#ax.set_axis_bgcolor('#a0a0a0')
		ax.plot(edf_container.data[0][:2000], color='#000000', alpha=0.3)
		ax.plot(chn_0_sm, color='#008000', alpha=0.5)
		
		ax.vlines([idx for idx, _ in chn_0_mins], -100.0, 150.0, color='#3B209C', alpha=0.5)
		ax.vlines([idx for idx, _ in chn_0_maxes], -100.0, 150.0, color='#ED9A00', alpha=0.5)
		
		plt.savefig(os.path.join('..', 'results', 'osas2002_chn0_find_extrema.png'))
		plt.close()
		print("[ OK ]")
	
	
	def testInterpolate(self):
		print("interpolate():")
		timed_y = [(1,5), (4,8), (8,4)]
		sample_x = [2,3,5,6,7,8]
		
		i_y = self.interpolate(timed_y, sample_x)
		
		print(i_y)
		
		assert i_y[0] == (2,6.0)
		assert i_y[1] == (3,7.0)
		assert i_y[3] == (6,6.0)
		
		print("[ OK ]")
	


if __name__ == "__main__":
	print("SenkPolyToolsTests() init:")
	spt = SenkPolyToolsTests()
	print("\t[ OK ]")
	
	#edf_container = spt.testCreateEdfContainer('osas2002plusqrs.edf')
	#spt.e = edf_container
	#
	#edf_container = spt.testLoadEdfHeaders(edf_container)
	#
	##spt.printEdfHeaders(edf_container)
	#
	#edf_container = spt.testLoadEdfData(edf_container)
	#
	#spt.testCurveSmoothing(edf_container)
	#
	#spt.testFindExtrema()
	#
	#print("\tfindExtrema() on each real data channel...", end=" ")
	#all_minima = {}
	#all_maxima = {}
	#for chn_idx, chn_d in enumerate(edf_container.data):
	#	all_minima[chn_idx], all_maxima[chn_idx] = spt.findExtrema(chn_d, 100)
	#print("[ OK ]")
	#
	#print("visualizeEdf()...", end=" ")
	#spt.visualizeEdf(edf_container, all_minima, all_maxima)
	#print("[ OK ]")
	#
	#spt.testSmoothedExtrema(edf_container)
	
	spt.testInterpolate()