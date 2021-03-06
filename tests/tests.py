import numpy as np
import matplotlib.pyplot as plt
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
#sys.path.append('..')

from senk_poly_tools.senk_poly_tools import SenkPolyTools

class SenkPolyToolsTests(SenkPolyTools):
	
	colors = [
		'#fff44f', # light yellow
		'#FF8A00', # orange
		'#AA0000', # red
		'#990000', # maroon
		'#11BEFF', # teal
		'#00AA00', # green
		'#813D00', # brown
		'#707070', # grey
		'#1F0064', # dark blue
		'#000000' # black
	]

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
	
	
	def visualizeEdf(self, edf, all_minima = {}, all_maxima = {}, channels = [], fig_dpi = 200):
		""" Visualize the signals from an EDF file.
		
			The displayed channels are derived from the data automatically or can be determined explicitly with the 'channels' parameter..
		
		"""
		ax_num = (len(channels) if channels else len(edf.data))
		# NOT safe to use explicit channel 0 - TODO: auto-detect data channels
		data_len = len(edf.data[0])
		x = range(data_len)
		
		fig, axes = plt.subplots(
			ax_num,
			1,
			figsize=(ax_num*5, 30),
			dpi=fig_dpi
		)
		fig.tight_layout()
		
		for chn in range(ax_num):
			if (not channels or (channels and (chn in channels))):
				# the recorded data length may vary per channel
				x = range(len(edf.data[chn]))
				#if edf.labels[chn] != self.annot_label:
				#print("chn: {} = {}, x: {}, y: {}".format(edf.labels[chn], chn, len(x), ))
				axes[chn].plot(
					x,
					edf.data[chn],
					color='#000000',
					alpha=0.5
				)
				
				axes[chn].set_title(edf.data_labels[chn], loc='left')
				
				chn_minima = all_minima.get(chn)
				if chn_minima:
					axes[chn].vlines([idx for idx, _ in chn_minima], -150.0, 150.0, color='#3B209C', alpha=0.5)
				
				chn_maxima = all_maxima.get(chn)
				if chn_maxima:
					axes[chn].vlines([idx for idx, _ in chn_maxima], -150.0, 150.0, color='#ED9A00', alpha=0.5)
		
		plt.savefig(os.path.join('..', 'results', '{}.png'.format(edf.file_basename)))
		plt.close()
	

	def drawRawCognData(self, edfc, chn_num, windows):
		""" Visualize the cognitive tasks from an EDF channel.
		
			edfc EdfContainer: an EdfContainer object containing the loaded EDF(+) data
			chn_num: the number of channel to extract from the container
			windows: data windows, typically
		"""
		lbl = edfc.data_labels[chn_num]
		data = edfc.data[chn_num]
		print("testCognTasks(): Channel '{}' loaded, data length: {}".format(lbl, len(data)))
		
		colors = [
			'#fff44f', # light yellow
			'#FF8A00', # orange
			'#AA0000', # red
			'#990000', # maroon
			'#11BEFF', # teal
			'#00AA00', # green
			'#813D00', # brown
			'#707070', # grey
			'#1F0064', # dark blue
			'#000000' # black
		]

		#window_matrix = []
		shortest_len = None
		
		f,ax = plt.subplots(1,1,figsize=(12,14))
		for wts, wte, cl in zip(windows, colors):
			w1 = int(wts*500)
			w2 = int(wte*500)
			
			w = data[w1:w2]
			window_matrix.append(w)
			
			#_, co2_maxima_w = spt.findExtrema(w, hwl=2000) # window: 150
			
			dl = len(w)
			
			w1_x = []
			w1_y = []
			for mx_x, mx_y in co2_interp_maxes: # co2_maxima_w
				w1_x.append(mx_x)
				w1_y.append(mx_y)
			
			#print("w1_x len: {}, x: \n{}, y: \n{}".format(len(w1_x), w1_x, w1_y))
			
			if shortest_len is not None:
				shortest_len = min(len(w1_x), shortest_len)
			else:
				shortest_len = len(w1_x)
			
			maximas.append(w1_y)
			
			ax.plot(w1_x, w1_y, color=cl) #, alpha=0.5
			#ax.vlines(w1_x, 20.0, 40.0, color=cl, alpha=0.15)
		
		#print(len(maximas[0]))
		
		maxima_matrix = []
		for r in maximas:
			maxima_matrix.append(r[:shortest_len])
		
		maxima_np_matrix = np.matrix(maxima_matrix)
		print("matrix dim: {}".format(maxima_np_matrix.shape))
		
		w_np_mtx_mean = np.array(maxima_np_matrix.mean(0))
		print("mean dim: {}".format(w_np_mtx_mean.shape))
		
		print(len([x for x in range(0, 170*250, 250)]))
		print(len(w_np_mtx_mean[0]))
		
		#w_np_mtx_mean[0].shape[0]*250
		ax.plot([x for x in range(0, 170*250, 250)], w_np_mtx_mean[0], color='#FF69B4', linewidth=2.0)
		
		plt.savefig(os.path.join('results', '{}_interp_w.png'.format(edfc.file_basename)))
		#plt.savefig(os.path.join('results', '{}_w.png'.format(edfc.file_basename)))
		plt.close()
	

	def drawInterpCognData(self, edfc, chn_num, windows):
		""" Visualize the cognitive tasks from an EDF channel.
		
			edfc EdfContainer: an EdfContainer object containing the loaded EDF(+) data
			chn_num: the number of channel to extract from the container
			windows: data windows, typically
		"""
		lbl = edfc.data_labels[chn_num]
		data = edfc.data[chn_num]
		print("testCognTasks(): Channel '{}' loaded, data length: {}".format(lbl, len(data)))
		


		#window_matrix = []
		maximas = []
		shortest_len = None
		
		f,ax = plt.subplots(1,1,figsize=(12,14))
		for wts, wte, cl in windows:
			w1 = int(wts*500)
			w2 = int(wte*500)
			
			w = D[w1:w2]
			#window_matrix.append(w)
			
			_, co2_maxima_w = spt.findExtrema(w, hwl=2000) # window: 150
			
			dl = len(w)
			#print("w len: {}".format(dl))
			co2_sampling = int(500*0.5) # 500 Hz * 0.5
			interp_x = [ix for ix in range(0, dl, co2_sampling)]
			
			co2_interp_maxes = spt.interpolate(co2_maxima_w, interp_x)
			
			w1_x = []
			w1_y = []
			for mx_x, mx_y in co2_interp_maxes: # co2_maxima_w
				w1_x.append(mx_x)
				w1_y.append(mx_y)
			
			#print("w1_x len: {}, x: \n{}, y: \n{}".format(len(w1_x), w1_x, w1_y))
			
			if shortest_len is not None:
				shortest_len = min(len(w1_x), shortest_len)
			else:
				shortest_len = len(w1_x)
			
			maximas.append(w1_y)
			
			ax.plot(w1_x, w1_y, color=cl) #, alpha=0.5
			#ax.vlines(w1_x, 20.0, 40.0, color=cl, alpha=0.15)
		
		#print(len(maximas[0]))
		
		maxima_matrix = []
		for r in maximas:
			maxima_matrix.append(r[:shortest_len])
		
		maxima_np_matrix = np.matrix(maxima_matrix)
		print("matrix dim: {}".format(maxima_np_matrix.shape))
		
		w_np_mtx_mean = np.array(maxima_np_matrix.mean(0))
		print("mean dim: {}".format(w_np_mtx_mean.shape))
		
		print(len([x for x in range(0, 170*250, 250)]))
		print(len(w_np_mtx_mean[0]))
		
		#w_np_mtx_mean[0].shape[0]*250
		ax.plot([x for x in range(0, 170*250, 250)], w_np_mtx_mean[0], color='#FF69B4', linewidth=2.0)
		
		plt.savefig(os.path.join('results', '{}_interp_w.png'.format(edfc.file_basename)))
		#plt.savefig(os.path.join('results', '{}_w.png'.format(edfc.file_basename)))
		plt.close()


	def drawRrCognData(self, edfc, chn_num, windows, fname = 'rr.png'):

		lbl = edfc.data_labels[chn_num]
		data = edfc.data[chn_num]
		print("drawRrCognData(): Channel '{}' loaded, data length: {}".format(lbl, len(data)))
		
		f,ax = plt.subplots(1,1,figsize=(12,14))
		#ax.set_ylim([30.0, 46.0])
		
		for wts, wte, cl in zip(windows, self.colors):
			w1 = int(wts*500)
			w2 = int(wte*500)
			
			w = data[w1:w2]
			
			_, maxima = spt.findExtrema(w, hwl=2000)
			max_idx, __ = zip(*maxima)
			idx_dists = spt.indexDists(max_idx)
			rr = [60/(dist/500) for dist in idx_dists]
			
			#print("RR: {}".format(rr[:10]))

			ax.plot(rr[1:], color=cl)  # leave out the first value as it is usually an artifact
		
		plt.savefig(os.path.join('results', '{}_{}'.format(edfc.file_basename, fname)))
		plt.close()
		print("OK")
	
	
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
	
	edf_container = spt.testCreateEdfContainer('osas2002plusqrs.edf')
	spt.e = edf_container
	
	edf_container = spt.testLoadEdfHeaders(edf_container)
	
	#spt.printEdfHeaders(edf_container)
	
	edf_container = spt.testLoadEdfData(edf_container)
	
	spt.testCurveSmoothing(edf_container)
	
	spt.testFindExtrema()
	
	print("\tfindExtrema() on each real data channel...", end=" ")
	all_minima = {}
	all_maxima = {}
	for chn_idx, chn_d in enumerate(edf_container.data):
		all_minima[chn_idx], all_maxima[chn_idx] = spt.findExtrema(chn_d, 100)
	print("[ OK ]")
	
	print("visualizeEdf()...", end=" ")
	spt.visualizeEdf(edf_container, all_minima, all_maxima)
	print("[ OK ]")
	
	spt.testSmoothedExtrema(edf_container)
	
	spt.testInterpolate()