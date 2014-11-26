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
import pprint
import datetime
import re
import numpy as np
from math import floor
#import matplotlib.pyplot as plt
import csv

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from edf_container import EdfContainer # senk_poly_tools.

class SenkPolyTools(object):
	
	annot_label = 'EDF Annotations'
	
	def createEdfContainer(self, file_name):
		""" Open an EDF(+) file and return an empty EdfContainer object. """
		
		abs_path = os.path.abspath(os.path.join(
			os.path.dirname(__file__),
			'..',
			'data',
			file_name
		))
		
		if not os.path.isfile(abs_path):
			raise FileNotFoundError("File '{}' not found!".format(abs_path))
		
		e = EdfContainer(abs_path)
		e.file_obj = open(abs_path, 'rb')
		
		return e
	

	def loadEdf(self, file_name):
		e = self.createEdfContainer(file_name)
		e = self.loadEdfHeaders(e)
		e = self.loadEdfData(e)
		
		return e
	
	
	def loadEdfHeaders(self, e):
		""" Load the headers of an EDF(+) file into an EdfContainer.
		
			Existing headers are overwritten.
			
			param e: (empty) EdfContainer object,
				see SenkPolyTools.createEdfContainer()
			returns EdfContainer with loaded headers
		"""
		
		if not isinstance(e, EdfContainer):
			raise Exception("The provided container is not an EdfContainer! Use the container returned by SenkPolyTools.openEdf().")
		
		f = e.file_obj
		
		f.seek(0)
		
		e.version = f.read(8)
		e.local_patient_id = f.read(80)
		e.local_recording_id = f.read(80)
		
		# parse timestamp to standard ISO datetime format
		start_date = f.read(8)
		start_time = f.read(8)
		
		start_date = start_date.decode(encoding='ASCII').split('.')
		start_time = start_time.decode(encoding='ASCII').split('.')
		
		year = int(start_date[2])
		
		e.date_time = datetime.datetime(
			(2000+year if year<=84 else 1900+year),
			int(start_date[1]),
			int(start_date[0]),
			int(start_time[0]),
			int(start_time[1]),
			int(start_time[2])
		)
		
		e.num_of_bytes_in_header = int(f.read(8))
		edf_type 				= f.read(44)[:5]
		e.edf_type 				= edf_type.decode(encoding='ASCII')
		e.num_of_records 		= int(f.read(8))
		e.record_duration 		= float(f.read(8))
		e.num_of_signals 		= int(f.read(4))
		
		nsr = range(e.num_of_signals)
		#print(nsr)
		e.labels 				= [f.read(16).strip().decode(encoding='ASCII') for _ in nsr]
		e.transducer_types 		= [f.read(80).strip().decode(encoding='ASCII') for _ in nsr]
		
		e.physical_dimension 	= [f.read(8).strip().decode(encoding='ASCII') for _ in nsr] # physical_dimensions: uV, cm/s, %, mV, ...
		e.physical_min 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.physical_max 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.digital_min 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.digital_max 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.gain 					= (e.physical_max - e.physical_min) / (e.digital_max - e.digital_min) # numpy arrays
		e.prefiltering 			= [f.read(80).strip() for _ in nsr]
		e.num_of_samples_per_record = [int(f.read(8)) for _ in nsr]
		
		e.sample_freq			= [ns/e.record_duration for ns in e.num_of_samples_per_record] # in Hertz
		
		# reserved bytes for each signal
		f.read(32 * e.num_of_signals)
		
		if f.tell() != e.num_of_bytes_in_header:
			raise ValueError("The number of bytes in the header does not match the file object cursor. Header length mismatch during reading?");
		
		return e
	
	
	def loadEdfData(self, e):
		""" Load the data from an EDF(+) file based on the already loaded headers
			in an EdfContainer object.
		
			param e: EdfContainer object with already loaded headers,
				see SenkPolyTools.loadEdfHeaders()
			returns EdfContainer with loaded data
		"""
		
		if e.file_obj is None:
			raise AttributeError("EdfContainer.file_obj is missing, use SenkPolyTools.createEdfContainer() to create a file stream.")
		
		err_msg = "EdfContainer.{} is missing, call SenkPolyTools.loadEdfHeaders()."
		if e.num_of_records is None:
			raise AttributeError(err_msg.format('num_of_records'))
		
		if e.num_of_samples_per_record is None:
			raise AttributeError(err_msg.format('num_of_samples_per_record'))
		
		if e.labels is None:
			raise AttributeError(err_msg.format('labels'))
		
		if e.digital_min is None:
			raise AttributeError(err_msg.format('digital_min'))
		
		if e.physical_min is None:
			raise AttributeError(err_msg.format('physical_min'))
		
		if e.gain is None:
			raise AttributeError(err_msg.format('gain'))
		
		data = []
		for i in range(e.num_of_records):
			# 'k' is the index in the data list
			# it is equal to j until the annotations channel is reached,
			# k == j-1 after the annotation channel
			k = 0
			for j, curr_num_of_samples in enumerate(e.num_of_samples_per_record):
				buff = e.file_obj.read(curr_num_of_samples * 2)
				
				if len(buff) != curr_num_of_samples * 2:
					raise InputError("Unexpected end of EDF file!")
				
				if e.labels[j] == self.annot_label:
					e.annotations.extend(self.parseEdfAnnotation(buff))
				else:
					# if there is no entry for the channel data, create one
					if len(data) == j:
						data.append([])
					
					# 2-byte little endian integer per channel
					dig = np.fromstring(buff, '<i2').astype(np.float32)
					phys = (dig - e.digital_min[j]) * e.gain[j] + e.physical_min[j]

					# @TODO: consider numpy.r_
					phys = phys.tolist()
					data[k].extend(phys) # note the 'k' index
					
					k = k+1
		
		e.data_labels = [l for l in e.labels if l != self.annot_label]
		
		# @TODO: more elegant numpy conversion
		for chn_idx, chn_data in enumerate(data):
			data[chn_idx] = np.array(chn_data)
		e.data = np.array(data) # e.data is a numpy 2D matrix
		del data
		
		return e
	

	def parseEdfAnnotation(self, annot):
		""" Parse annotations of an EDF+ file.
		
			The EDF+ annotations are unicode encoded bytes. The onset and duration are separated by \x14 (single byte with value 21), these are followed by one or more annotation strings, which are separated by \x15 (single byte with value 20).
			See the EDF+ specification for details:
				http://www.edfplus.info/specs/edfplus.html#edfplusannotations
			See also SenkPolyTools.loadEdfData()
		
			param annot: bytes, the annotations in a record
			returns: list of dicts, each dict is like {'onset': <float>, 'duration': <float>, 'annotation': <string>}
			
			Example
			-------
			>>> ann = SenkPolyTools.parseEdfAnnotation(bytes)
			>>> ann
			... [{'onset': '+0', 'duration': None, 'annotation': ''}, {'onset': '+94.8038', 'duration': '0.0000', 'annotation': '1'}]
		"""
		
		exp = '(?P<onset>[+\-]\d+(?:\.\d*)?)' + \
			'(?:\x15(?P<duration>\d+(?:\.\d*)?))?' + \
			'(\x14(?P<annotation>[^\x00]*))?' + \
			'(?:\x14\x00)'

		return [m.groupdict() for m in re.finditer(exp, annot.decode('utf-8'))]
		

	
	def printEdfHeaders(self, e):
		edf_h = vars(e)
		pprint.pprint(edf_h)
	

	#def visualizeEdf(self, edf, all_minima = {}, all_maxima = {}, channels = [], fig_dpi = 200):
	#	""" Visualize the signals from an EDF file.
	#	
	#		The displayed channels are derived from the data automatically or can be determined explicitly with the 'channels' parameter..
	#	
	#	"""
	#	ax_num = (len(channels) if channels else len(edf.data))
	#	# NOT safe to use explicit channel 0 - TODO: auto-detect data channels
	#	data_len = len(edf.data[0])
	#	x = range(data_len)
	#	
	#	fig, axes = plt.subplots(
	#		ax_num,
	#		1,
	#		figsize=(ax_num*5, 30),
	#		dpi=fig_dpi
	#	)
	#	fig.tight_layout()
	#	
	#	for chn in range(ax_num):
	#		if (not channels or (channels and (chn in channels))):
	#			# the recorded data length may vary per channel
	#			x = range(len(edf.data[chn]))
	#			#if edf.labels[chn] != self.annot_label:
	#			#print("chn: {} = {}, x: {}, y: {}".format(edf.labels[chn], chn, len(x), ))
	#			axes[chn].plot(
	#				x,
	#				edf.data[chn],
	#				color='#000000',
	#				alpha=0.5
	#			)
	#			
	#			axes[chn].set_title(edf.data_labels[chn], loc='left')
	#			
	#			chn_minima = all_minima.get(chn)
	#			if chn_minima:
	#				axes[chn].vlines([idx for idx, _ in chn_minima], -150.0, 150.0, color='#3B209C', alpha=0.5)
	#			
	#			chn_maxima = all_maxima.get(chn)
	#			if chn_maxima:
	#				axes[chn].vlines([idx for idx, _ in chn_maxima], -150.0, 150.0, color='#ED9A00', alpha=0.5)
	#	
	#	plt.savefig(os.path.join('..', 'results', '{}.png'.format(edf.file_basename)))
	#	plt.close()
	
	
	def smoothByAvg(self, x, wl = 3):
		""" Smooth a curve with moving average.
		
			A series of overlapping, wl length sections are created from the data (these are the moving windows), and they are convoluted against ones.
			
			Smoothing by moving average is a basic smoothing function for non-periodic functions (as TCD data is typically not periodic, for example; see FFT for periodic data).
			
			param x numpy array containing the data (typically a channel)
			param wl int, the length of the moving window
			returns a numpy array of the smoothed data
		"""
		if len(x) < wl:
			raise Exception("SenkPolyTools.smoothByAvg: data length < wl")
		
		s = np.r_[ x[wl-1:0:-1], x, x[-1:-wl:-1] ]
		w = np.ones(wl,'d')
		c = np.convolve(w/w.sum(), s, mode='valid')
		
		# the smoothed array will have the length of (the original array + wl -1),
		# but the first wl-1 points are not averaged
		# -> cut wl/2 points at the start of the results
		return c[floor(wl/2):]
	
		
	
	def findExtrema(self, x, hwl = 50):
		""" Find the local extrema (minima and maxima) in a vector within a running window.
		
			The method takes a running window and finds one local minimum and local maximum in it. The window serves as a combined lookback and lookahead frame, the currently examined value (x[i]) is always at the middle between the two frames. The window itself can be different length:
			1. It starts at index 0, its length is hwl (==only the lookahead frame),
				|[*---]-------|
			2. as it moves, more and more values can be included in the lookback frame,
				|[--*---]-----|	
			3. when the current index is equal or greater than the half window length (i>=hwl), but there are more indices remaining than the lookahead frame (len(x)>i+hwl), it's a full window,
				|-[---*---]---|
			4. the lookahead frame starts to shrink as the end of the data is approached,
				|------[---*-]|
			5. and finally it stops when the last value is examined
				|-------[---*]|
			
			If the currently examined value is the smallest/largest in the window,
			it is added to the output, if not, the window moves on.
			
			The result is a tuple of 2 lists, each list containing (index, value) tuples. Every time a local minimum or maximum was found, it is added to the respective list with its index.
			
			The first value is used from plateaus.
			
			Note that scipy.signal.argrelextrema() does *not* always find the smallest local minima, therefore a custom method was created.
			
			Parameters
			----------
			param x: 1D numpy array of data
			param hwl: integer, the half length of the running window
			returns: a (minima, maxima) tuple, both minima and maxima are lists of (index, value) tuples
			
			Example
			-------
			>>> x
			... array([1.0, .6, .8, .5, .7, .8, 1.0, 1.0, 1.0, .8, .6, .2, .5, .3, .4, .5, .8, 1.0, .7])
			>>> minima, maxima = findExtrema(x, 3)
			>>> minima
			... [(3, .5), (11, .2)]
			>>> maxima
			... [(0, 1.0), (6, 1.0), (17, 1.0)]
		"""
		x_len = x.shape[0]

		if x_len < hwl:
			raise Exception("SenkPolyTools.findExtrema: data length < hwl")
		
		minima = []
		maxima = []
		for i in range(x_len):
			ws = i-hwl if i-hwl>=0 else 0 # running window start
			we = i+hwl if i+hwl<=x_len else x_len # running window end

			w = x[ws:we]
			wi = i-ws
			
			wmin_i = np.argmin(w)
			if wi == wmin_i:
				minima.append( (i, x[i]) )
			
			wmax_i = np.argmax(w)
			if wi == wmax_i:
				maxima.append( (i, x[i]) )
		
		return (minima, maxima)
	
	
	def indexDists(self, ind):
		""" Get the arithmetic distance among indices in a list.
		
			Example
			-------
			Find out the RR distances from an ECG channel (let's say that channel 0 is the ECG):
			>>> minima, maxima = SenkPolyTools.findExtrema(EdfContainer.data[0])
			>>> rr_indices, rr_amplitudes = zip(*maxima)
			>>> dists = SenkPolyTools.indexDists(rr_indices)
			>>> heart_rate = [60/(dist/EdfContainer.sample_freq[0]) for dist in dists]
		"""
		
		prev_idx_val = 0
		dist = []
		for idx in ind:
			dist.append(idx-prev_idx_val)
			prev_idx_val = idx
		
		return dist


	def downsample(self, x, orig_freq, freq):
		""" Downsample signal data.
		
			The downsampling ratio is calculated as the floor(original frequency/new freqency), both frequency values should be in the same dimension (preferably Hertz).
			Note that this is a simple downsampling, not decimation, therefore high frequency data may cause aliasing.
			
			param x: numpy array or list, the signal data
			param orig_freq: float, the original frequency
			param freq: float, the new frequency
			returns: a numpy array of the downsampled data
		"""
		stepping = floor(orig_freq/freq)
		
		return np.array(x[:stepping:])
	
	



if __name__ == '__main__':
	argp = argparse.ArgumentParser(description="Build the global ComPPI network and export various subnetworks of it.")
	argp.add_argument(
		'-edf',
		help="The name of the input EDF file (must be in the ./data folder)")
	args = argp.parse_args()
	
	spt = SenkPolyTools()
	edfc = spt.loadEdf(args.edf)
	
	channels = ['TCD  1', 'TCD 2', 'Tonometry', 'CO2', 'EKG1', 'EKG2']
	results_base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results'))
	
	#spt.visualizeEdf(edfc)
	
	print("Annotations... ", end='')
	annot_fname = os.path.join(results_base_path, "{}_annotations.txt".format(edfc.file_basename))
	with open(annot_fname, 'w', newline='') as annfp:
		csvw = csv.writer(annfp)
		
		csvw.writerow(['Onset', 'Duration', 'Annotation(s)'])
		for ann in edfc.annotations:
			if ann['annotation']:
				csvw.writerow([
					ann['onset'],
					ann['duration'],
					ann['annotation']
				])
	print("OK")
	
	for chn_i, chn_data in enumerate(edfc.data):
		chn_n = edfc.data_labels[chn_i]
		chn_data = edfc.data[chn_i]
		print("Channel '{}'... ".format(chn_n), end="")
		
		if 'EKG' in chn_n or 'ECG' in chn_n:
			ecg_freq = edfc.sample_freq[chn_i]
			
			sm_data = spt.smoothByAvg(chn_data, 8)
			minima, maxima = spt.findExtrema(sm_data, 150)
			ecg_max_ind, __ = zip(*maxima)
			ecg_dists = spt.indexDists(ecg_max_ind)
			hr = [60/(dist/ecg_freq) for dist in ecg_dists]
			times = [curr_ecg_idx/ecg_freq for curr_ecg_idx in ecg_max_ind]
			times.insert(0, 0.0)
			
			fname = os.path.join(results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp)
				
				csvw.writerow(['Time', 'Heart Rate'])
				for ecg_t, ecg_hr in zip(times, hr):
					csvw.writerow([ecg_t, ecg_hr])
			
		
			#fig, ax = plt.subplots(1, 1, figsize=(80, 8), dpi=200)
			#fig.tight_layout()
			#
			##x = range(len(chn_data))
			#ax.plot(chn_data[:10000], color='#000000', alpha=0.3 )
			#ax.plot(sm_data[:10000], color='#990000', alpha=0.7 )
			#ax.vlines([idx for idx, _ in minima[:25]], -150.0, 150.0, color='#3B209C', alpha=0.5)
			#ax.vlines([idx for idx, _ in maxima[:25]], -150.0, 150.0, color='#ED9A00', alpha=0.5)
			#
			#plt.savefig(os.path.join('..', 'results', '{}--{}.png'.format(edfc.file_basename, chn_n)))
			#plt.close()
			
			
			print("OK")
		
		elif 'TCD' in chn_n:
			sm_data = spt.smoothByAvg(chn_data)
			sm_data_ds = spt.downsample(sm_data, edfc.sample_freq[chn_i], 0.5)
			
			fname = os.path.join(results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp)
	
				for d in sm_data_ds:
					csvw.writerow([d])
			
			#fig, ax = plt.subplots(1, 1, figsize=(80, 8), dpi=200)
			#fig.tight_layout()
			#
			#ax.plot(chn_data[:10000], color='#000000', alpha=0.3 )
			#ax.plot(sm_data[:10000], color='#990000', alpha=0.7 )
			##ax.vlines([idx for idx, _ in minima[:25]], -150.0, 150.0, color='#3B209C', alpha=0.5)
			##ax.vlines([idx for idx, _ in maxima[:25]], -150.0, 150.0, color='#ED9A00', alpha=0.5)
			#
			#plt.savefig(os.path.join('..', 'results', '{}--{}.png'.format(edfc.file_basename, chn_n)))
			#plt.close()
			
			print("OK")
			
		elif 'Tonometry' in chn_n: # ABP
			sm_data = spt.smoothByAvg(chn_data)
			sm_data_ds = spt.downsample(sm_data, edfc.sample_freq[chn_i], 0.5)
			
			fname = os.path.join(results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp)
				
				for d in sm_data_ds:
					csvw.writerow([d])
			
			#fig, ax = plt.subplots(1, 1, figsize=(80, 8), dpi=200)
			#fig.tight_layout()
			#
			#ax.plot(chn_data[:10000], color='#000000', alpha=0.3 )
			#ax.plot(sm_data[:10000], color='#990000', alpha=0.7 )
			##ax.vlines([idx for idx, _ in minima[:25]], -150.0, 150.0, color='#3B209C', alpha=0.5)
			##ax.vlines([idx for idx, _ in maxima[:25]], -150.0, 150.0, color='#ED9A00', alpha=0.5)
			#
			#plt.savefig(os.path.join('..', 'results', '{}--{}.png'.format(edfc.file_basename, chn_n)))
			#plt.close()
			
			print("OK")
		
		elif 'CO2' in chn_n:
			_, maxima = spt.findExtrema(chn_data)
			__, co2_maxes = zip(*maxima)
			sm_data = spt.smoothByAvg(co2_maxes)
			sm_data_ds = spt.downsample(sm_data, edfc.sample_freq[chn_i], 0.5)
			
			fname = os.path.join(results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp)
				
				for d in sm_data_ds:
					csvw.writerow([d])
			
			#fig, ax = plt.subplots(1, 1, figsize=(80, 8), dpi=200)
			#fig.tight_layout()
			#
			#ax.plot(chn_data[:1000], color='#000000', alpha=0.3 )
			#ax.plot(sm_data[:1000], color='#990000', alpha=0.7 )
			##ax.plot(sm_data_ds[:5], color='#009900', alpha=0.7 )
			#
			#plt.savefig(os.path.join('..', 'results', '{}--{}.png'.format(edfc.file_basename, chn_n)))
			#plt.close()
			
			print("OK")
		
		else:
			print("skipped")