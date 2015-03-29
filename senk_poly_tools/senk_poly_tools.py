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
import csv

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from edf_container import EdfContainer # senk_poly_tools.

class SenkPolyTools(object):
	
	
	def __init__(self, sampling_freq = 0.5, export_raw = False, annot_label = 'EDF Annotations'):
		""" Sampling frequency (in Hz) determines the sampling rate of the processed data. """
		
		print("SenkPolyTools initiating...", end = " ")
		
		self.results_base_path 	= os.path.abspath(os.path.join(
			os.path.dirname(__file__),
			'..',
			'results'
		))
		self.annot_label		= annot_label
		self.export_nr_raw		= (True if export_raw else False)
		self.sampling_freq		= self.filterSamplingFreq(sampling_freq)
		
		print("OK")
	

	def createEdfContainer(self, file_name):
		""" Open an EDF(+) file and return an empty EdfContainer object with the file handler. """
		
		print("Creating EDF(+) container... ", end="")
		
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
		
		print("OK")
		
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
		
		print("Loading EDF(+) headers... ", end = "")
		
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
		
		print("OK")
		
		return e
	
	
	def loadEdfData(self, e):
		""" Load the data from an EDF(+) file based on the already loaded headers
			in an EdfContainer object.
		
			param e: EdfContainer object with already loaded headers,
				see SenkPolyTools.loadEdfHeaders()
			returns EdfContainer with loaded data
		"""
		
		print("Loading EDF(+) data... ", end = "")
		
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
		
		print("OK")
		
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
		#return c[floor(wl/2):]
		return c
	
		
	
	def findExtrema(self, x, hwl = 150):
		""" Find the local extrema (minima and maxima) in a running window of a vector.
		
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
			msg = "SenkPolyTools.findExtrema: data length ({}) < hwl ({})".format(x_len, hwl)
			raise Exception(msg)
		
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
	

	def interpolate(self, timed_y, sample_x):
		""" Interpolate linearly sparsely sampled data to get more frequently sampled output.
		
			The method takes a list of values (channel data) that are sparser than the resampling frequency, and calculates the intermediate y values at each 'sample_x' coordinate. The calculation is a linear interpolation (weighted average).
		
			param timed_y: list of tuples, each tuple is (x_coordinate, y_value)
			param sample_x: list of integers, each integer is an x_coordinate
			returns: list of tuples, (sample_x_i, interpolated_y)
			
			Example
			-------
			Find the systolic blood pressure values (the maxima of the arterial blood pressure) in data sampled at 500 Hz. The time positions of the maxima will be less frequent and irregular, therefore interpolate them to constant 0.5 Hz.
			>>> _, abp_maxes = spt.findExtrema(channel_data)
			>>> # 500 Hz * 0.5 = 250
			>>> resampling_x = [ix for ix in range(0, len(channel_data), 250)]
			>>> abp_interp_maxes = spt.interpolate(abp_maxes, interp_x)
			>>> abp_interp_maxes
			... [(0, 0.0), (250, 121.0), (500, 125.0), (750, 117.0)]
		"""
		interpolated_y = []
		tyi = iter(timed_y)
		prev_tx, prev_ty = (0, 0.0)
		next_tx, next_ty = next(tyi)
		
		for sx in sample_x:
			if sx >= next_tx:
				try:
					prev_tx = next_tx
					prev_ty = next_ty
					next_tx, next_ty = next(tyi)
				except StopIteration:
					break;
			
			iy = ((next_ty - prev_ty) * (sx - prev_tx)) / (next_tx - prev_tx) + prev_ty
			
			interpolated_y.append( (sx, iy) )
		
		return interpolated_y
	

	def getSamplingFreq(self):
		if self.sampling_freq is None:
			raise ValueError("Sampling frequency has not been set yet! [ getSamplingFreq() ]")
	
		return self.sampling_freq
	

	def filterSamplingFreq(self, f):
		""" Get or set the global sampling frequency parameter. """
		
		if not isinstance(f, float):
			f = float(f)
		
		if f <= 0.0:
			raise ValueError("Sampling frequency must be larger than zero! [setSamplingFreq({})]".format(f))
		
		return f
	

	def calcSamplingRateFromFreq(self, chn_freq, new_freq):
		""" Calculate the resampling rate from frequency.
		
			The channel data is recorded at a given 'chn_freq' frequency, but the user may want to specify a 'new_freq' sampling frequency. This method returns the index stepping for the new sampling.
			
			For example, if
				chn_freq = 500 (data was recorded at 500 Hz)
				new_freq = 0.5 (user wants to sample the data at 0.5 Hz)
			then
				resampling rate = 1000
			because 0.5 Hz means 1 measurement per 2 seconds, and since the original data contains 500 records (measurements) per second, we have to take every 1000th record to achieve 0.5 Hz.
			
			See also
			--------
			SenkPolyTools.setSamplingFreq()
		"""
		
		if new_freq > chn_freq:
			raise ValueError("Requested sampling frequency can't be higher than the original channel record frequency! [calcSamplingRateFromFreq({}, {})]".format(chn_freq, new_freq))
		
		if chn_freq == 0:
			raise ValueError("Channel record frequency and sampling frequency can't be zero! [calcSamplingRateFromFreq({}, {})]".format(chn_freq, new_freq))
		
		return int(round(chn_freq/new_freq))


	def exportAnnotations(self, edfc):
		print("Annotations... ", end='')
		
		annot_fname = os.path.join(self.results_base_path, "{}_annotations.txt".format(edfc.file_basename))
		
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
	
	
	def processSenkEdf(self, edfc):
		""" Process an EDF file specific to the SE Neurology Clinic.
		
			This is a high level method which takes an EdfContainer object and processes and exports the ECG, TCD, ABP and CO2 channels.
			The output is a bunch of text files in the /results directory.
		"""
		
		for chn_i, chn_data in enumerate(edfc.data):
			chn_n = edfc.data_labels[chn_i]
			
			if 'EKG' in chn_n or 'ECG' in chn_n:
				spt.exportSenkEcg(edfc, chn_i)
			
			elif 'TCD' in chn_n:
				spt.exportSenkTcd(edfc, chn_i)
				
			elif 'Tonometry' in chn_n: # ABP
				spt.exportSenkAbp(edfc, chn_i)
			
			elif 'CO2' in chn_n:
				spt.exportSenkCo2(edfc, chn_i)
			
			else:
				print("Channel '{}'... skipped".format(chn_n))
	
	
	def exportSenkEcg(self, edfc, chn_num):
		chn_n = edfc.data_labels[chn_num]
		data = edfc.data[chn_num]
		ecg_freq = edfc.sample_freq[chn_num]
		
		print("Channel '{}'... ".format(chn_n), end="")
			
		sm_data = spt.smoothByAvg(data, 8)
		
		# export data after minimal noise reduction
		if self.export_nr_raw:
			fname = os.path.join(self.results_base_path, "{}_{}_nr_raw.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp, dialect='excel', delimiter=';')
				
				csvw.writerow(['Time (sec)', 'Data Value'])
				for sm_d_i, dv in enumerate(sm_data):
					csvw.writerow([sm_d_i/ecg_freq, dv])
		
		# find minima, maxima
		_, ecg_maxima = spt.findExtrema(sm_data)
		ecg_max_ind, __ = zip(*ecg_maxima)
		ecg_dists = spt.indexDists(ecg_max_ind)
		
		hr = [60/(dist/ecg_freq) for dist in ecg_dists]
		times = [curr_ecg_idx/ecg_freq for curr_ecg_idx in ecg_max_ind]
		#times.insert(0, 0.0)
		
		fname = os.path.join(self.results_base_path, "{}_{}_heart_rate.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
		with open(fname, "w", newline='') as fp:
			csvw = csv.writer(fp, dialect='excel', delimiter=';')
			
			csvw.writerow(["Time (sec)", "Heart Rate", "R-R Distance"])
			for ecg_t, ecg_hr, rr_dist in zip(times, hr, ecg_dists):
				csvw.writerow([ecg_t, ecg_hr, rr_dist])
		
		print("OK")
	
	
	def exportSenkTcd(self, edfc, chn_num):
		chn_n = edfc.data_labels[chn_num]
		data = edfc.data[chn_num]
		tcd_freq = edfc.sample_freq[chn_num]
		
		print("Channel '{}'... ".format(chn_n), end="")
		
		sm_data = spt.smoothByAvg(data)
		
		# export data after minimal noise reduction
		if self.export_nr_raw:
			fname = os.path.join(self.results_base_path, "{}_{}_nr_raw.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp, dialect='excel', delimiter=';')
				
				csvw.writerow(['Time (sec)', 'Data Value'])
				for sm_d_i, dv in enumerate(sm_data):
					csvw.writerow([sm_d_i/tcd_freq, dv])
		
		# find minima, maxima
		tcd_mins, tcd_maxes = spt.findExtrema(sm_data)
		
		# 0.5 Hz sampling (after interpolation)
		sm_data_len = len(sm_data)
		tcd_sampling = self.calcSamplingRateFromFreq(tcd_freq, self.getSamplingFreq())
		interp_x = [ix for ix in range(0, sm_data_len, tcd_sampling)]
		
		tcd_interp_mins = spt.interpolate(tcd_mins, interp_x)
		tcd_interp_maxes = spt.interpolate(tcd_maxes, interp_x)
		
		# 5 seconds moving average on interpolated minima
		_, tcd_imins_vals = zip(*tcd_interp_mins)
		tcd_imins_vals = np.array(tcd_imins_vals)
		tcd_imins_avg = spt.smoothByAvg(tcd_imins_vals, 10)
		
		# 5 seconds moving average on interpolated maxima
		_, tcd_imaxes_vals = zip(*tcd_interp_maxes)
		tcd_imaxes_vals = np.array(tcd_imaxes_vals)
		tcd_immaxes_avg = spt.smoothByAvg(tcd_imaxes_vals, 10)
		
		fname = os.path.join(self.results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
		with open(fname, "w", newline='') as fp:
			csvw = csv.writer(fp, dialect='excel', delimiter=';')
			csvw.writerow([
				"Time (sec)",
				"TCD Maxima",
				"TCD Minima",
				"(Max + 2*Min)/3",
				"TCD Maxima 5sec Moving Avg",
				"TCD Minima 5sec Moving Avg",
				"Moving Avg (Max + 2*Min)/3",
				"Resistance Index of Moving Avg"
			])
			
			for tcd_mn, tcd_mx, tcd_imn_a, tcd_imx_a in zip(tcd_interp_mins, tcd_interp_maxes, tcd_imins_avg, tcd_immaxes_avg):
				csvw.writerow([
					tcd_mn[0]/tcd_freq,
					tcd_mx[1],
					tcd_mn[1],
					(tcd_mx[1] + 2*tcd_mn[1])/3,
					tcd_imx_a,
					tcd_imn_a,
					(tcd_imx_a + 2*tcd_imn_a)/3,
					(tcd_imx_a - tcd_imn_a)/tcd_imx_a
				])
		
		print("OK")
	
	
	def exportSenkAbp(self, edfc, chn_num):
		chn_n = edfc.data_labels[chn_num]
		data = edfc.data[chn_num]
		abp_freq = edfc.sample_freq[chn_num]
		
		print("Channel '{}'... ".format(chn_n), end="")
		
		sm_data = spt.smoothByAvg(data, 8)
		
		# export data after minimal noise reduction
		if self.export_nr_raw:
			fname = os.path.join(self.results_base_path, "{}_{}_nr_raw.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp, dialect='excel', delimiter=';')
				
				csvw.writerow(['Time (sec)', 'Data Value'])
				for sm_d_i, dv in enumerate(sm_data):
					csvw.writerow([sm_d_i/abp_freq, dv])
		
		# find minima, maxima
		abp_mins, abp_maxes = spt.findExtrema(sm_data)
		
		# sampling settings
		sm_data_len = len(sm_data)
		abp_sampling = self.calcSamplingRateFromFreq(abp_freq, self.getSamplingFreq())
		interp_x = [ix for ix in range(0, sm_data_len, abp_sampling)]
		
		# heart rate from ABP for every half second, because ECG channel is often missing...
		abp_max_ind, __ = zip(*abp_maxes)
		abp_max_dists = spt.indexDists(abp_max_ind)
		hr = [60/(dist/abp_freq) for dist in abp_max_dists]
		
		abp_hr_interp = spt.interpolate(zip(abp_max_ind, hr), interp_x)
		
		# Half Hertz sampling of interpolated minima and maxima
		abp_interp_mins = spt.interpolate(abp_mins, interp_x)
		abp_interp_maxes = spt.interpolate(abp_maxes, interp_x)
		
		# 5 seconds moving average on interpolated minima
		_, abp_imins_vals = zip(*abp_interp_mins)
		abp_imins_vals = np.array(abp_imins_vals)
		abp_imins_avg = spt.smoothByAvg(abp_imins_vals, 10)
		
		# 5 seconds moving average on interpolated maxima
		_, abp_imaxes_vals = zip(*abp_interp_maxes)
		abp_imaxes_vals = np.array(abp_imaxes_vals)
		abp_imaxes_avg = spt.smoothByAvg(abp_imaxes_vals, 10)
		
		fname = os.path.join(self.results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
		
		with open(fname, "w", newline='') as fp:
			csvw = csv.writer(fp, dialect='excel', delimiter=';')
			csvw.writerow([
				"Time (sec)",
				"ABP Systole",
				"ABP Diastole",
				"(Syst + 2*Diast)/3",
				"ABP Systole 5sec Moving Avg",
				"ABP Diastole 5sec Moving Avg",
				"Moving Avg (Syst + 2*Diast)/3",
				"Heart Rate (Derived)"
			])
			
			for abp_mn, abp_mx, abp_imn_a, abp_imx_a, abp_hr in zip(abp_interp_mins, abp_interp_maxes, abp_imins_avg, abp_imaxes_avg, abp_hr_interp):
				csvw.writerow([
					abp_mn[0]/abp_freq,
					abp_mx[1],
					abp_mn[1],
					(abp_mx[1] + 2*abp_mn[1])/3,
					abp_imx_a,
					abp_imn_a,
					(abp_imx_a + 2*abp_imn_a)/3,
					abp_hr[1]
				])
		
		print("OK")
	
	
	def exportSenkCo2(self, edfc, chn_num):
		chn_n = edfc.data_labels[chn_num]
		data = edfc.data[chn_num]
		co2_freq = edfc.sample_freq[chn_num]
		
		print("Channel '{}'... ".format(chn_n), end="")
		
		sm_data = spt.smoothByAvg(data, 8)
		
		# export data after minimal noise reduction
		if self.export_nr_raw:
			fname = os.path.join(self.results_base_path, "{}_{}_nr_raw.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
			with open(fname, "w", newline='') as fp:
				csvw = csv.writer(fp, dialect='excel', delimiter=';')
				
				csvw.writerow(['Time (sec)', 'Data Value'])
				for sm_d_i, dv in enumerate(sm_data):
					csvw.writerow([sm_d_i/co2_freq, dv])
		
		# 0.5 Hz sampling parameters
		sm_data_len = len(sm_data)
		co2_sampling = self.calcSamplingRateFromFreq(co2_freq, self.getSamplingFreq())
		interp_x = [ix for ix in range(0, sm_data_len, co2_sampling)]
		
		# respiratory peaks
		_, co2_maxima = spt.findExtrema(sm_data, hwl=2000)
		
		# respiratory rate
		co2_max_ind, __ = zip(*co2_maxima)
		co2_dists = spt.indexDists(co2_max_ind)
		rr = [60/(dist/co2_freq) for dist in co2_dists]
		
		co2_rr_interp = spt.interpolate(zip(co2_max_ind, rr), interp_x)
		
		# .5 Hz sampling of respiratory peaks
		co2_interp_maxes = spt.interpolate(co2_maxima, interp_x)
		
		# 5 seconds moving average on interpolated maxima
		_, co2_imaxes_vals = zip(*co2_interp_maxes)
		co2_imaxes_vals = np.array(co2_imaxes_vals)
		co2_immaxes_avg = spt.smoothByAvg(co2_imaxes_vals, 10)
		
		fname = os.path.join(self.results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
		with open(fname, "w", newline='') as fp:
			csvw = csv.writer(fp, dialect='excel', delimiter=';')
			csvw.writerow([
				"Time (sec)",
				"CO2 Maxima",
				"CO2 Maxima 5sec Moving Avg",
				"Respiratory Rate"
			])
			
			for co2_mx, co2_imx_a, curr_rr in zip(co2_interp_maxes, co2_immaxes_avg, co2_rr_interp):
				csvw.writerow([
					(co2_mx[0]/co2_freq) - 7, # respiratory channel lag: 7 seconds
					co2_mx[1],
					co2_imx_a,
					curr_rr[1]
				])
		
		print("OK")
	



if __name__ == '__main__':
	argp = argparse.ArgumentParser(description="Build the global ComPPI network and export various subnetworks of it.")
	
	argp.add_argument(
		'-edf',
		required = True,
		help="The name of the input EDF file (only the file name, and it must be in the ./data/ folder)")
	
	argp.add_argument(
		'-s', '--sampling_frequency',
		help="Optional parameter to specify the sampling frequency (in Hz) at which the processed/interpolated data will be exported.\nThe value should be a number larger than zero (float or integer).\nDefault sampling frequency = 0.5")
	
	argp.add_argument(
		'-r', '--raw',
		action='store_true',
		help="Optional flag to export the raw data with minimal noise reduction (processed data will also be exported).")
	
	args = argp.parse_args()
	
	spt = SenkPolyTools(args.sampling_frequency, args.raw)
	edfc = spt.loadEdf(args.edf)
	
	spt.exportAnnotations(edfc)
	
	spt.processSenkEdf(edfc)