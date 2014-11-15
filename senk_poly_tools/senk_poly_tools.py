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
import numpy as np
from math import floor
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from senk_poly_tools.edf_container import EdfContainer

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
					e.annotations.append(buff)
				else:
					# if there is no entry for the channel data, create one
					if len(data) == j:
						data.append([])
					
					# 2-byte little endian integer per channel
					dig = np.fromstring(buff, '<i2').astype(np.float32)
					phys = (dig - e.digital_min[j]) * e.gain[j] + e.physical_min[j]
					
					# list extension + numpy arrayconversion at the end
					# is faster than converting each record to numpy array
					# @TODO: consider numpy.r_
					phys = phys.tolist()
					data[k].extend(phys) # note the 'k' index
					
					k = k+1
		
		e.data_labels = [l for l in e.labels if l != self.annot_label]
		e.data = np.array(data) # e.data is a numpy 2D matrix
		del data
		
		return e

	
	def printEdfHeaders(self, e):
		edf_h = vars(e)
		pprint.pprint(edf_h)
	

	def visualizeEdf(self, edf, channels = [], fig_dpi = 200):
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
			figsize=(ax_num*4, 30),
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
					color='#000000'
				)
				
				axes[chn].set_title(edf.data_labels[chn], loc='left')
		
		plt.savefig(os.path.join('..', 'results', '{}.png'.format(edf.file_basename)))
		plt.close()
	
	
	def smoothByAvg(self, x, window_len = 3):
		""" Smooth a curve with moving average.
		
			A series of overlapping, window_len length sections are created from the data (these are the moving windows), and they are convoluted against ones.
			
			Smoothing by moving average is a basic smoothing function for non-periodic functions (as TCD data is typically not periodic, for example; see FFT for periodic data).
			
			param x numpy array containing the data (typically a channel)
			param window_len int, the length of the moving window
			returns a numpy array of the smoothed data
		"""
		s = np.r_[ x[window_len-1:0:-1], x, x[-1:-window_len:-1] ]
		w = np.ones(window_len,'d')
		c = np.convolve(w/w.sum(), s, mode='valid')
		
		# the smoothed array will have the length of (the original array + window_len -1),
		# but the first window_len-1 points are not averaged
		# -> cut window_len/2 points at the start of the results
		return c[floor(window_len/2):]
	
		
	
	#def findMin(self, x):
	#	"""
	#		
	#		@TODO: consider http://docs.scipy.org/doc/scipy/reference/signal.html
	#	"""
	#	arr_len = len(arr)
	#	window = 50 # 
	#	
	#	for i in range(arr_len):
	#		
		
	



if __name__ == '__main__':
	argp = argparse.ArgumentParser(description="Build the global ComPPI network and export various subnetworks of it.")
	argp.add_argument(
		'-edf',
		help="The name of the input EDF file (must be in the ./data folder)")
	args = argp.parse_args()
	
	spt = SenkPolyTools()
	
	edfc = spt.createEdfContainer(os.path.join('..', 'data', args.edf))
	edfc = spt.loadEdfHeaders(edfc)
	#spt.printEdfHeaders(edfc)
	edfc = spt.loadEdfData(edfc)
	
	#sm_ch = spt.smooth(edfc.data[23])
	print(edfc.data_labels)
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 8), dpi=300)
	fig.tight_layout()
	
	y1 = edfc.data[19][0:1000]
	x1 = np.arange(1000)
	ax1.plot(x1, y1, color='#404040', alpha=0.5)
	
	y2 = spt.smoothByAvg(edfc.data[19][0:1000])
	x2 = np.arange(y2.shape[0])
	ax1.plot(x2, y2, color='#AA0000', alpha=0.5)
	
	ax1.set_title(edfc.data_labels[19] + ' Original (grey) vs Smooth (red)', loc='left')
	
	y3 = edfc.data[23][0:1000]
	ax2.plot(x1, y3, color='#404040', alpha=0.5)
	
	y4 = spt.smoothByAvg(edfc.data[23][0:1000])
	x4 = np.arange(y4.shape[0])
	ax2.plot(x4, y4, color='#AA0000', alpha=0.5)
	
	ax2.set_title(edfc.data_labels[23] + ' Original (grey) vs Smooth (red)', loc='left')
	
	plt.savefig(os.path.join('..', 'results', '{}_orig_vs_smoothbyavg.png'.format(edfc.file_basename)))
	plt.close()
	
	#spt.visualizeEdf(edfc)