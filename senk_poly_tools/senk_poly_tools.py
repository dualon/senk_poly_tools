""" Polygraph Tools, Semmelweis Egyetem Neurologiai Klinika

	The code is an experimental parser and analyzer for the EDF files that
	are results of various polygraph recordings.
	
	author: David M. 'dualon' Gyurko gyurko dot david at e-arc dot hu
	
	The EDF IO was inspired by Boris Reuderink's EEGTools:
	https://github.com/breuderink/eegtools/tree/master/eegtools

"""

import os
import argparse
import pprint
import datetime
import numpy as np
import matplotlib.pyplot as plt
from edf_container import EdfContainer

class SenkPolyTools(object):
	
	annot_label = 'EDF Annotations'
	
	def openEdfFile(self, file_path):
		""" Open an EDF(+) file and return an empty EdfContainer object. """
		
		if not os.path.isfile(file_path):
			raise FileNotFoundError("File '{}' not found!".format(file_path))
		
		e = EdfContainer(file_path)
		e.file_obj = open(file_path, 'rb')
		
		return e
	
	
	def loadEdfHeaders(self, e):
		""" Load the headers of an EDF(+) file into an EdfContainer.
		
			Existing headers are overwritten.
			
			param e: (empty) EdfContainer object,
				see SenkPolyTools.openEdfFile()
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
		edf_type = f.read(44)[:5]
		e.edf_type = edf_type.decode(encoding='ASCII')
		e.num_of_records = int(f.read(8))
		e.record_duration = float(f.read(8))
		e.num_of_signals = int(f.read(4))
		
		nsr = range(e.num_of_signals)
		#print(nsr)
		e.labels = [f.read(16).strip().decode(encoding='ASCII') for _ in nsr]
		e.transducer_types 		= [f.read(80).strip().decode(encoding='ASCII') for _ in nsr]
		
		e.physical_dimension 	= [f.read(8).strip().decode(encoding='ASCII') for _ in nsr] # physical_dimensions: uV, cm/s, %, mV, ...
		e.physical_min 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.physical_max 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.digital_min 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.digital_max 			= np.asarray([float(f.read(8)) for _ in nsr])
		e.gain 					= (e.physical_max - e.physical_min) / (e.digital_max - e.digital_min) # numpy arrays
		e.prefiltering 			= [f.read(80).strip() for _ in nsr]
		e.num_of_samples_per_record = [int(f.read(8)) for _ in nsr]
		
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
			raise AttributeError("EdfContainer.file_obj is missing, use SenkPolyTools.openEdfFile() to create a file stream.")
		
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
					phys = phys.tolist()
					data[k].extend(phys) # note the 'k' index
					
					k = k+1
		
		e.data_labels = [l for l in e.labels if l != self.annot_label]
		e.data = np.array(data) # e.data is a numpy 2D matrix
		del data
		
		return e

	
	def printEdfHeader(self, e):
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



if __name__ == '__main__':
	argp = argparse.ArgumentParser(description="Build the global ComPPI network and export various subnetworks of it.")
	argp.add_argument(
		'-edf',
		help="The name of the input EDF file (must be in the ./data folder)")
	args = argp.parse_args()
	
	spt = SenkPolyTools()
	
	edfc = spt.openEdfFile(os.path.join('..', 'data', args.edf))
	edfc = spt.loadEdfHeaders(edfc)
	spt.printEdfHeader(edfc)
	edfc = spt.loadEdfData(edfc)
	
	spt.visualizeEdf(edfc)