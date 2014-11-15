""" A structured data container for EDF(+) header and data records.
	
	An EdfContainer instance represents a single EDF file with the headers and
	all its data records.
	Note that the data container should *not* require all fields to be set and
	it should not enforce any validation.
	
	For the EDF(+) specification see:
	http://www.edfplus.info/specs/index.html
"""

import os

class EdfContainer(object):
	def __init__(self, file_path):
		""" Class attributes are defined in constructor to ensure instance-specificity. """
		if not os.path.isfile(file_path):
			raise InputError("File '{}' not found!".format(file_path))
		
		f_basename = os.path.basename(file_path)
		
		self.file_path 					= file_path
		self.file_basename				= os.path.splitext(f_basename)[0]
		self.file_obj 					= None
		self.version 					= None
		self.local_patient_id 			= None
		self.local_recording_id 		= None
		self.date_time 					= None
		self.num_of_bytes_in_header 	= None
		self.edf_type 					= None
		self.num_of_records 			= None
		self.record_duration 			= None
		self.num_of_signals 			= None
		self.labels 					= None
		self.data_labels				= None
		self.transducer_types 			= None
		self.physical_dimension 		= None
		self.physical_min 				= None
		self.physical_max 				= None
		self.digital_min 				= None
		self.digital_max 				= None
		self.prefiltering 				= None
		self.num_of_samples_per_record 	= None
		self.sample_freq				= None

		self.annotations 				= []
		self.data 						= []