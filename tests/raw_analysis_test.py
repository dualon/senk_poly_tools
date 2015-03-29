""" Polygraph Tools, Semmelweis Egyetem Neurologiai Klinika

	The code is an experimental parser and analyzer for the EDF files that
	are results of various polygraph recordings.
	
	author: David M. 'dualon' Gyurko gyurko dot david at e-arc dot hu
	
	The EDF IO was inspired by Boris Reuderink's EEGTools:
	https://github.com/breuderink/eegtools/tree/master/eegtools

"""

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from senk_poly_tools.senk_poly_tools import SenkPolyTools


argp = argparse.ArgumentParser(description="Build the global ComPPI network and export various subnetworks of it.")
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

print(edfc.data_labels[22])
D = edfc.data[22]
print(len(D))

#windows = [ # Cs. G. counting
#	(717, 807, '#fff44f'),
#	(1992.5, 2082.5, '#FF8A00'),
#	(226.5, 316.5, '#AA0000'),
#	(423, 513, '#990000'),
#	(521, 611, '#11BEFF'),
#	(1109.5, 1199.5, '#00AA00'),
#	(1305.5, 1395.5, '#813D00'),
#	(1502, 1592, '#707070'),
#	(1698, 1788, '#1F0064'),
#	(1894, 1984, '#000000')
#]

# intervals: bogi, calc
#windows = [
#	#(-15, 75, ''),
#	#(384.5, 474.5, '#fff44f'), # light yellow
#	#(580.5, 670.5, '#FF8A00'), # orange
#	#(776.5, 866.5, '#AA0000'), # red
#	#(973, 1063, '#990000'), # maroon
#	#(1071, 1161, '#11BEFF'), # teal
#	(1561.5, 1651.5, '#00AA00'), # green
#	(1757.5, 1847.5, '#813D00'), # brown
#	(2052.0, 2142.0, '#707070'), # grey
#	(2150.0, 2240.0, '#1F0064'), # dark blue
#	(2248, 2338, '#000000')
#]

# intervals: bogi, verbal
windows = [
	#(-15, 75, ''),
	#(188, 278, '#fff44f'), # light yellow
	#(286, 376, '#FF8A00'), # orange
	#(482.5, 572.5, '#AA0000'), # red
	#(678.5, 768.5, '#990000'), # maroon
	#(874.5, 964.5, '#11BEFF'), # teal
	#(1169, 1259, '#00AA00'), # green
	(1267, 1357, '#813D00'), # brown
	(1365, 1455, '#707070'), # grey
	(1463.5, 1553.5, '#1F0064'), # dark blue
	(1659.5, 1749.5, '#000000')
]

import matplotlib.pyplot as plt
from operator import itemgetter

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


#w1_i1 = *500
#w1_i2 = *500
#w1 = D[w1_i1:w1_i2]
#print("w1 l: {}".format(len(w1)))
#
#w2_i1 = int(*500)
#w2_i2 = int(*500)
#w2 = D[w2_i1:w2_i2]
#print("w2 l: {}".format(len(w2)))
#
#_, co2_maxima_w1 = spt.findExtrema(w1, hwl=2000) # window: 150
#_, co2_maxima_w2 = spt.findExtrema(w2, hwl=2000) # window: 150
#
#from operator import itemgetter
#w1_mx = max(co2_maxima_w1, key=itemgetter(1))
#w1_mn = min(co2_maxima_w1, key=itemgetter(1))
#print("w1 min: {}, max: {}".format(w1_mn, w1_mx))
#
#print()
#
#print("w1 maxima")
#print(co2_maxima_w1)
#
#w1_x = []
#w1_y = []
#for mx_x, mx_y in co2_maxima_w1:
#	w1_x.append(mx_x)
#	w1_y.append(mx_y)
#	
#w2_x = []
#w2_y = []
#for mx_x, mx_y in co2_maxima_w2:
#	w2_x.append(mx_x)
#	w2_y.append(mx_y)

#import matplotlib.pyplot as plt
#		
#f,ax = plt.subplots(1,1,figsize=(30,15))
##ax.plot(w1, color='#AA0000')
#ax.plot(w1_x, w1_y, color='#660000')
#
##ax.plot(w2, color='#00AA00')
#ax.plot(w2_x, w2_y, color='#006600')
#
#ax.vlines(w1_x, -50.0, 50.0, color='#990000', alpha=0.5)
#ax.vlines(w2_x, -50.0, 50.0, color='#009900', alpha=0.5)
#plt.savefig(os.path.join('results', '{}_w1_w2.png'.format(edfc.file_basename)))
#plt.close()

#for chn_i, chn_data in enumerate(edfc.data):
#	chn_n = edfc.data_labels[chn_i]
#	chn_data = edfc.data[chn_i]
#	print("Channel '{}'... ".format(chn_n), end="")
#	
#	if 'EKG' in chn_n or 'ECG' in chn_n:
#		ecg_freq = edfc.sample_freq[chn_i]
#		
#		sm_data = spt.smoothByAvg(chn_data, 8)
#		_, ecg_maxima = spt.findExtrema(sm_data)
#		ecg_max_ind, __ = zip(*ecg_maxima)
#		ecg_dists = spt.indexDists(ecg_max_ind)
#		hr = [60/(dist/ecg_freq) for dist in ecg_dists]
#		times = [curr_ecg_idx/ecg_freq for curr_ecg_idx in ecg_max_ind]
#		times.insert(0, 0.0)
#		
#		fname = os.path.join(results_base_path, "{}_{}_heart_rate.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
#		with open(fname, "w", newline='') as fp:
#			csvw = csv.writer(fp, dialect='excel', delimiter=';')
#			
#			csvw.writerow(['Time (sec)', 'Heart Rate'])
#			for ecg_t, ecg_hr in zip(times, hr):
#				csvw.writerow([ecg_t, ecg_hr])
#		
#		print("OK")
#	
#	elif 'TCD' in chn_n:
#		sm_data = spt.smoothByAvg(chn_data)
#		tcd_mins, tcd_maxes = spt.findExtrema(sm_data)
#		
#		sm_data_len = len(sm_data)
#		tcd_sampling = int(edfc.sample_freq[chn_i]*0.5) # 500 Hz * 0.5
#		interp_x = [ix for ix in range(0, sm_data_len, tcd_sampling)]
#		
#		tcd_interp_mins = spt.interpolate(tcd_mins, interp_x)
#		tcd_interp_maxes = spt.interpolate(tcd_maxes, interp_x)
#		
#		# 5 seconds moving average on interpolated minima
#		_, tcd_imins_vals = zip(*tcd_interp_mins)
#		tcd_imins_vals = np.array(tcd_imins_vals)
#		tcd_imins_avg = spt.smoothByAvg(tcd_imins_vals, 10)
#		
#		# 5 seconds moving average on interpolated maxima
#		_, tcd_imaxes_vals = zip(*tcd_interp_maxes)
#		tcd_imaxes_vals = np.array(tcd_imaxes_vals)
#		tcd_immaxes_avg = spt.smoothByAvg(tcd_imaxes_vals, 10)
#		
#		fname = os.path.join(results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
#		with open(fname, "w", newline='') as fp:
#			csvw = csv.writer(fp, dialect='excel', delimiter=';')
#			csvw.writerow(["Time (sec)", "TCD Maxima", "TCD Minima", "(Max + 2*Min)/3", "TCD Maxima 5sec Moving Avg", "TCD Minima 5sec Moving Avg", "Moving Avg (Max + 2*Min)/3", "Resistance Index of Moving Avg"])
#			
#			for tcd_mn, tcd_mx, tcd_imn_a, tcd_imx_a in zip(tcd_interp_mins, tcd_interp_maxes, tcd_imins_avg, tcd_immaxes_avg):
#				csvw.writerow([
#					tcd_mn[0]/edfc.sample_freq[chn_i],
#					tcd_mx[1],
#					tcd_mn[1],
#					(tcd_mx[1] + 2*tcd_mn[1])/3,
#					tcd_imx_a,
#					tcd_imn_a,
#					(tcd_imx_a + 2*tcd_imn_a)/3,
#					(tcd_imx_a - tcd_imn_a)/tcd_imx_a
#				])
#		
#		print("OK")
#		
#	elif 'Tonometry' in chn_n: # ABP
#		sm_data = spt.smoothByAvg(chn_data, 8)
#		abp_mins, abp_maxes = spt.findExtrema(sm_data)
#		
#		# sampling settings
#		abp_freq = edfc.sample_freq[chn_i]
#		sm_data_len = len(sm_data)
#		abp_sampling = int(abp_freq*0.5) # 500 Hz * 0.5
#		interp_x = [ix for ix in range(0, sm_data_len, abp_sampling)]
#		
#		# heart rate from ABP for every half second, because ECG channel is often missing...
#		abp_max_ind, __ = zip(*abp_maxes)
#		abp_max_dists = spt.indexDists(abp_max_ind)
#		hr = [60/(dist/abp_freq) for dist in abp_max_dists]
#		
#		abp_hr_interp = spt.interpolate(zip(abp_max_ind, hr), interp_x)
#		
#		# Half Hertz sampling of interpolated minima and maxima
#		abp_interp_mins = spt.interpolate(abp_mins, interp_x)
#		abp_interp_maxes = spt.interpolate(abp_maxes, interp_x)
#		
#		# 5 seconds moving average on interpolated minima
#		_, abp_imins_vals = zip(*abp_interp_mins)
#		abp_imins_vals = np.array(abp_imins_vals)
#		abp_imins_avg = spt.smoothByAvg(abp_imins_vals, 10)
#		
#		# 5 seconds moving average on interpolated maxima
#		_, abp_imaxes_vals = zip(*abp_interp_maxes)
#		abp_imaxes_vals = np.array(abp_imaxes_vals)
#		abp_imaxes_avg = spt.smoothByAvg(abp_imaxes_vals, 10)
#		
#		fname = os.path.join(results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
#		
#		with open(fname, "w", newline='') as fp:
#			csvw = csv.writer(fp, dialect='excel', delimiter=';')
#			csvw.writerow(["Time (sec)", "ABP Systole", "ABP Diastole", "(Syst + 2*Diast)/3", "ABP Systole 5sec Moving Avg", "ABP Diastole 5sec Moving Avg", "Moving Avg (Syst + 2*Diast)/3", "Heart Rate (Derived)"])
#			
#			for abp_mn, abp_mx, abp_imn_a, abp_imx_a, abp_hr in zip(abp_interp_mins, abp_interp_maxes, abp_imins_avg, abp_imaxes_avg, abp_hr_interp):
#				csvw.writerow([
#					abp_mn[0]/abp_freq,
#					abp_mx[1],
#					abp_mn[1],
#					(abp_mx[1] + 2*abp_mn[1])/3,
#					abp_imx_a,
#					abp_imn_a,
#					(abp_imx_a + 2*abp_imn_a)/3,
#					abp_hr[1]
#				])
#		
#		print("OK")
#	
#	elif 'CO2' in chn_n:
#		sm_data = spt.smoothByAvg(chn_data, 8)
#		
#		sm_data_len = len(sm_data)
#		co2_freq = edfc.sample_freq[chn_i]
#		co2_sampling = int(co2_freq*0.5) # 500 Hz * 0.5
#		interp_x = [ix for ix in range(0, sm_data_len, co2_sampling)]
#		
#		# respiratory peaks
#		_, co2_maxima = spt.findExtrema(sm_data) # window: 150
#		
#		# respiratory rate
#		co2_max_ind, __ = zip(*co2_maxima)
#		co2_dists = spt.indexDists(co2_max_ind)
#		rr = [60/(dist/co2_freq) for dist in co2_dists]
#		
#		co2_rr_interp = spt.interpolate(zip(co2_max_ind, rr), interp_x)
#		
#		#import matplotlib.pyplot as plt
#		#
#		#f,ax = plt.subplots(1,1,figsize=(80,8))
#		#ax.plot(sm_data[0:50000])
#		#ax.vlines([idx for idx, _ in co2_maxima[0:21]], -150.0, 150.0, color='#ED9A00', alpha=0.5)
#		#plt.savefig(os.path.join('results', '{}.png'.format(edfc.file_basename)))
#		#plt.close()
#		
#		# .5 Hz sampling of respiratory peaks
#		co2_interp_maxes = spt.interpolate(co2_maxima, interp_x)
#		
#		# 5 seconds moving average on interpolated maxima
#		_, co2_imaxes_vals = zip(*co2_interp_maxes)
#		co2_imaxes_vals = np.array(co2_imaxes_vals)
#		co2_immaxes_avg = spt.smoothByAvg(co2_imaxes_vals, 10)
#		
#		fname = os.path.join(results_base_path, "{}_{}_0.5hz.txt".format(edfc.file_basename, chn_n.replace(' ', '')))
#		with open(fname, "w", newline='') as fp:
#			csvw = csv.writer(fp, dialect='excel', delimiter=';')
#			csvw.writerow(["Time (sec)", "CO2 Maxima", "CO2 Maxima 5sec Moving Avg", "Respiratory Rate"])
#			
#			for co2_mx, co2_imx_a, curr_rr in zip(co2_interp_maxes, co2_immaxes_avg, co2_rr_interp):
#				csvw.writerow([
#					(co2_mx[0]/co2_freq) - 7, # respiratory channel lag: 7 seconds
#					co2_mx[1],
#					co2_imx_a,
#					curr_rr[1]
#				])
#		
#		print("OK")
#	
#	else:
#		print("skipped")