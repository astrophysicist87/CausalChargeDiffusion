#!/usr/bin/env python

from numpy import *
from numpy.core import numeric as _nx
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.patches import Ellipse
import sys, os

mpl.rcParams['pdf.fonttype'] = 42
labelsize = 20
mpl.rcParams['xtick.labelsize'] = labelsize
mpl.rcParams['ytick.labelsize'] = labelsize

#################################################################
# Parameters to run script with
#################################################################
#file information
#filename = sys.argv[1]

#grid information
nDY = 51

panelLabels = ['(a)', '(b)']
panelCounter = 0

CBFIndices = ['\pi\pi','p \\bar{p}','K K']
#filenames = ['ecf_pion_T0_250MeV.out', 'ecf_proton_T0_250MeV.out', 'ecf_kaon_T0_250MeV.out']
filenames = ['ecf_protonKaon_T0_350MeV.out']
#filenames = ['ecf_pion_T0_350MeV.out', 'ecf_proton_T0_350MeV.out', 'ecf_kaon_T0_350MeV.out']
ExpDataFilenames = ['star_pK.dat']
#ExpDataFilenames = ['star_pipi.dat', 'star_ppbar.dat', 'star_KK.dat']

#snapshotFractions = ['0.0','0.2','0.4','0.6','0.8','1.0']
snapshotFractions = ['0.01','0.025','0.1','0.25','0.5','1.0']
lineColors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
plotTitles = ['Lattice', r'$2\pi D_Q T = 0.5$', r'$2\pi D_Q T = 1.0$', r'$2\pi D_Q T = 1.5$']

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

'''def stack(arrays, axis=0):
	arrays = [np.asanyarray(arr) for arr in arrays]
	if not arrays:
		raise ValueError('need at least one array to stack')
	
	shapes = set(arr.shape for arr in arrays)
	if len(shapes) != 1:
		raise ValueError('all input arrays must have the same shape')
	
	result_ndim = arrays[0].ndim + 1
	if not -result_ndim <= axis < result_ndim:
		msg = 'axis {0} out of bounds [-{1}, {1})'.format(axis, result_ndim)
		raise IndexError(msg)
	if axis < 0:
		axis += result_ndim
	
	sl = (slice(None),) * axis + (_nx.newaxis,)
	expanded_arrays = [arr[sl] for arr in arrays]
	return _nx.concatenate(expanded_arrays, axis=axis)'''

#################################################################
# Plot correlations
def plotCorrelations(index):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0
	filename = filenames[index]

	nCols = 2
	chosenCols = [0, 1, 3, 5, 7]

	# read in file
	data = loadtxt(filename, usecols=tuple(chosenCols))
	expdata = loadtxt(ExpDataFilenames[index])

	#print 'Loading data file...'
	ax.plot(data[:,0], data[:,1], linestyle='-', color='red', linewidth=lw, label='Lattice')
	ax.plot(data[:,0], data[:,2], linestyle='-', color='blue', linewidth=lw, label=r'$2\pi D_Q T = 0.5$')
	ax.plot(data[:,0], data[:,3], linestyle='-', color='green', linewidth=lw, label=r'$2\pi D_Q T = 1.0$')
	ax.plot(data[:,0], data[:,4], linestyle='-', color='purple', linewidth=lw, label=r'$2\pi D_Q T = 1.5$')
	#ax.plot(expdata[:,0], expdata[:,1], linestyle='None', color='black', marker='s', markersize=5, label='STAR')
	ax.errorbar(expdata[:,0], expdata[:,1], yerr=expdata[:,2], linestyle='None', color='black', marker='s', markersize=5, label='STAR')
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$B_{%(idx)s}$' % {'idx': CBFIndices[index]}, fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	#text(0.9, 0.15, r'(b)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	#plt.title(pathname)
	
	#plt.show(block=False)
	#plt.show()
	outfilename = os.path.splitext(filename)[0] + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename


#################################################################
# Plot snapshots
def plotSnapshots(columnToPlot):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0

	index = 0
	nSnapshots = len(snapshotFractions)
	chosenCols = [0, 1, 3, 5, 7]

	for isnap in xrange(nSnapshots):
		# read in file
		filename = 'ecf_pion_snapshot_%(sf)s.out' % {'sf': snapshotFractions[isnap]}
		data = loadtxt(filename, usecols=tuple(chosenCols))

		#print 'Loading data file...'
		ax.plot(data[:,0], data[:,columnToPlot+1], linestyle='-', color=lineColors[isnap], linewidth=lw, label=r'$\tau_f = %(sf)s \times \tau_{FO}$' % {'sf': snapshotFractions[isnap]})

	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$B_{%(idx)s}$' % {'idx': CBFIndices[index]}, fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=2, prop={'size': plotfontsize+5})
	ax.set_title(plotTitles[columnToPlot],fontsize=plotfontsize+5)
	
	#plt.show(block=False)
	#plt.show()
	outfilename = 'ecf_pion_snapshots_ParamV' + str(columnToPlot) + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename



def generate_all_plots():
	for i in xrange(len(filenames)):
		plotCorrelations(i)
	#for i in xrange(len(plotTitles)):
	#	plotSnapshots(i)
	#pause()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
