#!/usr/bin/env python

from numpy import *
#from numpy.core import numeric as _nx
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
#from scipy import interpolate
#from scipy.interpolate import griddata
#from matplotlib.patches import Ellipse
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

#CBFIndices = ['\pi\pi','p \\bar{p}','K K']
#filenames = ['ecf_pion_T0_250MeV.out', 'ecf_proton_T0_250MeV.out', 'ecf_kaon_T0_250MeV.out']
#filenames = ['ecf_protonKaon_T0_350MeV.out']
#filenames = ['ecf_pion_T0_350MeV.out', 'ecf_proton_T0_350MeV.out', 'ecf_kaon_T0_350MeV.out']
#ExpDataFilenames = ['star_pK.dat']
#ExpDataFilenames = ['star_pipi.dat', 'star_ppbar.dat', 'star_KK.dat']
vQ2_labels = ['0_1', '0_333', '0_5', '1_0', '2_0', '5_0', '10_0']
noiseType = 'white'
GreenType = 'color'
filenames = ['results_%(NT)s_noise_%(GT)s_Green_vQ2_%(vQ2s)s.out' % {'NT': noiseType, 'GT': GreenType, 'vQ2s': vQ2} for vQ2 in vQ2_labels]


lineColors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
#plotTitles = ['Lattice', r'$2\pi D_Q T = 0.5$', r'$2\pi D_Q T = 1.0$', r'$2\pi D_Q T = 1.5$']

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

#################################################################
# Plot correlations
def plotCorrelations():
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0
	subtractSelfCorrelations=True

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(filenames))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	
	idx=0
	for filename in filenames:
		columnToPlot = 1
		if subtractSelfCorrelations:
			columnToPlot = 2
		chosenCols = [0, 1, 3]
		
		colorVal = scalarMap.to_rgba(len(filenames)-idx-1)
		#colorVal = scalarMap.to_rgba(idx)

		# read in file
		data = loadtxt(filename, usecols=tuple(chosenCols))

		#print 'Loading data file...'
		ax.plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw)
		idx+=1
	
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$B_{\pi\pi}$', fontsize = labelsize + 10)
	#ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	#text(0.9, 0.15, r'(b)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	#plt.title(pathname)
	
	#plt.show(block=False)
	plt.show()
	#outfilename = os.path.splitext(filename)[0] + '.pdf'
	#plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	#print 'Saved to', outfilename


def generate_all_plots():
	plotCorrelations()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
