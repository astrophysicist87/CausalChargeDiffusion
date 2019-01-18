#!/usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import sys, os

mpl.rcParams['pdf.fonttype'] = 42
labelsize = 20
mpl.rcParams['xtick.labelsize'] = labelsize
mpl.rcParams['ytick.labelsize'] = labelsize

#################################################################
# Parameters to run script with
#################################################################
panelLabels = ['(a)', '(b)']
panelCounter = 0

#resultsDir='results/'

CBFIndices = ['\pi\pi','p \\bar{p}','K K']
particleLabels = ['pi','p','K']
#vQ2Labels = ['0_100', '0_333', '1_000', '3_000']
#vQ2LineLabels = ['1/10', '1/3', '1', '3']
vQ2Labels = ['3_000', '1_000', '0_333', '0_100']
vQ2LineLabels = ['3', '1', '1/3', '1/10']
lineStyles = {'0_100' : ':', '0_333' : '-.', '1_000' : '--' , '3_000' : '-'}
#vQ2Labels = ['0_100', '0_333', '1_000']
#vQ2LineLabels = ['1/10', '1/3', '1']
colorVals = ['red', 'blue', 'green', 'purple']
#vQ2Labels = ['0_333']
#vQ2LineLabels = ['1/3']
arrowheadLocs = {'0_333' : [0.6,0.55], '1_000' : [0.55,0.5], '0_100' : [0.5,0.4] , '3_000' : [0.4,0.4]}

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def file_len(fname):
	try:
		with open(fname) as f:
			for i, l in enumerate(f):
			    pass
		return i + 1
	except UnboundLocalError:
		return 0
		



#################################################################
def generate_plot():
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0

	dndn_filenames = 'SC_Dxi_vQ2_%(vQ2)s.out'

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(vQ2Labels)+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	arrowLength = 0.055
	for vQ2 in vQ2Labels:
		filename = dndn_filenames % {'vQ2': vQ2}
				
		#colorVal = scalarMap.to_rgba(len(vQ2Labels)+1-idx-1)
		#colorVal = scalarMap.to_rgba(idx)
		colorVal = colorVals[idx]

		if os.path.isfile(filename) and file_len(filename) > 0:
			# read in file
			data = loadtxt(filename)
			lineLabel = r'$v_Q^2 = \,$' + vQ2LineLabels[idx]

			#print 'Loading data file...'
			ax.plot(data[:,0], data[:,1], color=colorVal, linestyle='-', linewidth=lw, label=lineLabel)
			arrowheadLoc = [0.06, 0.9*max(data[:,1])]
			ax.annotate('', xy=(arrowheadLoc[0],arrowheadLoc[1]),\
                            xycoords='data', xytext=(arrowheadLoc[0]+arrowLength,arrowheadLoc[1]),\
                            textcoords='data', arrowprops=dict(facecolor='black', width=1, headwidth=5, frac=0.15))
			plt.text(arrowheadLoc[0]+0.07, arrowheadLoc[1],\
                     lineLabel, horizontalalignment='left', verticalalignment='center',\
                     transform = ax.transData, size=20)
		idx+=1
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta \xi$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$ \left< \delta n(\Delta \xi) \delta n(0) \right> _{\rm self}$', fontsize = labelsize + 10)
	#ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5}, frameon=False)
	#ax.text(0.25, 0.85, r'$v_Q^2=1/3$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	ax.set_xlim(-0.25, 0.25)
	ax.set_ylim(bottom=-0.1)
	#plt.title(fileStem + ': vQ2 comparison')
	
	plt.show(block=False)
	#outfilename = 'SelfCorrelationsVsVQ2_PRCAPPROVED.pdf'
	#plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	#print 'Saved to', outfilename


#############################################

if __name__ == "__main__":
	#plt.ion()
	generate_plot()
	pause()
	#print 'Finished.'

# End of file
