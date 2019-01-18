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

CBFIndices = ['\pi^+\pi^-','p \\bar{p}','K^+ K^-']
particleLabels = ['pi','p','K']
#vQ2Labels = ['0_100', '0_333', '1_000', '3_000', '100_000']
#vQ2LineLabels = ['1/10', '1/3', '1', '3', 'White']
#vQ2Labels = ['100_000', '3_000', '1_000', '0_333', '0_100']
#vQ2LineLabels = ['White', '3', '1', '1/3', '1/10']
vQ2Labels = ['100_000', '1_000', '0_333', '0_100']
vQ2LineLabels = ['White', '1', '1/3', '1/10']
colorVals = ['black', 'blue', 'green', 'red']
#vQ2Labels = ['0_333']
#vQ2LineLabels = ['1/3']
vQ2LineStyles = {'0_100' : ':', '0_333' : '-.', '1_000' : '--', '100_000' : '-'}
vQ2DashStyles = {'0_100' : (100,100), '0_333' : (100,100), '1_000' : (100,100), '100_000' : (100,100)}

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
def generate_plot_panelA(fraction):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	dndn_filenames = 'dndn_Dxi_f_%(frac)s_vQ2_%(vQ2)s.out'

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(vQ2Labels)+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	for vQ2 in ['0_333']:
		filename = dndn_filenames % {'vQ2': vQ2, 'frac': fraction}
				
		#colorVal = scalarMap.to_rgba(len(vQ2Labels)+1-idx-1)
		#colorVal = scalarMap.to_rgba(idx)
		colorVal = 'black'

		if os.path.isfile(filename) and file_len(filename) > 0:
			# read in file
			data = loadtxt(filename)
			lineLabel = r'$v_Q^2 = 1/3$'

			#print 'Loading data file...'
			ax.plot(data[:,0], data[:,1], color=colorVal, linestyle='-', linewidth=lw, label=lineLabel)
		idx+=1
	
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta \xi$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$\left< \delta n(\Delta \xi) \delta n(0) \right>$', fontsize = labelsize + 10)
	#ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5}, frameon=False)
	ax.text(0.8, 0.85, r'$v_Q^2=1/3$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	ax.set_xlim(-1.0, 1.0)
	ax.set_ylim(bottom=0.0, top=31.0)
	#plt.title(fileStem + ': vQ2 comparison')
	
	#plt.show(block=False)
	outfilename = 'dndn_Dxi_f_0_05_vQ2_0_333.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename





#################################################################
def generate_plot_panelB(fraction):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	dndn_filenames = 'dndn_Dxi_f_%(frac)s_vQ2_%(vQ2)s.out'

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(vQ2Labels)+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	for vQ2 in vQ2Labels:
		filename = dndn_filenames % {'vQ2': vQ2, 'frac': fraction}
				
		#colorVal = scalarMap.to_rgba(len(vQ2Labels)+1-idx-1)
		#colorVal = scalarMap.to_rgba(idx)
		colorVal = colorVals[idx]

		if os.path.isfile(filename) and file_len(filename) > 0:
			# read in file
			data = loadtxt(filename)
			if vQ2 == '100_000':
				lineLabel = vQ2LineLabels[idx]
			else:
				lineLabel = r'$v_Q^2 = $' + vQ2LineLabels[idx]

			#print 'Loading data file...'
			ax.plot(data[:,0], data[:,1], color=colorVal, linestyle=vQ2LineStyles[vQ2], linewidth=lw, label=lineLabel)
		idx+=1
	
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta \xi$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$\left< \delta n(\Delta \xi) \delta n(0) \right>$', fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5}, frameon=False)
	#ax.text(0.25, 0.85, r'$v_Q^2=1/3$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	ax.set_xlim(-1.0, 1.0)
	ax.set_ylim(bottom=0.0, top=70.0)
	#plt.title(fileStem + ': vQ2 comparison')
	
	#plt.show(block=False)
	outfilename = 'dndn_Dxi_f_0_05_vQ2_comparison_PRCAPPROVED.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename


#############################################

if __name__ == "__main__":
	#plt.ion()
	#generate_plot_panelA('0_05')
	generate_plot_panelB('0_05')
	#pause()
	print 'Finished.'

# End of file
