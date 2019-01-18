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
mpl.rcParams['xtick.major.pad']='12'
mpl.rcParams['ytick.major.pad']='8'

#################################################################
# Parameters to run script with
#################################################################
panelLabels = ['(a)', '(b)']
panelCounter = 0

CBFIndices = ['\pi^+ \pi^-','p \\bar{p}','K^+ K^-']
particleLabels = ['pi','p','K']
vQ2Labels = ['0_333', '1_000', '10_000', '100_000']
vQ2LineLabels = {'0_333' : r'$v_Q^2=1/3$', '1_000' : r'$v_Q^2=1$', '10_000' : r'$v_Q^2=10$' , '100_000' : 'White'}
vQ2LineStyles = {'0_333' : ':', '1_000' : '-.', '10_000' : '--', '100_000' : '-'}
snapshotLabels = ['0_05', '0_10', '0_15', '0_25', '0_50', '1_00']
colorVals = {'0_333' : 'blue', '1_000' : 'green', '10_000' : 'red' , '100_000' : 'black'}
lineWidthVals = {'0_333' : 2.0, '1_000' : 2.0, '10_000' : 2.0 , '100_000' : 2.0}

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
def generate_plot_vQ2_comparison(particleIndex):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	subtractSelfCorrelations = True

	dndn_color_filenames = 'results/twoPC_%(PL)s%(PL)s_color_noise_color_Green_vQ2_%(vQ2)s_snapshot_1_00.dat'
	dndn_white_filenames = 'results/twoPC_%(PL)s%(PL)s_white_noise_white_Green_snapshot_1_00.dat'

	chosenCols = [0, 1, 3]

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(vQ2Labels)+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	for vQ2 in vQ2Labels:
		if vQ2=='100_000':
			filename = dndn_white_filenames % {'PL': particleLabels[particleIndex]}
		else:
			filename = dndn_color_filenames % {'PL': particleLabels[particleIndex], 'vQ2': vQ2}
		
		columnToPlot = 1
		SCstring = '_withSC'
		if subtractSelfCorrelations:
			columnToPlot = 2
			SCstring = '_withoutSC'
		
		colorVal = scalarMap.to_rgba(len(vQ2Labels)+1-idx-1)
		#colorVal = scalarMap.to_rgba(idx)

		# read in file
		data = loadtxt(filename, usecols=tuple(chosenCols))

		#print 'Loading data file...'
		ax.plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw, label=vQ2LineLabels[vQ2])
		idx+=1
	
	
	ax.axhline(0.0, color='black', linewidth=1)
	ax.tick_params(length=5, width=1.5, top='off', right='off')
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$B_{%(CBF)s}$' % {'CBF': CBFIndices[particleIndex]}, fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize + 5})

	#outfilename = 'results/twoPC_%(PL)s%(PL)s_vQ2_comparison' % {'PL': particleLabels[particleIndex]} + SCstring + '.pdf'
	#plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	#print 'Saved to', outfilename



#################################################################
def generate_plot_vQ2_comparison_TOGETHER():
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(3, 1, figsize=(6,10))
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0

	subtractSelfCorrelations = True

	dndn_color_filenames = 'results/twoPC_%(PL)s%(PL)s_color_noise_color_Green_vQ2_%(vQ2)s_snapshot_1_00.dat'
	dndn_white_filenames = 'results/twoPC_%(PL)s%(PL)s_white_noise_white_Green_snapshot_1_00.dat'

	chosenCols = [0, 1, 3]

	idx=0
	for pidx in xrange(3):
		for vQ2 in vQ2Labels:
			if vQ2=='100_000':
				filename = dndn_white_filenames % {'PL': particleLabels[pidx]}
			else:
				filename = dndn_color_filenames % {'PL': particleLabels[pidx], 'vQ2': vQ2}
		
			columnToPlot = 1
			SCstring = '_withSC'
			if subtractSelfCorrelations:
				columnToPlot = 2
				SCstring = '_withoutSC'
		
			colorVal = colorVals[vQ2]

			# read in file
			data = loadtxt(filename, usecols=tuple(chosenCols))

			#print 'Loading data file...'
			axs[pidx].plot(data[:,0], data[:,columnToPlot], color=colorVal,\
                           linestyle=vQ2LineStyles[vQ2],\
                           linewidth=lineWidthVals[vQ2],\
                           label=vQ2LineLabels[vQ2])
			axs[pidx].set_ylim(bottom=-0.001)
			idx+=1
	
	
		axs[pidx].axhline(0.0, color='black', linewidth=1)
		axs[pidx].tick_params(length=5, width=1.5, top='off', right='off')
	
		axs[pidx].set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
		axs[pidx].set_ylabel(r'$B_{%(CBF)s}$' % {'CBF': CBFIndices[pidx]}, fontsize = labelsize + 10)
	
	
	axs[1].legend(loc=0, ncol=1, prop={'size': plotfontsize + 5}, frameon=False)
	for pidx in xrange(3):
		if pidx!=2:
			plt.setp( axs[pidx].get_xticklabels(), visible=False)
		(ymin, ymax) = axs[pidx].get_ylim()
		axs[pidx].set_ylim(top=0.99*ymax)

	#plt.show(block=False)
	#outfilename = 'results/twoPC_grid_vQ2_comparison' + SCstring + '.pdf'
	outfilename = 'twoPC_grid_vQ2_comparison' + SCstring + '_PRCAPPROVED.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename






#################################################################
def generate_plot_snapshots(particleIndex):
	# set-up
	plotfontsize = 12
	fig, axs = plt.subplots(2, 2, figsize=(15, 10))
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	twoPC_color_filenames = 'results/twoPC_%(PL)s%(PL)s_color_noise_color_Green_vQ2_%(vQ2)s_snapshot_%(frac)s.dat'
	twoPC_white_filenames = 'results/twoPC_%(PL)s%(PL)s_white_noise_white_Green_snapshot_%(frac)s.dat'

	subtractSelfCorrelations = True

	cm = plt.get_cmap('gist_rainbow') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(snapshotLabels))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	panelIndex = 0
	for vQ2 in vQ2Labels:
		rowIndex = panelIndex / 2
		columnIndex = panelIndex % 2
		vQ2 = vQ2Labels[panelIndex]
		idx=0
		for fraction in snapshotLabels:
			if vQ2=='100_000':
				filename = twoPC_white_filenames % {'PL': particleLabels[particleIndex], 'frac': fraction}
				plotLabel = vQ2LineLabels[panelIndex]
			else:
				filename = twoPC_color_filenames % {'PL': particleLabels[particleIndex], 'vQ2': vQ2, 'frac': fraction}
				plotLabel = r'$v_Q^2 = %(vQ2)s$' % {'vQ2': vQ2LineLabels[panelIndex]}
		
			columnToPlot = 1
			SCstring = '_withSC'
			if subtractSelfCorrelations:
				columnToPlot = 2
				SCstring = '_withoutSC'
			chosenCols = [0, 1, 3]
		
			colorVal = scalarMap.to_rgba(len(snapshotLabels)-idx-1)
			#colorVal = scalarMap.to_rgba(idx)

			# read in file
			data = loadtxt(filename, usecols=tuple(chosenCols))

			#print 'Loading data file...'
			axs[rowIndex,columnIndex].plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw, label=r'$f=$' + fraction.replace('_','.'))
			axs[rowIndex,columnIndex].set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
			axs[rowIndex,columnIndex].set_ylabel(r'$B_{%(CBF)s}(\tau)$' % {'CBF': CBFIndices[particleIndex]}, fontsize = labelsize + 10)
			if columnIndex==1:
				axs[rowIndex,columnIndex].yaxis.tick_right()
				axs[rowIndex,columnIndex].yaxis.set_label_position("right")
			if rowIndex==0:
				plt.setp( axs[rowIndex,columnIndex].get_xticklabels(), visible=False)
			if panelIndex==0:
				axs[rowIndex,columnIndex].legend(loc=0, ncol=2, prop={'size': plotfontsize+5})
			plt.text(0.8, 0.2, plotLabel, horizontalalignment='center', verticalalignment='center', transform = axs[rowIndex,columnIndex].transAxes, size=30)
	
			idx+=1
			
		panelIndex += 1

	for rowIndex in xrange(2):
		for columnIndex in xrange(2):
			(xmin, xmax) = axs[rowIndex,columnIndex].get_xlim()
			(ymin, ymax) = axs[rowIndex,columnIndex].get_ylim()
			axs[rowIndex,columnIndex].set_xlim(right=0.99*xmax)
			axs[rowIndex,columnIndex].set_ylim(top=0.99*ymax)

	
	#axs[rowIndex,columnIndex].axhline(0.0, color='black', linewidth=1)
	
	outfilename = 'results/twoPC_%(PL)s%(PL)s_vQ2_grid_snapshots' % {'PL': particleLabels[particleIndex]} + SCstring + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename



#############################################

if __name__ == "__main__":
	#plt.ion()
	#for pidx in xrange(1):
	#	generate_plot_vQ2_comparison(pidx)
	#generate_plot_snapshots(0)
	generate_plot_vQ2_comparison_TOGETHER()
	#pause()
	print 'Finished all.'

# End of file
