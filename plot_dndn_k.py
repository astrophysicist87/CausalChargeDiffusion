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
panelLabels = ['(a)', '(b)', '(c)', '(d)']
panelCounter = 0

CBFIndices = ['\pi\pi','p \\bar{p}','K K']
particleLabels = ['pi','p','K']
#vQ2Labels = ['0_333', '1_000', '10_000', '100_000']
vQ2Labels = ['0_333', '1_000', '100_000']
#vQ2PlotLabels = {'0_333' : r'$v_Q^2=1/3$', '1_000' : r'$v_Q^2=1$', '10_000' : r'$v_Q^2=10$' , '100_000' : 'White'}
vQ2PlotLabels = {'0_333' : r'$v_Q^2=1/3$', '1_000' : r'$v_Q^2=1$', '100_000' : 'White'}
snapshotLabels = ['0_05', '0_10', '0_15', '0_25', '0_50', '1_00']
colorVals = ['red', 'blue', 'green', 'purple']
arrowheadLocs = {'0_333' : [0.6,0.55], '1_000' : [0.55,0.5], '10_000' : [0.5,0.4] , '100_000' : [0.5,0.4]}
initCurveLabelLocs = {'0_333' : [4.2,1.0], '1_000' : [3.4,2.0], '10_000' : [0.5,0.4] , '100_000' : [2.4,2.0]}

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
def plotCorrelatorSnapshots_TOGETHER():
	# set-up
	plotfontsize = 12
	#fig, ax = plt.subplots(1, 1)
	fig, axs = plt.subplots(3, 1, figsize=(8,16))
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	subtractSelfCorrelations = True
	dndn_color_filenames = 'results/dndn_k_pipi_color_noise_color_Green_vQ2_%(vQ2)s_snapshot_%(frac)s.dat'
	dndn_white_filenames = 'results/dndn_k_pipi_white_noise_white_Green_snapshot_%(frac)s.dat'

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(snapshotLabels))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	pidx=0
	for vQ2 in vQ2Labels:
		idx=0
		for fraction in snapshotLabels:
			if fraction=='0_50':
				idx+=1
				continue
			if vQ2 == '100_000':
				filename = dndn_white_filenames % {'frac': fraction}
				#plotLabel = 'White'
			else:
				filename = dndn_color_filenames % {'vQ2': vQ2, 'frac': fraction}
				#plotLabel = r'$v_Q^2 = %(vQ2)s$' % {'vQ2': vQ2.replace('_','.')}
		
			columnToPlot = 1
			SCstring = '_withSC'
			if subtractSelfCorrelations:
				columnToPlot = 2
				SCstring = '_withoutSC'
			chosenCols = [0, 1, 2]
		
			#colorVal = scalarMap.to_rgba(len(snapshotLabels)-idx-1)
			colorVal = scalarMap.to_rgba(idx)

			# read in file
			#print filename
			data = loadtxt(filename, usecols=tuple(chosenCols))

			#print 'Loading data file...'
			axs[pidx].plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw, label=r'$f=$' + fraction.replace('_','.'))
			idx+=1
	
			arrowLength = 0.15
			arrowheadLoc = arrowheadLocs[vQ2]
			initCurveLabelLoc = initCurveLabelLocs[vQ2]
			axs[pidx].axhline(0.0, color='black', linewidth=1)

			axs[pidx].set_xlabel(r'$k$', fontsize = labelsize + 10)
			axs[pidx].set_ylabel(r'$\left< \delta \tilde n(k) \delta \tilde n(-k) \right>$', fontsize = labelsize + 10)
			axs[pidx].set_xlim(right=15.0)
			axs[pidx].tick_params(length=5, width=1.5, top='off', right='off')
			axs[pidx].annotate('', xy=(arrowheadLoc[0],arrowheadLoc[1]), xycoords='axes fraction',\
                               xytext=(arrowheadLoc[0],arrowheadLoc[1]+arrowLength),\
                               textcoords='axes fraction', arrowprops=dict(facecolor='black',\
                               width=1, headwidth=5, frac=0.15))
			plt.text(arrowheadLoc[0]+0.125, arrowheadLoc[1]+0.5*arrowLength,\
                     r'Increasing', horizontalalignment='center', verticalalignment='center',\
                     transform = axs[pidx].transAxes, size=15)
			plt.text(arrowheadLoc[0]+0.25, arrowheadLoc[1]+0.5*arrowLength,\
                     r'$\tau$', horizontalalignment='center', verticalalignment='center',\
                     transform = axs[pidx].transAxes, size=25)
			#axs[pidx].legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
			plt.text(0.7, 0.9, vQ2PlotLabels[vQ2], horizontalalignment='center',\
                     verticalalignment='center', transform = axs[pidx].transAxes, size=30)
			plt.text(initCurveLabelLoc[0], initCurveLabelLoc[1], r'$5\%$ of evolution',\
                     horizontalalignment='left', verticalalignment='center', transform = axs[pidx].transData, size=15)
		pidx+=1
	
	for pidx in xrange(3):
		if pidx!=2:
			plt.setp( axs[pidx].get_xticklabels(), visible=False)
		(ymin, ymax) = axs[pidx].get_ylim()
		axs[pidx].set_ylim(bottom=0.99*ymin,top=0.99*ymax)
		if pidx==2:
			axs[pidx].set_ylim(bottom=-0.1,top=3.2)

	#plt.show(block=False)
	outfilename = 'dndn_k_pipi_grid_vQ2_comparison' + SCstring + '_PRCAPPROVED.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename



#################################################################
def plotCorrelatorSnapshots(vQ2):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	subtractSelfCorrelations = True
	dndn_color_filenames = 'results/dndn_k_pipi_color_noise_color_Green_vQ2_%(vQ2)s_snapshot_%(frac)s.dat'
	dndn_white_filenames = 'results/dndn_k_pipi_white_noise_white_Green_snapshot_%(frac)s.dat'

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(snapshotLabels))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	for fraction in snapshotLabels:
		if fraction=='0_50':
			idx+=1
			continue
		if vQ2 == '100_000':
			filename = dndn_white_filenames % {'frac': fraction}
			#plotLabel = 'White'
		else:
			filename = dndn_color_filenames % {'vQ2': vQ2, 'frac': fraction}
			#plotLabel = r'$v_Q^2 = %(vQ2)s$' % {'vQ2': vQ2.replace('_','.')}
		
		columnToPlot = 1
		SCstring = '_withSC'
		if subtractSelfCorrelations:
			columnToPlot = 2
			SCstring = '_withoutSC'
		chosenCols = [0, 1, 2]
		
		#colorVal = scalarMap.to_rgba(len(snapshotLabels)-idx-1)
		colorVal = scalarMap.to_rgba(idx)

		# read in file
		#print filename
		data = loadtxt(filename, usecols=tuple(chosenCols))

		#print 'Loading data file...'
		ax.plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw, label=r'$f=$' + fraction.replace('_','.'))
		idx+=1
	
	arrowLength = 0.15
	arrowheadLoc = arrowheadLocs[vQ2]
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$k$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$\left< \delta \tilde n(k) \delta \tilde n(-k) \right>$', fontsize = labelsize + 10)
	ax.set_xlim(right=15.0)
	ax.tick_params(length=5, width=1.5, top='off', right='off')
	ax.annotate('', xy=(arrowheadLoc[0],arrowheadLoc[1]),  xycoords='axes fraction',  xytext=(arrowheadLoc[0],arrowheadLoc[1]+arrowLength), textcoords='axes fraction', arrowprops=dict(facecolor='black', width=1, headwidth=5, frac=0.15))
	plt.text(arrowheadLoc[0]+0.125, arrowheadLoc[1]+0.5*arrowLength, r'Increasing', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=15)
	plt.text(arrowheadLoc[0]+0.25, arrowheadLoc[1]+0.5*arrowLength, r'$\tau$', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=25)
	#ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	plt.text(0.5, 0.9, vQ2PlotLabels[vQ2], horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	
	#outfilename = 'dndn_k_pipi_vQ2_%(vQ2)s_snapshots' % {'vQ2': vQ2} + SCstring + '.pdf'
	#plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	#print 'Saved to', outfilename


#############################################

if __name__ == "__main__":
	#plt.ion()
	#for vQ2 in vQ2Labels:
	#	plotCorrelatorSnapshots(vQ2)
	plotCorrelatorSnapshots_TOGETHER()
	#pause()
	print 'Finished all.'

# End of file
