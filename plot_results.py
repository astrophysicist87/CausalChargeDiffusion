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
nDY = 501

panelLabels = ['(a)', '(b)']
panelCounter = 0

resultsDir='results/'

CBFIndices = ['\pi\pi','p \\bar{p}','K K']
particleLabels = ['pi','p','K']
#ExpDataFilenames = ['star_pK.dat']
#ExpDataFilenames = ['star_pipi.dat', 'star_ppbar.dat', 'star_KK.dat']
vQ2Labels = ['0_333', '1_000', '10_000', '100_000']
snapshotLabels = ['0_05', '0_10', '0_15', '0_25', '0_50', '1_00']

#lineColors = ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
#plotTitles = ['Lattice', r'$2\pi D_Q T = 0.5$', r'$2\pi D_Q T = 1.0$', r'$2\pi D_Q T = 1.5$']

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")


#################################################################
# Plot two-particle correlation snapshots
def plotTwoPCSnapshots(particleLabel, noiseType, GreenType, vQ2, subtractSelfCorrelations):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	twoPC_color_filenames = resultsDir + 'twoPC_%(PL)s%(PL)s_%(NT)s_noise_%(GT)s_Green_vQ2_%(vQ2)s_snapshot_%(frac)s.dat'
	twoPC_white_filenames = resultsDir + 'twoPC_%(PL)s%(PL)s_white_noise_white_Green_snapshot_%(frac)s.dat'

	cm = plt.get_cmap('gist_rainbow') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(snapshotLabels))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	for fraction in snapshotLabels:
		if noiseType=='white' and GreenType=='white':
			filename = twoPC_white_filenames % {'PL': particleLabel, 'frac': fraction}
		else:
			filename = twoPC_color_filenames % {'PL': particleLabel, 'NT': noiseType, 'GT': GreenType, 'vQ2': vQ2, 'frac': fraction}
		
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
		ax.plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw, label=r'$f=$' + fraction.replace('_','.'))
		idx+=1
	
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$B_{\pi\pi}(\tau = \tau_0 + f (\tau_f - \tau_0))$', fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	plt.title(noiseType + ' noise, ' + GreenType + ' Green function')
	text(0.8, 0.2, r'$v_Q^2 = %(vQ2)s$' % {'vQ2': vQ2.replace('_','.')}, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	
	#plt.show(block=False)
	outfilename = resultsDir + 'twoPC_%(PL)s%(PL)s_vQ2_%(vQ2)s_snapshots' % {'PL': particleLabel, 'vQ2': vQ2} + SCstring + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename



#################################################################
# Plot snapshots
def plotCorrelatorSnapshots(particleLabel, space, noiseType, GreenType, vQ2, subtractSelfCorrelations):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	dndn_color_filenames = resultsDir + 'dndn_%(space)s_%(PL)s%(PL)s_%(NT)s_noise_%(GT)s_Green_vQ2_%(vQ2)s_snapshot_%(frac)s.dat'
	dndn_white_filenames = resultsDir + 'dndn_%(space)s_%(PL)s%(PL)s_white_noise_white_Green_snapshot_%(frac)s.dat'

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(snapshotLabels))
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	for fraction in snapshotLabels:
		if noiseType=='white' and GreenType=='white':
			filename = dndn_white_filenames % {'PL': particleLabel, 'space': space, 'frac': fraction}
		else:
			filename = dndn_color_filenames % {'PL': particleLabel, 'space': space, 'NT': noiseType, 'GT': GreenType, 'vQ2': vQ2, 'frac': fraction}
		
		columnToPlot = 1
		SCstring = '_withSC'
		if subtractSelfCorrelations:
			columnToPlot = 2
			SCstring = '_withoutSC'
		chosenCols = [0, 1, 2]
		
		colorVal = scalarMap.to_rgba(len(snapshotLabels)-idx-1)
		#colorVal = scalarMap.to_rgba(idx)

		# read in file
		data = loadtxt(filename, usecols=tuple(chosenCols))

		#print 'Loading data file...'
		ax.plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw, label=r'$f=$' + fraction.replace('_','.'))
		idx+=1
	
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	spaceXLabel = r'$k$'
	spaceYLabel = r'$\left< \delta \tilde n(k) \delta \tilde n(-k) \right>$'
	if space == 'Dxi':
		spaceXLabel = r'$\Delta \xi$'
		spaceYLabel = r'$\left< \delta \tilde n(\Delta \xi) \delta \tilde n(0) \right>$'
	ax.set_xlabel(spaceXLabel, fontsize = labelsize + 10)
	ax.set_ylabel(spaceYLabel, fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	plt.title(noiseType + ' noise, ' + GreenType + ' Green function')
	text(0.15, 0.9, r'$v_Q^2 = %(vQ2)s$' % {'vQ2': vQ2.replace('_','.')}, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	
	#plt.show(block=False)
	outfilename = resultsDir + 'dndn_%(space)s_%(PL)s%(PL)s_vQ2_%(vQ2)s_snapshots' % {'PL': particleLabel, 'space': space, 'vQ2': vQ2} + SCstring + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename


#################################################################
# Plot chosen vQ2s
def plotChosenvQ2s(fileStem, subtractSelfCorrelations, my_y_label, chosenCols):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 2.0

	fraction = '1_00'

	dndn_color_filenames = resultsDir + fileStem + '_color_noise_color_Green_vQ2_%(vQ2)s_snapshot_%(frac)s.dat'
	dndn_white_filenames = resultsDir + fileStem + '_white_noise_white_Green_snapshot_%(frac)s.dat'

	cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=len(vQ2Labels)+1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

	idx=0
	for vQ2 in vQ2Labels:
		if vQ2=='100_000':
			filename = dndn_white_filenames % {'frac': fraction}
		else:
			filename = dndn_color_filenames % {'vQ2': vQ2, 'frac': fraction}
		
		columnToPlot = 1
		SCstring = '_withSC'
		if subtractSelfCorrelations:
			columnToPlot = 2
			SCstring = '_withoutSC'
		#chosenCols = [0, 1, 2]
		
		colorVal = scalarMap.to_rgba(len(vQ2Labels)+1-idx-1)
		#colorVal = scalarMap.to_rgba(idx)

		# read in file
		data = loadtxt(filename, usecols=tuple(chosenCols))

		#print 'Loading data file...'
		ax.plot(data[:,0], data[:,columnToPlot], color=colorVal, linestyle='-', linewidth=lw, label=r'$v_Q^2 = $' + vQ2.replace('_','.'))
		idx+=1
	
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta \xi$', fontsize = labelsize + 10)
	ax.set_ylabel(my_y_label, fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	plt.title(fileStem + ': vQ2 comparison')
	
	#plt.show(block=False)
	outfilename = resultsDir + fileStem + '_vQ2_comparison_' + SCstring + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename



#################################################################
# Plot self correlations
def plotSelfCorrelations():
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0
		
	# read in file
	data_vQ2_0_1 = loadtxt('dndn_Dxi_vQ2_0_1.out')
	data_vQ2_0_333 = loadtxt('dndn_Dxi_vQ2_0_333.out')
	data_vQ2_1_0 = loadtxt('dndn_Dxi_vQ2_1_0.out')
	data_vQ2_3_0 = loadtxt('dndn_Dxi_vQ2_3_0.out')

	#print 'Loading data file...'
	ax.plot(data_vQ2_0_1[:,0], data_vQ2_0_1[:,1], color='purple', linestyle='-', linewidth=lw, label=r'$v_Q^2 = 1/10$')
	ax.plot(data_vQ2_0_333[:,0], data_vQ2_0_333[:,1], color='blue', linestyle='-', linewidth=lw, label=r'$v_Q^2 = 1/3$')
	ax.plot(data_vQ2_1_0[:,0], data_vQ2_1_0[:,1], color='green', linestyle='-', linewidth=lw, label=r'$v_Q^2 = 1$')
	ax.plot(data_vQ2_3_0[:,0], data_vQ2_3_0[:,1], color='red', linestyle='-', linewidth=lw, label=r'$v_Q^2 = 3$')
	
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta \xi$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$\left< \delta n(\Delta \xi, \tau_f) \delta n (0, \tau_f) \right> _{\mathrm{self}}$', fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	plt.xlim(-0.099, 0.099) 
	#plt.title(fileStem + ': vQ2 comparison, subtractSC = ' + str(subtractSelfCorrelations))
	
	#plt.show(block=False)
	outfilename = 'SelfCorrelationsVsVQ2.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename



#################################################################
def generate_all_plots():
	plotTwoPCSnapshots('pi', 'white', 'white', '100', True)	#vQ2 value irrelevant for pure white noise
	plotTwoPCSnapshots('pi', 'color', 'color', '10_000', True)
	plotTwoPCSnapshots('pi', 'color', 'color', '1_000', True)
	plotTwoPCSnapshots('pi', 'color', 'color', '0_333', True)
	#plotTwoPCSnapshots('pi', 'white', 'white', '100', False)	#vQ2 value irrelevant for pure white noise
	#plotTwoPCSnapshots('pi', 'color', 'color', '0_333', False)
	#plotTwoPCSnapshots('pi', 'color', 'color', '1_000', False)
	#pause()
	plotCorrelatorSnapshots('pi', 'k', 'white', 'white', '100', True)	#vQ2 value irrelevant for pure white noise
	plotCorrelatorSnapshots('pi', 'k', 'color', 'color', '10_000', True)
	plotCorrelatorSnapshots('pi', 'k', 'color', 'color', '1_000', True)
	plotCorrelatorSnapshots('pi', 'k', 'color', 'color', '0_333', True)
	#plotCorrelatorSnapshots('pi', 'k', 'white', 'white', '100', False)	#vQ2 value irrelevant for pure white noise
	#plotCorrelatorSnapshots('pi', 'k', 'color', 'color', '1_000', False)
	#plotCorrelatorSnapshots('pi', 'k', 'color', 'color', '0_333', False)
	plotCorrelatorSnapshots('pi', 'Dxi', 'white', 'white', '100', True)
	plotCorrelatorSnapshots('pi', 'Dxi', 'color', 'color', '10_000', True)
	plotCorrelatorSnapshots('pi', 'Dxi', 'color', 'color', '1_000', True)
	plotCorrelatorSnapshots('pi', 'Dxi', 'color', 'color', '0_333', True)
	#pause()
	########
	fileStem = 'dndn_%(space)s_%(PL)s%(PL)s' % {'PL': 'pi', 'space': 'k'}
	plotChosenvQ2s(fileStem, True, r'$\left< \delta \tilde n(k) \delta \tilde n(-k) \right>$', [0, 1, 2])
	#fileStem = 'dndn_%(space)s_%(PL)s%(PL)s' % {'PL': 'p', 'space': 'k'}
	#plotChosenvQ2s(fileStem, True, r'$\left< \delta \tilde n(k) \delta \tilde n(-k) \right>$', [0, 1, 2])
	#fileStem = 'dndn_%(space)s_%(PL)s%(PL)s' % {'PL': 'K', 'space': 'k'}
	#plotChosenvQ2s(fileStem, True, r'$\left< \delta \tilde n(k) \delta \tilde n(-k) \right>$', [0, 1, 2])
	#pause()
	########
	plotChosenvQ2s('twoPC_%(PL)s%(PL)s' % {'PL': 'pi'}, True, r'$B_{\pi\pi}$', [0, 1, 3])
	#plotChosenvQ2s('twoPC_%(PL)s%(PL)s' % {'PL': 'p'}, True, r'$B_{p\bar p}$', [0, 1, 3])
	#plotChosenvQ2s('twoPC_%(PL)s%(PL)s' % {'PL': 'K'}, True, r'$B_{K\bar K}$', [0, 1, 3])
	#################
	#plotSelfCorrelations()
	#pause()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
