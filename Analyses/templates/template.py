########################################################################
# Selector for LHEF/LHCO formatted root files used for parton/
# detector level studies.
# Requires pyROOT, numpy and selector_tools.py
# Python version 2.7 to be safe 
from ROOT import TPySelector,TH1F
import ROOT
import random, os
import numpy as np
from RootAnalyser.event import Event
from RootAnalyser.root import *
from RootAnalyser.kinematics import *
from RootAnalyser.readers import read_tree
from itertools import product, combinations
########################################################################
# Template Selector
########################################################################
avail_objects='''
########################################################################
LHCO, LHEF selector
Available data members of Event instance, referenced by e.g. evt.photons
Containers:
    leptons - List of all lepton objects
    electrons, muons, taus, photons - Lists of each particle object
    jets - List of all jets
    bjets, ljets - Lists of b-tagged and non b-tagged jets
    exotics - Lists of all particles not identified above 
             (ONLY for LHE event format & they also contribute to MET)
Counters:
    npho, nlep, nele, nmu, ntau, njet, nljet, nbjet, nexo -
    numbers of each particle type i.e. length of corresponding list
Global event info:
    ht_tot - scalar sum of pT of all visible particles
    ht_jet - scalar sum of pT of all jets
    ee_jet - sum of all jet energies
    weight - evt weight, information only present in LHEF
Missing energy:
    MET - Missing transverse energy
    MET_phi - azimuthal direction of MET
########################################################################
'''
########################################################################
# Define your acceptance for different particles
acceptance={
    'pt_gam_min' : 20.,
    'pt_ele_min' : 20.,
    'pt_mu_min'  : 20.,
    'pt_jet_min' : 25.,
    'eta_gam_max': 9999999.,
    'eta_ele_max': 2.47,
    'eta_mu_max' : 2.5,
    'eta_tau_max': 9999999.,
    'eta_jet_max': 9999999.}
################################################################################
name='test'
################################################################################
# here you can define some global variables like cut values etc.
global_1 = 10.
global_2 = 20.
cuts = (100.,200.,300.,400.)
################################################################################
# write out histograms or not?
write_histos=True
# where to write out plot data
# plot_dir='/PATH/TO/PLOT/DIRECTORY/'
plot_dir = '/Users/Ken/GoogleDrive/Work/RootAnalyser/bums'
if not os.path.exists(plot_dir):
    print "plot_dir path, '{}' doesn't exist!".format(plot_dir)
    raise AssertionError
################################################################################
# Tuple of 1D histogram init arguments
# format: (NAME, TITLE, NBINS, XMIN, XMAX)
histargs = (
            ('hist1','Histogram I',50,0.,600.),
            ('hist2','Histogram II',50,0.,200.),
            )
################################################################################
class MyPySelector( TPySelector ):
    def Begin( self ):
        return 1
    def SlaveBegin( self, tree ):
        ############################# DON'T MODIFY #############################
        print '########################################################################'
        self.nEvents = tree.GetEntries()
        print 'Processing {} events...'.format(self.nEvents)
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        self.tot_wgt = 0.
        ########################################################################
        # here you can add stuff that 'belongs' to the analyser by assigning
        # data members to 'self', they can be retrieved inside the Process()
        # function below and referenced at the end of the analysis in 
        # Terminate().
        self.myvar = 20.
        self.npass = 0
        self.nfail = 0
        return 1
        
    @try_except
    def Process( self, entry ):
        ############################# DON'T MODIFY #############################
        tree = self.fChain
        if tree.GetEntry( entry ) <= 0:
            return 0
        ########################################################################
        # Read tree and fill info
        evt = read_tree(tree, acceptance=acceptance) # returns selector_tools.Event instance
        evt_wgt = evt.weight
        self.tot_wgt += evt_wgt
        # uncomment this is you want to kill after analysing 10 events
        # if entry ==10: self.Abort('ABORT!') 
        ########################################################################
        # Begin your analysis here!
        # have a browse through kinematics.py for what functions are available
        # Filling a histogram:
        if evt.njets >=2:
            jets = pt_sort(evt.jets) # sort jets according to pT
            
            # fill histograms referencing by NAME
            eta1, eta2 = jet[0].eta, jet[1].eta
            self.histos['hist1'].Fill(eta1*eta2) 
            
            mjj = Minv(jet1,jet2) # jet invariant mass
            self.histos['hist2'].Fill(mjj)
            
            # checking the efficiency of a cut
            if mjj > 100.:
                self.npass+=1
            else:
                self.nfail+=1
        
        # End analysis
        ########################################################################
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
        ############################# DON'T MODIFY #############################
        fname = self.fChain.GetCurrentFile().GetName()
        file_name= os.path.basename(fname).replace('.root','_analysis_cuts{}.root'.format(tag))
        cur_dir = ROOT.gDirectory
        out_file = ROOT.TFile(file_name,'RECREATE','test')
        out_file.cd()
        for h in histargs:
            hist= self.histos[h[0]]
            hist.Write()
            if write_histos: write_hist(hist,'{}/{}_{}.dat'.format(plot_dir,file_name.replace('.root',''),hist.GetName()))
        cur_dir.cd()
        out_file.Close()
        del out_file
        print ''
        print 'wrote root output to {}'.format(os.path.abspath(file_name))
        if write_histos: print 'wrote histogram data to {}'.format(plot_dir)
        ########################################################################
        # here you can do some diagnostics by printing things out
        print '########################################################################'
        print 'Sum of event weights: ', self.tot_wgt
        print 'Cut Efficiency:'
        print self.npass, ' events out of ', self.npass+self.nfail, 'passed the cut.'
        print '########################################################################'
        return 1
################################################################################

