########################################################################
# Template selector for LHEF/LHCO formatted root files used for parton/
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
avail_objects='''
########################################################################
LHCO,LHEF selector
Available data members of Event instance
Containers:
    leptons - List of all lepton objects
    electrons,muons,taus,photons - Lists of each particle object
    jets - List of all jets
    bjets,ljets - Lists of b-tagged and non b-tagged jets
Counters:
    npho, nlep, nele, nmu, ntau, njet, nljet, nbjet -
    numbers of each particle type i.e. length of corresponding list
Global event info:
    ht_tot - scalar sum of pT of all visible particles
    ht_jet - scalar sum of pT of all jets
    ee_jet - sum of all jet energies
Missing energy:
    MET - Missing transverse energy
    MET_phi - azimuthal direction of MET
########################################################################
'''
########################################################################
# Acceptance
acceptance={
    'pt_gam_min' : 0.,
    'pt_ele_min' : 20.,
    'pt_mu_min'  : 20.,
    'pt_tau_min' : 0.,
    'pt_jet_min' : 20.,
    'eta_gam_max': 9999999.,
    'eta_ele_max': 2.5,
    'eta_mu_max' : 2.7,
    'eta_tau_max': 9999999.,
    'eta_jet_max': 9999999.}
tag='test'
# Analysis cuts
pt_b1, pt_b2 = 45.,20.
pt_l1, pt_l2 = 40.,20.
Mll_min, Mll_max = 83.,99.
drbb_min, drbb_max = 0., 2.5
drll_min, drll_max = 0., 1.6
ptbb_min, ptbb_max = 0.,9999999.
ptb1_min, ptb1_max = 0.,9999999.
ptb2_min, ptb2_max = 0.,9999999.
ptll_min, ptll_max = 0.,9999999.
ptl1_min, ptl1_max = 0.,9999999.
ptl2_min, ptl2_max = 0.,9999999.
htbb_min, htbb_max = 120., 9999999.
htbbll_min, htbbll_max = 260., 9999999.
# signal region
Mbb_min, Mbb_max = 130.,190.
Mbbll_min, Mbbll_max = 340.,420.
########################################################################

########################################################################
write_histos=True
plot_dir='/Users/Ken/GoogleDrive/python/Analyses/plots' # where to write out plot data
assert os.path.exists(plot_dir),"'plot_dir path doesn't exist!"
########################################################################
# Tuple of histogram init arguments
histargs = (
            ('Mbb_2btag','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_sel','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_0l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_0l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_1l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_1l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_1l_MTW','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_MZ','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ1','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ2','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ3','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ4','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ5','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_0l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_0l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_1l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_1l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_1l_MTW','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_MZ','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ1','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ2','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ3','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ4','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ5','Inv. mass of two hardest b jets',12,10.,250.),
            )
########################################################################
class MyPySelector( TPySelector ):
    def Begin( self ):
        return 1
    def SlaveBegin( self, tree ):
        self.nEvents = tree.GetEntries()
        print 'Processing {} events...'.format(self.nEvents)
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        self.pTcuts = np.arange(0.,110.,10.)
        self.nPass = [0 for p in self.pTcuts]
        self.weights = [0. for p in self.pTcuts]
        self.tot_wgt = 0.
        return 1
        
    @try_except
    def Process( self, entry ):
        tree = self.fChain
        if tree.GetEntry( entry ) <= 0:
            return 0
        ################################################################
        # Read tree and fill info
        evt = read_tree(tree, acceptance=acceptance) # returns selector_tools.Event instance
        evt_wgt = evt.weight
        self.tot_wgt += evt_wgt
        if entry ==10: self.Abort('ABORT!')
        ################################################################
        # Begin Analysis
        two_btag = evt.nbjet == 2
        two_leptons = (evt.nele == 2 and evt.nmu < 1) or (evt.nmu == 2 and evt.nele < 1)
        one_lepton = (evt.nele == 1 and evt.nmu < 1) or (evt.nmu == 1 and evt.nele < 1)
        no_leptons = (evt.nele == 0 and evt.nmu == 0)
        if two_btag:
            # bjets
            sortbjets= pt_sort(evt.bjets)
            bjet1, bjet2 = sortbjets[0], sortbjets[1]
            Mbb = Minv(bjet1,bjet2)
            self.histos['Mbb_2btag'].Fill(Mbb)
            if evt.njet==2: # exactly 2 jet category
                no_additional_jets, one_additional_jet = True,False
            elif evt.nljet==1: # one additional jet category
                if evt.ljets[0].pt < min([b.pt for b in evt.bjets]):
                    no_additional_jets, one_additional_jet = False,True
                else: 
                    no_additional_jets, one_additional_jet = False,False
            else: no_additional_jets, one_additional_jet = False,False 
                     
            if (bjet1.pt > pt_b1) and (bjet2.pt > pt_b2):
                self.histos['Mbb_sel'].Fill(Mbb)
                if no_additional_jets: self.histos['Mbb_2j'].Fill(Mbb)
                if no_additional_jets and no_leptons:
                    self.histos['Mbb_2j_0l'].Fill(Mbb)
                    if evt.MET > 120.: self.histos['Mbb_2j_0l_MET'].Fill(Mbb) 
                if no_additional_jets and one_lepton:
                    self.histos['Mbb_2j_1l'].Fill(Mbb)
                    if evt.MET > 25.: 
                        self.histos['Mbb_2j_1l_MET'].Fill(Mbb)
                        lep = evt.leptons[0]
                        MET_px, MET_py=met_px_py(evt.MET,evt.MET_phi)
                        MTW = (evt.MET + lep.pt) - (lep.px + MET_px)**2- (lep.py + MET_py)**2
                        if MTW < 120.:self.histos['Mbb_2j_1l_MTW'].Fill(Mbb)
                if no_additional_jets and two_leptons:
                    self.histos['Mbb_2j_2l'].Fill(Mbb)
                    if evt.MET < 60.: 
                        self.histos['Mbb_2j_2l_MET'].Fill(Mbb)
                        sortleps=pt_sort(evt.leptons)
                        lep1, lep2 = sortleps[0], sortleps[1]
                        Mll =  Minv(lep1,lep2)
                        if (Mll > Mll_min) and (Mll < Mll_max):
                            self.histos['Mbb_2j_2l_MZ'].Fill(Mbb)
                            ptZ = pT(lep1,lep2)
                            dRbb = deltaR(bjet1,bjet2)
                            if ptZ <= 90.:
                                if dRbb >=0.7 and dRbb <=3.4:
                                     self.histos['Mbb_2j_2l_ptZ1'].Fill(Mbb)
                            if ptZ > 90. and ptZ <= 120.:
                                if dRbb >=0.7 and dRbb <=3.0:
                                     self.histos['Mbb_2j_2l_ptZ2'].Fill(Mbb)
                            if ptZ > 120. and ptZ <= 160.:
                                if dRbb >=0.7 and dRbb <=2.3:
                                     self.histos['Mbb_2j_2l_ptZ3'].Fill(Mbb)
                            if ptZ > 160. and ptZ <= 200.:
                                if dRbb >=0.7 and dRbb <=1.8:
                                     self.histos['Mbb_2j_2l_ptZ4'].Fill(Mbb)
                            if ptZ > 200.:
                                if dRbb < 1.4:
                                     self.histos['Mbb_2j_2l_ptZ5'].Fill(Mbb)
                                     
                if one_additional_jet: self.histos['Mbb_3j'].Fill(Mbb)
                if one_additional_jet and no_leptons:
                    self.histos['Mbb_3j_0l'].Fill(Mbb)
                    if evt.MET > 120.: self.histos['Mbb_3j_0l_MET'].Fill(Mbb)
                if one_additional_jet and one_lepton:
                    self.histos['Mbb_3j_1l'].Fill(Mbb)
                    if evt.MET > 25.: 
                        self.histos['Mbb_3j_1l_MET'].Fill(Mbb)
                        lep = evt.leptons[0]
                        MET_px, MET_py = met_px_py(evt.MET,evt.MET_phi)
                        MTW = (evt.MET + lep.pt) - (lep.px + MET_px)**2- (lep.py + MET_py)**2
                        if MTW < 120.:self.histos['Mbb_2j_1l_MTW'].Fill(Mbb)
                if one_additional_jet and two_leptons:
                    self.histos['Mbb_3j_2l'].Fill(Mbb)
                    if evt.MET < 60.: 
                        self.histos['Mbb_3j_2l_MET'].Fill(Mbb)
                        sortleps=pt_sort(evt.leptons)
                        lep1, lep2 = sortleps[0], sortleps[1]
                        Mll =  Minv(lep1,lep2)
                        if (Mll > Mll_min) and (Mll < Mll_max):
                            self.histos['Mbb_3j_2l_MZ'].Fill(Mbb)
                            ptZ = pT(lep1,lep2)
                            dRbb = deltaR(bjet1,bjet2)
                            if ptZ <= 90.:
                                if dRbb >=0.7 and dRbb <=3.4:
                                     self.histos['Mbb_3j_2l_ptZ1'].Fill(Mbb)
                            if ptZ > 90. and ptZ <= 120.:
                                if dRbb >=0.7 and dRbb <=3.0:
                                     self.histos['Mbb_3j_2l_ptZ2'].Fill(Mbb)
                            if ptZ > 120. and ptZ <= 160.:
                                if dRbb >=0.7 and dRbb <=2.3:
                                     self.histos['Mbb_3j_2l_ptZ3'].Fill(Mbb)
                            if ptZ > 160. and ptZ <= 200.:
                                if dRbb >=0.7 and dRbb <=1.8:
                                     self.histos['Mbb_3j_2l_ptZ4'].Fill(Mbb)
                            if ptZ > 200.:
                                if dRbb < 1.4:
                                     self.histos['Mbb_3j_2l_ptZ5'].Fill(Mbb) 
                                                                
        if evt.nmu+evt.nele!=evt.nlep or evt.nbjet+evt.nljet!=evt.njet: 
            self.Abort('ABORT!' )
            
        if not entry % 20000: print entry
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
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
        print '########################################################################'
        # print 'Total cross section (pb): ', self.tot_wgt
        # print 'pT Cut Efficiencies:'
        # # print self.pTcuts,self.nPass,self.weights
        # for cut,count,wgt in zip(self.pTcuts,self.nPass,self.weights):
        #     print '{}: {} (unweighted); {} (weighted) => {} pb'.format(cut, float(count)/float(self.nEvents),wgt/self.tot_wgt,wgt )
        return 1
####################################################################################################

