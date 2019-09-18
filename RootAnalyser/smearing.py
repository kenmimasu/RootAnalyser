import numpy as np
import copy
from operator import attrgetter
from scipy.stats import norm
from util import locate_2D
################################################################################
# Smearing functions for particle momenta  
# Parametrised functions obtained from fits in
# (b-)jets: https://arxiv.org/pdf/1706.04965.pdf (p34)
# taus: http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/PERF-2013-06/
# forward jets: https://indico.cern.ch/event/579660/contributions/2496130/attachments/1495285/2326310/20170719_JERC_CMS_BOOST.pdf (p3)
# def smear(particle, res):
#     '''Smear particle 4 momentum according to a Gaussian of width res.'''
#     # duplicate particle
#     smeared_particle = copy.deepcopy(particle)
#     # determine smearing factor
#     smear_factor =  norm.rvs(loc=1.,scale=res)
#     # rescale pt, mass
#     smeared_particle.pt *= smear_factor
#     smeared_particle.mass *= smear_factor
#     # reset four momentum
#     smeared_particle.set_four_momentum(self.pt,self.eta,self.phi,self.mass)
#
#     return smeared_particle

def smear_jet(*particles, **kwargs):
    '''Smear (b-)jet according to pt & eta dependent resolution.'''
    seed = kwargs.get('seed',None)
    result = []
    for p in particles:
        # get bjet resolution as a function of pt & eta
        pti, etai = p.pt, p.eta

        res = 0
        if (abs(etai)<=1.3):
            res = 0.215702 + -0.0266012*np.log(pti) + 1.62047e-05*pti
        elif (abs(etai)<=2.5):
            res = 0.25 + -0.0266012*np.log(pti) + 1.62047e-05*pti 
        elif (abs(etai)<=4.7):
            # fitted from low pileup curve of plot
            res = 0.0742 + 0.00357*np.log(pti) + 1.87/pti
        else:
            # do nothing
            result.append(p)
            continue

        result.append( p.smeared(res, seed=seed) )
    
    return result[0] if (len(result)==1) else tuple(result) 

def smear_tau_hadr(*particles, **kwargs):
    seed = kwargs.get('seed',None)
    #### Data for hadronic tau
    _tau_had_eta_bins = [0., 0.8, 1.3, 1.6, 2.4]

    _tau_had_pt_bins = [15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 
                        85.0, 100.0, 150.0, 200.0, 300.0, 7000.0]

    _tau_had_eff_bins = [ 
    # eta -->                         pt |
    #                                    V   
    #   bin 1, bin 2, bin 3, bin 4
        22.  , 24.5 , 24.5 , 24.5 , # bin 1
        16.5 , 20.75, 22.5 , 20.9 , # bin 2
        13.15, 17.  , 20.  , 21.  , # bin 3
        11.5 , 15.  , 19.  , 19.  , # bin 4
        10.5 , 13.5 , 17.  , 17.  , # bin 5
        9.5  , 13.5 , 14.5 , 14.5 , # bin 6
        8.5  , 10.75, 12.  , 13.  , # bin 7
        7.   , 9.   , 10.5 , 12.  , # bin 8
        6.5  , 8.   , 9.25 , 10.  , # bin 9
        6.   , 7.   , 7.3  , 8.5  , # bin 10
        5.   , 6.   , 6.   , 8.   , # bin 11
        4.75 , 6.   , 7.   , 5.     # bin 12
    ]
    
    result = []
    for p in particles:
        # locate pt, eta in efficiency bin
        pti, etai = p.pt, abs(p.eta)
        ieff = locate_2D(etai, pti, _tau_had_eta_bins, _tau_had_pt_bins)
        try:
            res= _tau_had_eff_bins[ieff]/100.
            # smear particle
            result.append( p.smeared(res, seed=seed) )
        except IndexError:
            # do nothing
            result.append(p)
        
    return result[0] if (len(result)==1) else tuple(result) 

def smear_tau_elec(*particles, **kwargs):
    seed = kwargs.get('seed',None)
    result = []
    for p in particles:
        pti, etai = p.pt, p.eta
        if (abs(etai) > 2.47):
            # do nothing
            result.append(p)
            continue            
        
        par0=0.235186
        par1=0.297411
        par2=0.0137397*(1.+etai*0.07)
        
        res = np.sqrt( (par0/pti)**2 + (par1/pti)**2 + par2**2 ) 
        
        # smear particle
        result.append( p.smeared(res, seed=seed) )
    
    return result[0] if (len(result)==1) else tuple(result) 
    

def smear_tau_muon(*particles, **kwargs):
    seed = kwargs.get('seed', None)
    result = []
    for p in particles:
        pti, etai = p.pt, p.eta
        # get resolution as a function of pt for different eta bins
        if (abs(etai) < 1.05):
            par0, par1 = 0.0158746, 0.000393809
        elif ( 1.05 <=  (abs(etai)) and (abs(etai) < 1.7)):
            par0, par1 = 0.0239047, 0.000632637
        elif ( 1.7 <=  (abs(etai)) and (abs(etai) < 2.4)):
            par0, par1 =0.0330916, 0.000933841
        else:
            # do nothing
            result.append(p)
            continue
        
        res = np.sqrt( par0**2 +(par1*pti)**2  )
        
        # smear particle
        result.append( p.smeared(res, seed=seed) )
    
    return result[0] if (len(result)==1) else tuple(result) 

if __name__=='__main__':
    _tau_had_eta_bins = [0., 0.8, 1.3, 1.6, 2.4]

    _tau_had_pt_bins = [15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 
                        85.0, 100.0, 150.0, 200.0, 300.0, 7000.0]

    _tau_had_eff_bins = [ 
    # eta -->                         pt |
    #                                    V   
    #   bin 1, bin 2, bin 3, bin 4
        22.  , 24.5 , 24.5 , 24.5 , # bin 1
        16.5 , 20.75, 22.5 , 20.9 , # bin 2
        13.15, 17.  , 20.  , 21.  , # bin 3
        11.5 , 15.  , 19.  , 19.  , # bin 4
        10.5 , 13.5 , 17.  , 17.  , # bin 5
        9.5  , 13.5 , 14.5 , 14.5 , # bin 6
        8.5  , 10.75, 12.  , 13.  , # bin 7
        7.   , 9.   , 10.5 , 12.  , # bin 8
        6.5  , 8.   , 9.25 , 10.  , # bin 9
        6.   , 7.   , 7.3  , 8.5  , # bin 10
        5.   , 6.   , 6.   , 8.   , # bin 11
        4.75 , 6.   , 7.   , 5.     # bin 12
    ]
    
    pti, etai = 50., abs(-0.01)
    ieff = locate_2D(etai, pti, _tau_had_eta_bins, _tau_had_pt_bins)
    try:
        res= _tau_had_eff_bins[ieff]/100.
        # smear particle
        print res
    except IndexError:
        # do nothing
        pass