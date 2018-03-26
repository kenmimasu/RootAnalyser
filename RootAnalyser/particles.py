import numpy as np
import copy
from scipy.stats import norm
################################################################################
# Particle classes with alternate constructors from TRootXXXX (LHCO) and TRootLHEFParticle objects.
# See http://madgraph.phys.ucl.ac.be/Downloads/ExRootAnalysis/RootTreeDescription.html
# for description of data members in TRootLHEFParticle and TRootPhoton/Jet/Muon... classes
################################################################################
# Container classes
class Particle: # Base particle class

    def __init__(self, pt, eta, phi, mass=0., PID=None):
        self.pt  = pt # Transverse momentum
        self.phi = phi if phi > 0 else phi+2.*np.pi # Azimuthal angle
        self.mass = mass
        self.PID = PID
        self.eta = eta # Pseudorapidity
        # Function to set cartesian momenta and energy
        self.set_four_momentum(self.pt,self.eta,self.phi,self.mass) 
        self.modp = np.sqrt( self.px**2 + self.py**2 + self.pz**2 )
        self.costheta=self.pz/self.modp
        self.theta = np.arccos(self.costheta)
        
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        return instance(TRootLHEFParticle.PT, TRootLHEFParticle.Eta, 
                        TRootLHEFParticle.Phi, mass=TRootLHEFParticle.M)

    def set_four_momentum(self,pt,eta,phi,mass):
        self.modp = pt*np.cosh(eta)
        self.ee = np.sqrt(self.modp**2+mass**2)
        self.px = pt*np.cos(phi)
        self.py = pt*np.sin(phi)
        self.pz = np.sqrt(self.modp**2-pt**2)*(abs(eta)/eta if eta!=0. else 0.)
        
    def get_y(self,eta,pt,mass):
        modp =  pt*np.cosh(eta)
        ee = np.sqrt(modp**2+mass**2)
        pz = np.sqrt(modp**2-pt**2)*(abs(eta)/eta if eta!=0. else 0.)
        return np.log((ee+pz)/(ee-pz))/2.
        
    def set_pt_eta_phi(self,ee,px,py,pz):
        ct = pz/np.sqrt(px**2+py**2+pz**2)
        self.pt = np.sqrt(px**2+py**2) 
        self.eta = np.log( (1. + ct)/(1. - ct) )/2.
        self.phi = np.arctan(py/px) if px > 0 else np.arctan(py/px) + np.pi

    def zboost(self, rapidity=None, beta=None):
        if rapidity is not None:
            ee = np.cosh(rapidity)*self.ee + np.sinh(rapidity)*self.pz
            pz = np.sinh(rapidity)*self.ee + np.cosh(rapidity)*self.pz
        elif beta is not None:
            gamma = 1./np.sqrt(1.-beta**2)
            ee = gamma*(self.ee + beta*self.pz)
            pz = gamma*(self.pz + beta*self.ee)
        else:
            ee, pz = self.ee,self.pz
        
        self.ee, self.pz = ee, pz
        self.costheta=self.pz/self.modp
        self.theta = np.arccos(self.costheta)
        self.eta = np.log( (1. + self.costheta)/(1. - self.costheta) )/2.   
        return self
    
    def smeared(self, res, seed=None):
        '''Smear particle 4 momentum according to a Gaussian of width res.'''
        # duplicate particle
        smeared_particle = copy.deepcopy(self)
        # determine smearing factor
        smear_factor =  norm.rvs(loc=1.,scale=res, random_state=seed)
        # rescale pt, mass
        smeared_particle.pt *= smear_factor
        smeared_particle.mass *= smear_factor

        # reset four momentum
        smeared_particle.set_four_momentum(self.pt,self.eta,self.phi,self.mass)
    
        return smeared_particle
    
    def __repr__(self):
        return repr(self.__dict__)
##########
class Photon(Particle):
    def __init__(self, pt, eta, phi, hadem):
        Particle.__init__(self, pt, eta, phi)
        self.hadem = hadem # hadronic over electromagnetic energy
    @classmethod # Initialise from a TRoot particle object
    def TRoot(instance, TRootPhoton):
        pt, eta, phi = TRootPhoton.PT, TRootPhoton.Eta, TRootPhoton.Phi
        hadem = TRootPhoton.EhadOverEem
        return instance( pt, eta, phi, hadem)
    @classmethod # Initialise from a TRootLHEFParticle object
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        hadem = -999.
        return instance( pt, eta, phi, hadem)
    @classmethod # Initialise fake from a Jet object
    def FakeJet(instance, jet):
        pt, eta, phi, hadem = jet.pt, jet.eta, jet.phi, jet.hadem
        return instance( pt, eta, phi, hadem)
##########
class Lepton(Particle): # Base Lepton class
    def __init__(self, pt, eta, phi, ntrk, charge):
        Particle.__init__(self, pt, eta, phi)
        self.charge = charge # EM charge
        self.ntrk   = ntrk # Number of tracks
    @classmethod
    def TRoot(instance, TRootLepton):
        pt, eta, phi = TRootLepton.PT, TRootLepton.Eta, TRootLepton.Phi
        ntrk, charge = TRootLepton.Ntrk, TRootLepton.Charge
        return instance( pt, eta, phi, ntrk, charge)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        ntrk, charge = -999., -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, ntrk, charge)
##########
class Electron(Lepton):
    def __init__(self, pt, eta, phi, hadem, ntrk, charge):
        Lepton.__init__(self, pt, eta, phi, ntrk, charge)
        self.flavour = 'e'
        self.hadem = hadem
    @classmethod
    def TRoot(instance, TRootElectron):
        pt, eta, phi = TRootElectron.PT, TRootElectron.Eta, TRootElectron.Phi
        hadem, ntrk, charge = TRootElectron.EhadOverEem, TRootElectron.Ntrk, TRootElectron.Charge
        return instance( pt, eta, phi, hadem, ntrk, charge)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        hadem, ntrk, charge = -999.,-999., -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, hadem, ntrk, charge)
##########
class Muon(Lepton):
    def __init__(self, pt, eta, phi, charge, ntrk, pTiso, ETiso, jindx):
        Lepton.__init__(self, pt, eta, phi, ntrk, charge)
        self.flavour = 'mu'
        self.pTiso = pTiso # Track pT isolation variable
        self.ETiso = ETiso # Calorimeter ET isolation variable
        self.jindx = jindx
    @classmethod
    def TRoot(instance, TRootMuon):
        pt, eta, phi = TRootMuon.PT, TRootMuon.Eta, TRootMuon.Phi
        ntrk, charge = TRootMuon.Ntrk, TRootMuon.Charge
        pTiso, ETiso, jindx = TRootMuon.PTiso, TRootMuon.ETiso, TRootMuon.JetIndex
        return instance( pt, eta, phi, charge, ntrk, pTiso, ETiso, jindx)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        ntrk, charge = -999., -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        pTiso, ETiso, jindx = -999., -999., -999.
        return instance( pt, eta, phi, charge, ntrk, pTiso, ETiso, jindx)
##########
class Tau(Lepton):
    def __init__(self, pt, eta, phi, hadem, ntrk, charge):
        Lepton.__init__(self, pt, eta, phi, ntrk, charge)
        self.flavour = 'tau'
        self.hadem = hadem
    @classmethod
    def TRoot(instance, TRootTau):
        pt, eta, phi = TRootTau.PT, TRootTau.Eta, TRootTau.Phi
        hadem, ntrk, charge = TRootTau.EhadOverEem, TRootTau.Ntrk, TRootTau.Charge
        return instance( pt, eta, phi, hadem, ntrk, charge)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        hadem, ntrk, charge = -999.,-999., -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, hadem, ntrk, charge)
##########
class Jet(Particle):
    def __init__(self, pt, eta, phi, mass, hadem, ntrk, btag):
        Particle.__init__(self, pt, eta, phi, mass=mass)
        self.mass  = mass
        self.hadem = hadem
        self.ntrk  = ntrk
        self.btag  = bool(btag)
    @classmethod
    def TRoot(instance, TRootJet):
        pt, eta, phi, mass = TRootJet.PT, TRootJet.Eta, TRootJet.Phi, TRootJet.Mass
        hadem, ntrk, btag = TRootJet.EhadOverEem, TRootJet.Ntrk, TRootJet.BTag
        return instance(pt, eta, phi, mass, hadem, ntrk, btag)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi, mass = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi, TRootLHEFParticle.M
        hadem, ntrk, btag = -999.,-999., 1 if abs(TRootLHEFParticle.PID)==5 else 0
        return instance( pt, eta, phi, mass, hadem, ntrk, btag)
##########
class Top(Particle):
    def __init__(self, pt, eta, phi, mass=175.):
        Particle.__init__(self, pt, eta, phi, mass=mass)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi, mass = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi , TRootLHEFParticle.M
        return instance( pt, eta, phi, mass)
################################################################################
