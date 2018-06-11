import numpy as np
import copy
from scipy.stats import norm
################################################################################
# Particle classes with alternate constructors from TRootXXXX (LHCO) and TRootLHEFParticle objects.
# See http://madgraph.phys.ucl.ac.be/Downloads/ExRootAnalysis/RootTreeDescription.html
# for description of data members in TRootLHEFParticle and TRootPhoton/Jet/Muon... classes
# https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription for Delphes
################################################################################
# Container classes
class Particle: # Base particle class

    def __init__(self, pt, eta, phi, mass=0., PID=-999.):
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

    def zboost(self, rapidity=-999., beta=-999.):
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
    
    def smeared(self, res, seed=-999.):
        '''Smear particle 4 momentum according to a Gaussian of width res.'''
        # duplicate particle
        smeared_particle = copy.deepcopy(self)
        # determine smearing factor
        
        if type(seed) is str:
            if seed.lower()=='auto': 
                seed = abs(hash('{:.3f}'.format(self.pt))) % 2147483647
        
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
    def __init__(self, pt, eta, phi, hadem=-999., isol=-999.):
        Particle.__init__(self, pt, eta, phi)
        self.hadem = hadem # hadronic over electromagnetic energy
        self.isol = isol # isolation variable
    @classmethod # Initialise from a TRoot particle object
    def TRoot(instance, TRootPhoton):
        pt, eta, phi = TRootPhoton.PT, TRootPhoton.Eta, TRootPhoton.Phi
        hadem = TRootPhoton.EhadOverEem
        return instance( pt, eta, phi, hadem)
    @classmethod # Initialise from a Delphes particle object
    def Delphes(instance, DelphesPhoton):
        pt, eta, phi = DelphesPhoton.PT, DelphesPhoton.Eta, DelphesPhoton.Phi
        hadem = DelphesPhoton.EhadOverEem
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
    def __init__(self, pt, eta, phi, charge, ntrk=-999.):
        Particle.__init__(self, pt, eta, phi)
        self.charge = charge # EM charge
        self.ntrk   = ntrk # Number of tracks
    @classmethod
    def TRoot(instance, TRootLepton):
        pt, eta, phi = TRootLepton.PT, TRootLepton.Eta, TRootLepton.Phi
        ntrk, charge = TRootLepton.Ntrk, TRootLepton.Charge
        return instance( pt, eta, phi, charge, ntrk)
    @classmethod
    def Delphes(instance, DelphesLepton):
        pt, eta, phi = DelphesLepton.PT, DelphesLepton.Eta, DelphesLepton.Phi
        charge = DelphesLepton.Charge
        return instance( pt, eta, phi, charge)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        charge = -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, charge)
##########
class Electron(Lepton):
    def __init__(self, pt, eta, phi, charge, hadem=-999., ntrk=-999., isol=-999.):
        Lepton.__init__(self, pt, eta, phi, charge, ntrk=ntrk)
        self.flavour = 'e'
        self.hadem = hadem
        self.isol = isol
    @classmethod
    def TRoot(instance, TRootElectron):
        pt, eta, phi = TRootElectron.PT, TRootElectron.Eta, TRootElectron.Phi
        hadem, ntrk, charge = TRootElectron.EhadOverEem, TRootElectron.Ntrk, TRootElectron.Charge
        return instance( pt, eta, phi, charge, hadem=hadem, ntrk=ntrk)
    @classmethod
    def Delphes(instance, DelphesElectron):
        pt, eta, phi = DelphesElectron.PT, DelphesElectron.Eta, DelphesElectron.Phi
        hadem, charge = DelphesElectron.EhadOverEem, DelphesElectron.Charge
        return instance( pt, eta, phi, charge, hadem=hadem)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        hadem, ntrk, charge = -999.,-999., -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, hadem, ntrk, charge)
##########
class Muon(Lepton):
    def __init__(self, pt, eta, phi, charge, ntrk=-999., pTiso=-999., ETiso=-999., jindx=-999):
        Lepton.__init__(self, pt, eta, phi, charge, ntrk=ntrk)
        self.flavour = 'mu'
        self.pTiso = pTiso # Track pT isolation variable
        self.ETiso = ETiso # Calorimeter ET isolation variable
        self.jindx = jindx
    @classmethod
    def TRoot(instance, TRootMuon):
        pt, eta, phi = TRootMuon.PT, TRootMuon.Eta, TRootMuon.Phi
        ntrk, charge = TRootMuon.Ntrk, TRootMuon.Charge
        pTiso, ETiso, jindx = TRootMuon.PTiso, TRootMuon.ETiso, TRootMuon.JetIndex
        return instance( pt, eta, phi, charge, ntrk=ntrk, pTiso=pTiso, ETiso=ETiso, jindx=jindx)
    @classmethod
    def Delphes(instance, DelphesMuon):
        pt, eta, phi = DelphesMuon.PT, DelphesMuon.Eta, DelphesMuon.Phi
        charge = DelphesMuon.Charge
        return instance( pt, eta, phi, charge)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        charge = -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, charge)
##########
class Tau(Lepton):
    def __init__(self, pt, eta, phi, charge, hadem=-999., ntrk=-999.):
        Lepton.__init__(self, pt, eta, phi, charge, ntrk=ntrk)
        self.flavour = 'tau'
        self.hadem = hadem
    @classmethod
    def TRoot(instance, TRootTau):
        pt, eta, phi = TRootTau.PT, TRootTau.Eta, TRootTau.Phi
        hadem, ntrk, charge = TRootTau.EhadOverEem, TRootTau.Ntrk, TRootTau.Charge
        return instance( pt, eta, phi, hadem, ntrk, charge)
    @classmethod
    def Delphes(instance, DelphesTau):
        pt, eta, phi = DelphesTau.PT, DelphesTau.Eta, DelphesTau.Phi
        return instance( pt, eta, phi, 0.)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi
        hadem, ntrk, charge = -999.,-999., -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, hadem, ntrk, charge)
##########
class Jet(Particle):
    def __init__(self, pt, eta, phi, mass, btag, hadem=-999., ntrk=-999.):
        Particle.__init__(self, pt, eta, phi, mass=mass)
        self.mass  = mass
        self.hadem = hadem
        self.ntrk  = ntrk
        self.btag  = bool(btag)
    @classmethod
    def TRoot(instance, TRootJet):
        pt, eta, phi, mass = TRootJet.PT, TRootJet.Eta, TRootJet.Phi, TRootJet.Mass
        hadem, ntrk, btag = TRootJet.EhadOverEem, TRootJet.Ntrk, TRootJet.BTag
        return instance(pt, eta, phi, mass, btag, hadem=hadem, ntrk=ntrk)
    @classmethod
    def Delphes(instance, DelphesJet):
        pt, eta, phi, mass = DelphesJet.PT, DelphesJet.Eta, DelphesJet.Phi, DelphesJet.Mass
        hadem, ntrk, btag = DelphesJet.EhadOverEem, DelphesJet.NCharged, DelphesJet.BTag
        return instance(pt, eta, phi, mass, btag, hadem=hadem, ntrk=hadem)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi, mass = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi, TRootLHEFParticle.M
        btag = 1 if abs(TRootLHEFParticle.PID)==5 else 0
        return instance( pt, eta, phi, mass, btag)
##########
class Top(Particle):
    def __init__(self, pt, eta, phi, mass=175.):
        Particle.__init__(self, pt, eta, phi, mass=mass)
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        pt, eta, phi, mass = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi , TRootLHEFParticle.M
        return instance( pt, eta, phi, mass)
################################################################################
