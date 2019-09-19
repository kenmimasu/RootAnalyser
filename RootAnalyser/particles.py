import numpy as np
import copy
from scipy.stats import norm
from collections import MutableSequence
from RootAnalyser.kinematics import fourmom
################################################################################
# Particle classes with alternate constructors from TRootXXXX (LHCO) and TRootLHEFParticle objects.
# See http://madgraph.phys.ucl.ac.be/Downloads/ExRootAnalysis/RootTreeDescription.html
# for description of data members in TRootLHEFParticle and TRootPhoton/Jet/Muon... classes
# https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/RootTreeDescription for Delphes
################################################################################
# Lorentz Classes
class FourVector(object):

    def __init__(self, x0, x1, x2, x3):
        self._entries = np.array([x0, x1, x2, x3])

    # abstract methods
    def __len__(self):
        return 4
    
    def __getitem__(self, index):
        return self._entries[index]
    
    def __setitem__(self, index, value):
        self._entries[index] = value
        return 
    
    def dot(self, fourvector):
        return (self[0]*fourvector[0] - self[1]*fourvector[1]
               -self[2]*fourvector[2] - self[3]*fourvector[3])
    
    def mod(self):
        return np.sqrt(self.dot(self))
    
    def _boost(self, bx, by, bz):
        be = np.sqrt(bx**2 + by**2 + bz**2)
        ga = 1./np.sqrt(1.-be**2)
        ff = ga**2/(1. + ga)
        
        LL = np.array([[ga    , -ga*bx      , -ga*by      , -ga*bz      ],
                       [-ga*bx, 1.+ ff*bx*bx, ff*bx*by    , ff*bx*bz    ],
                       [-ga*by, ff*by*bx    , 1.+ ff*by*by, ff*by*bz    ],
                       [-ga*bz, ff*bz*bx    , ff*bz*by    , 1.+ ff*bz*bz]])
        
        return np.dot(LL, self._entries)
    
    def boost(self, *args):
        self._entries = self.boost(*args)
    
    def boosted(self, *args):
        boosted_fourvector = self._boost(*args)
        return FourVector(*boosted_fourvector)

class FourMomentum(FourVector):
    
    @classmethod
    def init2(self, pt, eta, phi, mass=0.):
        return instance(*self.cartesian(pt, eta, phi, mass))
        
    @staticmethod
    def cartesian(pt, eta, phi, mass=0.):
        modp = pt*np.cosh(eta)
        ee = np.sqrt(modp**2+mass**2)
        px = pt*np.cos(phi)
        py = pt*np.sin(phi)
        pz = np.sqrt(modp**2-pt**2)*(abs(eta)/eta if eta!=0. else 0.)
        return np.array([ee, px, py, pz])
    
    def ee(self):
        return self[0]
        
    def px(self):
        return self[1]
    
    def py(self):
        return self[2]
        
    def pz(self):
        return self[3]
        
    def pt(self):
        return np.sqrt(self[1]**2+self[2]**2)
    
    def eta(self):
        modp = self.mod()
        if modp!=0.:
            ct = self.pz()/modp 
            return np.log( (1. + ct)/(1. - ct) )/2.
        else:
            return 0.
    
    def phi(self):
        pt = self.pt()
        px, py = self.px(), self.py()
        if pt!=0.:            
            if px != 0.:
                return np.arctan(py/px) if px > 0 else np.arctan(py/px) + np.pi
            else:
                return py/abs(py)*np.pi/2.
        else:
            return 0.
    
    def costheta(self):
        modp = self.mod()
        return 0. if modp==0. else self.pz()/modp
    
    def theta(self):
        return np.arccos(self.costheta())
    
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
        self.calc_four_momentum(self.pt,self.eta,self.phi,self.mass) 
        self.modp = np.sqrt( self.px**2 + self.py**2 + self.pz**2 )
        self.costheta=self.pz/self.modp
        self.theta = np.arccos(self.costheta)
        
    @classmethod
    def from_cartesian(instance, px, py, pz, mass=0., PID=-999.):
        pt, eta ,phi = instance.pt_eta_phi_from_cartesian(px, py, pz, mass)
        return instance(pt, eta, phi, mass=mass, PID=PID)
        
    @classmethod
    def LHEF(instance, TRootLHEFParticle):
        return instance(TRootLHEFParticle.PT, TRootLHEFParticle.Eta, 
                        TRootLHEFParticle.Phi, mass=TRootLHEFParticle.M)
                        
    @classmethod
    def system(instance, *parts):
        ee, px, py, pz = 0., 0., 0., 0.
        for p in parts:
            ee+= p.ee
            px+= p.px
            py+= p.py
            pz+= p.pz
        modpsq = px**2+py**2+pz**2
        modp = np.sqrt(modpsq)
        ct = pz/modp
        mass = np.sqrt(ee**2 - modpsq)
        pt = np.sqrt(px**2+py**2)
        eta = np.log( (1. + ct)/(1. - ct) )/2.
        phi = np.arctan(py/px) if px > 0 else np.arctan(py/px) + np.pi        
        return instance(pt, eta, phi, mass)

    def set_four_momentum(self, ee=None, px=None, py=None, pz=None):
        if ee is not None:
            self.ee = ee
        if px is not None:
            self.px = px
        if py is not None:
            self.py = py
        if pz is not None:
            self.pz = pz
        
        self.modp = np.sqrt( self.px**2 + self.py**2 + self.pz**2 )
        self.costheta = self.pz/self.modp
        self.theta = np.arccos(self.costheta)

        self.recalc_pt_eta_phi()
        

    def calc_four_momentum(self,pt,eta,phi,mass):
        self.modp = pt*np.cosh(eta)
        self.ee = np.sqrt(self.modp**2+mass**2)
        self.px = pt*np.cos(phi)
        self.py = pt*np.sin(phi)
        self.pz = np.sqrt(self.modp**2-pt**2)*(abs(eta)/eta if eta!=0. else 0.)
    
    def recalc_four_momentum(self):
        self.calc_four_momentum(self.pt,self.eta,self.phi,self.mass)

    def get_four_momentum(self):
        return [self.ee, self.px, self.py, self.pz]
    
    def get_y(self,eta,pt,mass):
        modp =  pt*np.cosh(eta)
        ee = np.sqrt(modp**2+mass**2)
        pz = np.sqrt(modp**2-pt**2)*(abs(eta)/eta if eta!=0. else 0.)
        return np.log((ee+pz)/(ee-pz))/2.
    
    def set_pt_eta_phi(self, pt=None, eta=None, phi=None, mass=None):
        if pt is not None:
            self.pt = pt
        if eta is not None:
            self.eta = eta
        if phi is not None:
            self.phi = phi
        if mass is not None:
            self.mass = mass
        
        self.recalc_four_momentum()
        self.modp = np.sqrt( self.px**2 + self.py**2 + self.pz**2 )
        self.costheta = self.pz/self.modp
        self.theta = np.arccos(self.costheta)

    @staticmethod 
    def pt_eta_phi_from_cartesian(px, py, pz, mass):
        modp = np.sqrt(px**2+py**2+pz**2)
        if modp!=0.:
            ct = pz/modp 
            eta = np.log( (1. + ct)/(1. - ct) )/2.
        else:
            eta=0.
            
        pt = np.sqrt(px**2+py**2) 
        if pt!=0.:            
            if px != 0.:
                phi = np.arctan(py/px) if px > 0 else np.arctan(py/px) + np.pi
            else:
                phi = py/abs(py)*np.pi/2.
        else:
            phi = 0.
        
        return pt, eta, phi
            
    def calc_pt_eta_phi(self,ee,px,py,pz):
        self.modp = np.sqrt(px**2+py**2+pz**2)
        if self.modp!=0.:
            ct = pz/self.modp 
            self.eta = np.log( (1. + ct)/(1. - ct) )/2.
        else:
            self.eta=0.
            
        self.pt = np.sqrt(px**2+py**2) 
        if self.pt!=0.:            
            if self.px != 0.:
                self.phi = np.arctan(self.py/self.px) if self.px > 0 else np.arctan(self.py/self.px) + np.pi
            else:
                self.phi = self.py/abs(self.py)*np.pi/2.
        else:
            self.phi = 0.
            
    def recalc_pt_eta_phi(self):
        self.calc_pt_eta_phi(self.ee, self.px, self.py, self.pz)

    def zboost(self, rapidity=None, beta=None):
        '''Boost particle in z-direction'''
        if rapidity is not None:
            ee = np.cosh(rapidity)*self.ee + np.sinh(rapidity)*self.pz
            pz = np.sinh(rapidity)*self.ee + np.cosh(rapidity)*self.pz
        elif beta is not None:
            gamma = 1./np.sqrt(1.-beta**2)
            ee = gamma*(self.ee + beta*self.pz)
            pz = gamma*(self.pz + beta*self.ee)
        else:
            # do nothing
            return self
        
        boost_particle = copy.deepcopy(self)
        boost_particle.set_four_momentum(ee=ee, pz=pz)

        return boost_particle
        
    def boost(self, ee, px, py, pz):
        '''boost particle momenta into the rest frame of 4 momentum
        (ee, px, py, pz ).'''
        pp = np.sqrt(px**2 + py**2 + pz**2)
        # mm = np.sqrt(ee**2 - pp**2)
        be = pp/ee
        bx, by, bz = px/ee, py/ee, pz/ee
        ga = 1./np.sqrt(1.-be**2)
        ff = ga**2/(1. + ga)
        
        LL = np.array([[ga    , -ga*bx      , -ga*by      , -ga*bz      ],
                       [-ga*bx, 1.+ ff*bx*bx, ff*bx*by    , ff*bx*bz    ],
                       [-ga*by, ff*by*bx    , 1.+ ff*by*by, ff*by*bz    ],
                       [-ga*bz, ff*bz*bx    , ff*bz*by    , 1.+ ff*bz*bz]])
        
        fv = np.array(self.get_four_momentum())
        ee1, px1, py1, pz1 = np.dot(LL, fv)
        
        boost_particle = copy.deepcopy(self)
        boost_particle.set_four_momentum(ee1, px1, py1, pz1)
        
        return boost_particle
    
    def boost_to_restframe(self, *parts):
        '''boost particle to the rest frame of a (system of) particle(s)'''
        return self.boost(*fourmom(*parts))
        
    def smeared(self, res, seed=999):
        '''Smear particle 4 momentum according to a Gaussian of width res.'''
        # determine smearing factor
        
        if type(seed) is str:
            if seed.lower()=='auto': 
                seed = abs(hash('{:.3f}'.format(self.pt))) % 2147483647
        
        smear_factor =  norm.rvs(loc=1.,scale=res, random_state=seed)
        
        # duplicate particle
        smeared_particle = copy.deepcopy(self)
        # rescale pt, mass
        smeared_particle.pt *= smear_factor
        smeared_particle.mass *= smear_factor

        # reset four momentum
        smeared_particle.recalc_four_momentum()
    
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
        charge = -TRootLHEFParticle.PID/abs(TRootLHEFParticle.PID)
        return instance( pt, eta, phi, charge)
        
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
    def __init__(self, pt, eta, phi, anti=False, mass=175.):
        self.anti=anti
        Particle.__init__(self, pt, eta, phi, mass=mass)
    @classmethod
    def LHEF(instance, TRootLHEFParticle, anti=False):
        pt, eta, phi, mass = TRootLHEFParticle.PT, TRootLHEFParticle.Eta, TRootLHEFParticle.Phi , TRootLHEFParticle.M
        return instance( pt, eta, phi, anti=False, mass=mass)
################################################################################
