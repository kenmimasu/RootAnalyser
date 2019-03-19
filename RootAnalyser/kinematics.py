import numpy as np
import copy
from operator import attrgetter
################################################################################
# Useful functions  

def sort(particle_list, attr):
    '''
    Takes a list of particles and returns another list sorted by specified 
    attribute
    '''
    return sorted(particle_list, key=attrgetter(attr),reverse=True)

def pt_sort(particle_list):
    '''Takes a list of particles and returns another list sorted by pT'''
    return sorted(particle_list, key=attrgetter('pt'),reverse=True)
    
def ee_sort(particle_list):
    '''Takes a list of particles and returns another list sorted by Energy'''
    return sorted(particle_list, key=attrgetter('ee'),reverse=True)
    
def pdot(part1,part2):
    '''Returns the scalar product of the 3-momenta of two particles'''
    return part1.px*part2.px + part1.py*part2.py + part1.pz*part2.pz

def fourmom(*particles):
    '''Takes a list of particles and returns the four momentum of the system'''
    ee, px, py, pz = 0., 0., 0., 0.
    for p in particles:
        ee += p.ee
        px += p.px
        py += p.py
        pz += p.pz
    return ee, px, py, pz 
        
def Minv(*particles):  # Minv(part1,part2,...)
    '''Takes a list of particles and returns the invariant mass of the system'''
    ee, px, py, pz = fourmom(*particles)
    return np.sqrt(ee**2 - px**2 - py**2 - pz**2)
    
def rap(*particles):
    '''Takes a list of particles and returns the rapidity of the system'''
    ee, px, py, pz = fourmom(*particles)
    return np.log( (ee + pz)/(ee - pz) )/2.
    
def costheta(*particles):
    '''
    Returns cosine of polar angle of particle system w.r.t z-direction.
    '''
    _, px, py, pz = fourmom(*particles)
    modp = np.sqrt(px**2+py**2+pz**2)
    costheta = pz/modp if modp!=0. else np.nan 
    return costheta

def costhetastar(*particles):
    '''
    Returns cosine of polar angle of the first particle w.r.t direction defined by the 
    boost direction of the whole system.
    '''
    _, px, py, pz = fourmom(*particles)
    p1 = particles[0]
    modp = np.sqrt(p1.px**2+p1.py**2+p1.pz**2)
    costheta = p1.pz/modp*pz/abs(pz) if (modp!=0.) else np.nan 
    return costheta


def psrap(*particles):
    '''Takes a list of particles and returns the pseudo-rapidity of the system'''
    costh = costheta(*particles)
    return np.log( (1 + costh)/(1 - costh) )/2.    
    
def pT(*particles):# pT(part1,part2,...)
    '''Takes a list of particles and returns the pT of the system'''
    px, py = 0., 0.
    for p in particles:
        px += p.px
        py += p.py
    return np.sqrt(px**2 + py**2) 
    

def HT(*particles):# HT(part1,part2,...)
    '''Takes a list of particles and returns the HT (scalar sum of pTs) of the system'''
    return sum([p.pt for p in particles])  
    
def phi(*particles):# phi(part1,part2,...)
    '''Takes a list of particles and returns the azimuthal angle of the 3-momentum of the system'''
    px, py = 0.,0.
    for p in particles:
        px += p.px; py += p.py
    phi = np.arctan(py/px) if px > 0 else np.arctan(py/px) + np.pi
    return phi if phi > 0 else phi + 2*np.pi
    
def dphi(*args,**kwargs):# dphi(part1,part2,...,ref_phi=XXX)
    '''Takes a list of particles and returns the azimuthal angle of the 3-momentum of the system optionally w.r.t keyword argument ref_phi'''
    ref_phi = kwargs.get('ref_phi',0.) # look for ref_phi keyword argument, if not use zero
    phi_tot = phi(*args)
    dphi = abs(phi_tot-ref_phi)
    return dphi if dphi < np.pi else 2.*np.pi-dphi

def dphi_0_2pi(*args,**kwargs):# dphi(part1,part2,...,ref_phi=XXX)
    '''Takes a list of particles and returns the azimuthal angle of the 3-momentum of the system optionally w.r.t keyword argument ref_phi'''
    ref_phi = kwargs.get('ref_phi',0.) # look for ref_phi keyword argument, if not use zero
    phi_tot = phi(*args)
    return abs(phi_tot-ref_phi)

def dtheta(part1,part2,degrees=False):
    '''Returns the opening angle between two particles' 3-momenta'''
    costheta = pdot(part1,part2)/part1.modp/part2.modp
    theta = np.arccos(costheta)
    return theta*180./np.pi if degrees else theta
    
def deltaR(part1, part2):
    '''Returns the spatial separation dR = sqrt(dEta^2 + dPhi^2) of two particles'''
    eta1, phi1 = part1.eta, part1.phi
    eta2, phi2 = part2.eta, part2.phi
    deta, dphi = abs(eta1-eta2), abs(phi1-phi2)
    if dphi > np.pi: dphi = 2*np.pi - dphi
    return np.sqrt( (deta)**2 + (dphi)**2  )

def delta(part1, part2, attr, absval=True):
    '''
    Returns the difference of attr data member of two particles
    '''
    a1, a2 = getattr(part1, attr), getattr(part2, attr)    
    return (a1-a2) if not absval else abs(a1-a2)
    
def met_px_py(MET,MET_phi):
    '''Translates missing energy and its azimuthal angle to Cartesian components'''
    return MET*np.cos(MET_phi), MET*np.sin(MET_phi)

def zboost(particle,rapidity=None,beta=None):
    '''return a new particle instance boosted in the z direction'''
    newpart = copy.deepcopy(particle)
    return newpart.zboost(rapidity=rapidity,beta=beta)

def CoM(*particles):
    newparts = []
    rapidity = rap(*particles)
    for p in particles:
        newparts.append(zboost(p,rapidity=-rapidity))
    return tuple(newparts)
    
    
    

################################################################################
