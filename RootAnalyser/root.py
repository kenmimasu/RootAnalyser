import ROOT
from itertools import product
################################################################################
# Functions for writing histograms
def write_hist(hist, filename=None, normed=True):
    '''Function to write out ROOT histogram object to text file'''
    if filename is None: filename = '{}.dat'.format(hist.GetName())
    with open(filename,'w') as data:

        if isinstance(hist, ROOT.TH2):
            write_TH2(hist, data, normed=normed)
            
        elif isinstance(hist, ROOT.TH1):
            write_TH1(hist, data, normed=normed)

def write_TH1(hist, handle, normed=True):
    nevents = hist.GetEntries()
    integral = hist.Integral()
    norm = integral if (normed and integral > 0.) else 1.
    if normed:
        handle.write('# {} Entries, integral = 1}\n'.format(nevents))
    else:
        handle.write('# {} Entries, integral = {}\n'.format(nevents,integral))
    handle.write('# x\ty\tdy\n')
    for i in range(1, hist.GetNbinsX()+2):
        val = hist.GetBinContent(i)
        err = hist.GetBinError(i)
        x, y, dy = hist.GetBinCenter(i), val/norm, err/norm
        handle.write('{}\t{}\t{}\n'.format( x, y, dy ))

def write_TH2(hist, handle, normed=True):
    nevents = hist.GetEntries()
    integral = hist.Integral()
    norm = integral if (normed and integral > 0.) else 1.
    if normed:
        handle.write('# {} Entries, integral = 1}\n'.format(nevents))
    else:
        handle.write('# {} Entries, integral = {}\n'.format(nevents,integral))
    handle.write('# x\ty\tz\tdz\n')
    nxbins, nybins = hist.GetNbinsX()+2, hist.GetNbinsY()+2
    xaxis, yaxis = hist.GetXaxis(), hist.GetYaxis()
    for ix, iy in product(range(1, nxbins),range(1, nybins)):
        val = hist.GetBinContent(ix, iy)
        err = hist.GetBinError(ix,iy)
        x, y = xaxis.GetBinCenter(ix), yaxis.GetBinCenter(iy) 
        z, dz = val/norm, err/norm
        handle.write('{}\t{}\t{}\t{}\n'.format( x, y, z, dz ))
        
################################################################################
# Decorator to provide traceback when Process() loop raises an exception
def try_except(fn):
    """decorator for extra debugging output"""
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except:
            import traceback
            traceback.print_exc()
            assert(0)
    return wrapped