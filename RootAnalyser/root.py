def write_hist(hist,filename=None, normed=True):
    '''Function to write out ROOT histogram object to text file'''
    if filename is None: filename = '{}.dat'.format(hist.GetName())
    # print "writing ",filename
    with open(filename,'w') as data:
        nevents = hist.GetEntries()
        integral = hist.Integral()
        norm = integral if (normed and integral > 0.) else 1.
        data.write('# {} Entries\n'.format(nevents))
        data.write('# x\ty\tdy\n')
        for i in range(1, hist.GetNbinsX()+1):
            val = hist.GetBinContent(i)
            err = hist.GetBinError(i)
            x, y, dy = hist.GetBinCenter(i), val/norm, err/norm
            data.write('{}\t{}\t{}\n'.format( x, y, dy ))

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