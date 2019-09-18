
def locate_2D(x, y, xedges, yedges):
    # locate x, y in bins specified by edges xedges & yedges
    # returns 1st bin for underflow and last bin for overflow
    for ix in xrange(len(xedges)-1):
        if x < xedges[ix]: 
            break
        elif ((x >= xedges[ix]) and (x < xedges[ix+1])): 
            break
    for iy in xrange(len(yedges)-1):
        if x < xedges[ix]: 
            break
        elif ((y >= yedges[iy]) and (y < yedges[iy+1])): 
            break
    
    # get element of z
    zbin  =  ix + (len(xedges)-1)*iy
    # zbin  =  iy + len(yedges)*ix
    
    return zbin