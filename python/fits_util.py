from astropy.io import fits
import numpy as np
import sys
import os

def _merge_pass0(fnames,count_hdrs,time_hdr=None):
    """ Iterate through FITS files to determine the start and end times and
        to count the total number of rows needed.

    Parameters
    ----------
    fnames : the list of FITS file names
    count_hdrs : names of HDUs to count the rows of
    time_hdr : name of HDU to use to access TSTART,TSTOP,DATA, and DATE-END
        If not provided, use the first of count_hdrs.

    Returns
    -------
    dtypes : dtype for each HDR to be added
    nrows : number of rows for each of count_hdrs
    (tstart,tstop,dstart,dend) : entries for the first TSTART, etc.
    """

    nrows = [0 for hdr in count_hdrs]
    fnames = sorted(fnames)
    if time_hdr is None:
        time_hdr = count_hdrs[0]

    # Get dtype for each HDR to be added
    with fits.open(fnames[0],memmap=True) as f:
        dtypes = [f[hdr].data.dtype for hdr in count_hdrs]

    tstart = np.inf # MET
    tstop = -np.inf

    for fname in fnames:

        with fits.open(fname,memmap=True) as f:

            for ihdr,hdr in enumerate(count_hdrs):
                hdu = f[hdr]
                nrows[ihdr] += len(hdu.data)

            hdr = f[time_hdr].header
            if hdr['TSTART'] <= tstart:
                tstart = hdr['TSTART']
                dstart = hdr['DATE-OBS']
            if hdr['TSTOP'] >= tstop:
                tstop = hdr['TSTOP']
                dend = hdr['DATE-END']

    return dtypes,nrows,(tstart,tstop,dstart,dend)

def _merge_pass1(fnames,merge_hdrs,sort_fields,time_hdr=None):
    """ Combine the rows in the list of file names into a single array.

    Parameters
    ----------
    fnames : the list of FITS file names
    merge_hdrs : names of HDUs to merge
    sort_fields : list of field by which to sort each HDU; None is valid
    time_hdr : optional HDR to use to get TSTART, DATE-START, etc.
    """

    dtypes,nrows,time_info = _merge_pass0(
            fnames,merge_hdrs,time_hdr=time_hdr)

    datas = [np.zeros(nrow,dtype=dtype) for dtype,nrow in zip(dtypes,nrows)]
    counts = [0] * len(merge_hdrs)
    
    for fname in sorted(fnames):
        with fits.open(fname,memmap=True) as f:
            for ihdr,hdr in enumerate(merge_hdrs):
                n = len(f[hdr].data)
                datas[ihdr][counts[ihdr]:counts[ihdr]+n] = f[hdr].data
                counts[ihdr] += n

    for ist,st in enumerate(sort_fields):
        if st is not None:
            idx = np.argsort(datas[ist][st])
            datas[ist] = datas[ist][idx]

    return datas,time_info


def _merge_pass2(fnames,merge_hdrs,datas,time_info,outfname):
    """ Write out the combined data arrays.

    Parameters
    ----------
    fnames : the list of FITS file names
    merge_hdrs : names of HDUs to merge
    datas : a combined array for each of the HDUs to merge
    time_info : tuple of TSTART, TSTOP, DATE-OBS, DATE-END
    outfname : location of output file
    """

    tstart,tstop,dstart,dend = time_info
    with fits.open(fnames[0]) as f:
        for ihdr,hdr in enumerate(merge_hdrs):
            hdu = f[hdr]
            hdu.data = datas[ihdr]
            hdu.header['TSTART'] = tstart
            hdu.header['TSTOP'] = tstop
            hdu.header['DATE-OBS'] = dstart
            hdu.header['DATE-END'] = dend
            hdu.header['TELAPSE'] = tstop - tstart

        f.writeto(outfname,overwrite=False)
    
def combine_ft1s(fnames,outfname):
    """ Concatenate the provided FT1 files into a single file.

    Parameters
    ----------
    fnames : a list of input filenames
    outfname : destination to which to write the merged file
    """

    merge_hdrs = ['EVENTS','GTI']
    sort_fields = ['TIME','START']

    # get a merged data array for each HDU
    datas,time_info = _merge_pass1(fnames,merge_hdrs,sort_fields)

    # special processing for FT1/GTI
    _,idx = np.unique(datas[1]['START'],return_index=True)
    datas[1] = datas[1][idx]

    _merge_pass2(fnames,merge_hdrs,datas,time_info,outfname)

def combine_ft2s(fnames,outfname):
    """ Concatenate the provided FT2 files into a single file.

    Parameters
    ----------
    fnames : a list of input filenames
    outfname : destination to which to write the merged file
    """

    merge_hdrs = ['SC_DATA']
    sort_fields = ['START']

    # get a merged data array for each HDU
    datas,time_info = _merge_pass1(fnames,merge_hdrs,sort_fields)

    _merge_pass2(fnames,merge_hdrs,datas,time_info,outfname)

# NB -- this version has diverged from the one in godot, sync?
def parse_dss(fname):
    """ Return some information about the data cuts in an FT1 file using
    the DSS keywords.

    Specifically, the relevant max radius, energy, and zenith angle cuts.

    Returns
    -------
    (ra,dec,maxrad), (emin,emax), (zmax)
    """
    rcuts = None
    ecuts = None
    zmax = None
    with fits.open(fname,lazy_load_hdus=True) as f:
        h = f[1]._header
        ndss = len([x for x in h.keys() if x.startswith('DSTYP')])
        for idss in range(1,ndss+1):
            dstyp = h[f'DSTYP{idss}'].strip()
            dsval = h[f'DSVAL{idss}'].strip()
            if dsval.startswith('CIRCLE'):
                rcuts = [float(x) for x in dsval[8:-1].split(',')]
                continue
            if dstyp == 'ZENITH_ANGLE':
                zmax = float(dsval.split(':')[-1])
                continue
            if dstyp == 'ENERGY':
                ecuts = [float(x) for x in dsval.split(':')]
                continue
    return rcuts,ecuts,zmax

def get_livetime(fname):
    with fits.open(fname,lazy_load_hdus=True) as f:
        t1 = f['gti'].data['stop']
        t0 = f['gti'].data['start']
        return np.sum(t1-t0)
