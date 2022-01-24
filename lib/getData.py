import os
import fire
import shutil

import numpy as np
import astroquery.mast as mast
import astropy.io.fits as pyfits

def getData (pi, target, pid, sciFil, imgFil, fileSearch='ima.fits', local='.', sortPath='.'):
    """
    Description
    ----------
    Downloads data from MAST for a given pi, target, and proposal id, and sorts into directories for WFC3.

    Parameters
    ----------
    pi: string
        Name of pi, wild cards are allowed
    target: string
        Name of target star
    pid: string
        Proposal id number.
    sciFil: string
        Filter used for science observations
    imgFil: string
        Filter used for direct images (0th order spectra)
    fileSearch: string (optional)
        Name subsearch for data files. Defaults to ima.fits
    local: string (optional)
        Directory to store mastDownload hierarchy of directories.
    sortPath: string (optional)
        Directory to store sorted files for WFC3

    Returns
    -------
    None, but makes lots of files.

    Warning
    -------
    This was intended to download scan files, and so hard coded to sort for SCAN_TYP='C'.
    Will duplicate files so you don't redownload data. Watch your hard disk space. 

    """
    local = os.path.expanduser(local)
    sortPath = os.path.expanduser(sortPath)

    # find all observations associated with a target in a given pi's proposal
    obs_table = mast.Observations.query_criteria(proposal_pi=pi, target_name=target, proposal_id=int(pid))

    # get a list of associated data products
    data_products_by_obs = mast.Observations.get_product_list(obs_table)

    # search for only the ima fits files
    indices = [row for row in np.arange(0, len(data_products_by_obs)) if fileSearch in data_products_by_obs['productFilename'][row]]

    # download only these data products and don't redownload if you already have them
    print(f'Downloading {len(indices)} Files:')
    manifest = mast.Observations.download_products(data_products_by_obs[indices], extension="fits", cache=True, download_dir=local)

    # sort through headers looking for image files, cals, and science files
    cals = []
    images = []
    science = []
    for file in manifest['Local Path']:
        with pyfits.open(file) as hdul:
            filt = hdul[0].header['FILTER']
            scan = hdul[0].header['SCAN_TYP']
        if (filt == sciFil) and (scan == 'C'):
            science.append(file)
        elif (filt == imgFil) and (scan == 'N'):
            images.append(file)
            cals.append(file)
        else:
            cals.append(file)

    masterPath = f"{sortPath}/{target}-{pid}"
    print(f"Sorting files into {masterPath}")

    # move cal files, checks if files exist first, does nothing if they do

    if not os.path.exists(f'{masterPath}/cals'):
        os.makedirs(f'{masterPath}/cals')
    for file in cals:
        if not os.path.exists(f"{masterPath}/cals/{file.split('/')[-1]}"):
            shutil.copy2(file, f"{masterPath}/cals/{file.split('/')[-1]}")

    # move science files, checks if files exist first, does nothing if they do
    if not os.path.exists(f'{masterPath}/sci'):
        os.makedirs(f'{masterPath}/sci')
    for file in science:
        if not os.path.exists(f"{masterPath}/sci/{file.split('/')[-1]}"):
            shutil.copy2(file, f"{masterPath}/sci/{file.split('/')[-1]}")

    # make directImg list, checks if file exists first, does nothing if it does
    if not os.path.exists('directImg.txt'):
        directImg = [file.split('/')[-1] for file in images]
        np.savetxt('directImg.txt', directImg, fmt='%s')

    return f"{sortPath}{target}-{pid}"

# fire functionality lets you run this from either a jupyter notebook, another .py file, or the command line
# python getData.py Desert* HAT-P-2 16194 G141 F126N --fileSearch=ima.fits --local=~/supportData --sortPath=~/supportData
if __name__ == '__main__':
    fire.Fire(getData)
