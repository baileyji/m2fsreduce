#from numpy import *
import numpy as np
from astropy.io import fits


#def mergeimage(frameno,side=None):
    frameno=84
    side='b'
    if side == None:
        mergeimage(frameno,side='r')
        mergeimage(frameno,side='b')
        return
    basename=side+'{frameno:04}'
    f=basename+'c{quad}.fits'

    quadrant1=fits.open(f.format(frameno=frameno, quad=1))
    quadrant2=fits.open(f.format(frameno=frameno, quad=2))
    quadrant3=fits.open(f.format(frameno=frameno, quad=3))
    quadrant4=fits.open(f.format(frameno=frameno, quad=4))


    #Grab the bias and crop region from the first quadrant
    biassec=quadrant1[0].header['BIASSEC']
    bias=[int(s) for s in re.findall(r'\d+',biassec)]

    trimsec=quadrant1[0].header['TRIMSEC']
    crop=[int(s) for s in re.findall(r'\d+',trimsec)]

    #Create an output array
    out=np.zeros((crop[1]*2,crop[3]*2),dtype=np.int16)

    #Define where the quadrants show go
    quadLoc=[(0, 2048, 0, 2056),
             (2048, 4096, 0,2056),
             (2048, 4096, 2056, 4112),
             (0, 2048, 2056, 4112)]

    #Create a list of the quadrant's data to loop through
    quadrantData=[quadrant1[0].data, quadrant2[0].data,
                 quadrant2[0].data, quadrant4[0].data]

    for i, qdata in enumerate(quadrantData):

        #Compute mean for overscan region rowwise
        biaslevels=np.mean(qdata[crop[0]-1:crop[1],bias[2]-1:bias[3]],axis=1)

        #Crop image & subtract bias levels row wise
        qdata=(qdata[crop[0]-1:crop[1],crop[2]-1:crop[3]] -
                          biaslevels[:,np.newaxis])


        out[quadLoc[i][0]:quadLoc[i][1],quadLoc[i][2]:quadLoc[i][3]]=qdata


    #Write out the merged file
    hdu = fits.PrimaryHDU(out)
    hdu.header=quadrant1[0].header
    hdu.header.pop('TRIMSEC')
    hdu.header.pop('BIASSEC')
    hdu.header.pop('DATASEC')
    hdu.header['FILENAME']=basename
    hdu.writeto(basename+'.fits')



def makesuperflat(files, name, superbias=None):
    #load first flat to get info
    with fits.open(files[0]) as flat
        pass
    
    #Number of flats
    nflat=len(files)

    #Make the datacube
    cube=np.zeros((nrow,ncol, nflats), dtype=np.float32)

    for f in enumerate(files):
        with fits.open(f) as flat:
        
            #Remove CRs & put in cube
            
            cube[:,:,i]=sigmaClip(flat[0].data)

    #Merge flat

    #Subtract the superbias from the superflat
    try:
        with fits.open(superbias) as bias:
               flat.data-=bias[0].data
    except IOError:
        pass
        
    #Write out the superflat
    hdu.header['COMMENT']=''.join([)
    hdu.writeto(name)


def makesuperbias(files, name):
    #load first bias to get info
    with fits.open(files[0]) as bias
        pass
    
    #Number of bias frames
    nbias=len(files)
    
    #Make the datacube
    cube=np.zeros((nrow,ncol, nbias), dtype=np.int16)
    
    for f in enumerate(files):
        with fits.open(f) as bias:
            cube[:,:,i]=bias[0].data
    
    #Merge bias

    hdu.writeto(name)
