# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 14:26:04 2014

@author: williams savero torres

Modifications by Piotr:
- adding image histogram equalization for better noise ration detection
- resolving system paths problem -> using op.path.join
- sorting input files by names
- adjusting (boosting) plot properties

Updated by Vinh Dang on 18 - Mar - 2015:
- fix the bug when the number of maxmas less than 4. Now, the number of maxmas will be sure 4
- make the curve smooth for better detection
- support to process 1 single file
"""
from __future__ import division
import numpy as np
from numpy import*
import scipy.ndimage as ndimage
#import matplotlib.pyplot as plt
import pylab
#import scipy.misc
import scipy
from pylab import *
from skimage.segmentation import random_walker
from scipy import ndimage
#for image adjustments:
from skimage.exposure import equalize_hist

def BATCH(kerr):
    #reading image
    img = ndimage.imread(kerr)
    #equalizing histogram - added by Piotr
    img = equalize_hist(img)

    # Note the 0 sigma for the last axis, we don't wan't to blurr the color planes together!
    img = ndimage.gaussian_filter(img, sigma=(2, 2), mode='nearest', order=0)

    # Make a line with "num" points...
    # Modify x0, x1, y0, y1 if needed to make sure the line cut two circles
    x0, y0 = 1080, 438 
    x1, y1 = 1080, 1388  
    num = 10000
    x, y = np.linspace(x0, x1, num), np.linspace(y0, y1, num)

    # Extract the values along the line, using cubic interpolation
    zi = scipy.ndimage.map_coordinates(np.transpose(img), np.vstack((x,y)), order=5)
    #print zi[0]

    # Store the original value of zi
    z0 = zi
    # Modify if needed to adjust the smooth function
    smooth_range = 50
    d1 = []
    for i in range (len (zi)):
        sum = 0.0
        count = 0
        for j in range (i - smooth_range, i + smooth_range + 1):
            if (j >= 0 and j < len (zi)):
                sum += zi [j]
                count += 1
        d1.append (sum / count)
    zi = d1

    # Ma = x[zi == max(zi)]
    d=diff(zi, n=1, axis=-1)
    mean_d=mean(abs(d))
    #print d

#!--------------peaks definition--------!#
    def peakdet(v, delta, x = None):

        maxtab = []
        mintab = []

        if x is None:
            x = arange(len(v))

        v = asarray(v)

        if len(v) != len(x):
            sys.exit('Input vectors v and x must have same length')

        if not isscalar(delta):
            sys.exit('Input argument delta must be a scalar')

        if delta <= 0:
            sys.exit('Input argument delta must be positive')

        mn, mx = Inf, -Inf
        mnpos, mxpos = NaN, NaN

        lookformax = True

        for i in arange(len(v)):
            this = v[i]
            if this > mx:
                mx = this
                mxpos = x[i]
            if this < mn:
                mn = this
                mnpos = x[i]

            if lookformax:
                if this < mx-delta:
                    maxtab.append((mxpos, mx))
                    mn = this
                    mnpos = x[i]
                    lookformax = False
            else:
                if this > mn+delta:
                    mintab.append((mnpos, mn))
                    mx = this
                    mxpos = x[i]
                    lookformax = True

        return maxtab, mintab

    #!---------------------------------------------------------------
    #! setting up threshold for detection
    treshold=3*mean_d
    det=peakdet(abs(d), treshold, x = None)
    print "detection =",det
    for i in range(len(det[0])):
        print "%i peak position = %.3f"%(i,(y[det[0][i][0]]))
    #print "number of maxmas found =",len(det[0])
    #print "number of minimas found =",len(det[1])

    #! testing number of maximas
    #!---------------------------------------------------------------
    while len(det[0])!= 4:
        if (len(det[0]) > 4):
            print "not sufficient threshold was used..."
            print "using 2% higher threshold.."
            treshold=treshold + 0.02*treshold
            det=peakdet(abs(d), treshold, x = None)
            print "new number of maxmas found with higher threshold =",len(det[0])
        else:
            print "not enough threshold was used..."
            print "using 5% lower threshold.."
            treshold=treshold - 0.05*treshold
            det=peakdet(abs(d), treshold, x = None)
            print "new number of maxmas found with lower threshold =",len(det[0])

    #=======================================Ploting
    figure(figsize=(16,8))

    subplot(311,xticklabels=[],yticklabels=[])
    pylab.gray()
    imshow(img,aspect='equal')
    plot(x, y, 'r-', linewidth=1)

    subplot(312,xticklabels=[])
    plot(y,zi,'b-',lw=3)
    plot (y, z0, 'y-', lw = 3)
    #if len(det[0])==4:
    d1=((y[det[0][1][0]]-y[det[0][0][0]]))
    d2=((y[det[0][3][0]]-y[det[0][2][0]]))

    figtext(0.7,0.95,"$d_{up}[\mu m]=%.6f$"%float(d1),size=20)
    figtext(0.7,0.85,"$d_{down}[\mu m]=%.6f$"%float(d2),size=20)
    dt1_list.append(d1)
    dt2_list.append(d2)
    #else:
     #   pass
    summary.write( '%s \t %i \t %.3f \t %.3f \n'%(str(os.path.basename(kerr[:-4])),count,d1,d2))
    summary.flush()

    for i in range(len(det[0])):
        axvline(x=y[det[0][i][0]],color='r',lw=3)

    subplot(313,yticklabels=[])
    plot(y[:-1],d,'g-',lw=3)
    axhline(y=(treshold),color='r',label='$treshold$',lw=3)
    axhline(y=(-treshold),color='r',label='$treshold$',lw=3)
    fill_between(y[:-1], 0, abs(d) ,alpha=0.3,facecolor='g', interpolate=True)
    tight_layout()
    savefig(os.path.join(wdir,'analyses','%s.png'%(str(os.path.basename(kerr[:-4])))))
    clf()
    close()
    #figure.close()

#!_________________________PROGRAMS LOOP___________________
#!---------------------------------------------------------------
print "starting program.."
import os
import glob
import sys
#! defining working directory
wdir=str(os.getcwd())

#! checking existance of the analyses directory
if not os.path.exists(os.path.join(wdir,'analyses')):
    os.makedirs(os.path.join(wdir,'analyses'))
    print "analyses directory created"

#! creating summary file with header
summary = open(os.path.join(wdir,"analyses","SUMMARY.data"), "w")
summary.write( '#filename \t count \t d1[um] \t d2[um]\n')
summary.flush()

#! creating empty  d times list
dt1_list=[]
dt2_list=[]
#! iterating on all files in the directory
if (len(sys.argv) == 1):
    count=1
    for infile in sorted(glob.glob(os.path.join(wdir,'*.png') )):
        print "\n current file is: " + infile
        BATCH(infile)
        count+=1

    print "no more files have been found"
elif (len(sys.argv) == 2):  #process a particular file
    BATCH (sys.argv[1])

summary.close()
print "\n ALL DONE \n"
