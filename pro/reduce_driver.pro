compile_opt idl2
ysas_startup_m2fs
common ysas_common

imfile='Desktop/M2FS Commissioning Run/image/b_superflat_84-86.fits'
fits_read,imfile,im,head

str={filename:'r0104.fits', header:head,$
     flux:float(im),err:sqrt(float(im>1)),$
     mask:long(im*0)}
findtrace,str,tracestr,dcr=4,halfwid=8,npts=40

getwindow,0,xs=1500,ys=1000,/show
loadct,0
display,im,min=0,max=500

colors
cols=rnggen(0,4095)
for i=0, n_elements(tracestr)-1 do begin
    pix=cols[tracestr[i].firstx:tracestr[i].lastx]
    tr=poly(pix,tracestr[i].fitcoef)
    oplot,pix,tr,col=3
    if tracestr[i].fixed eq 1 then xyouts,pix[0],tr[0],'Fixed'
endfor


touse=where(tracestr.fitcoef[0] lt 2500 and tracestr.fitcoef[0] gt 2000)
orders=tracestr[touse].fitcoef
colrange=transpose([[tracestr[touse].firstx],[tracestr[touse].lastx]])

SXADDPAR,head, 'E_READN',ysas_readNoise
SXADDPAR,head, 'E_BACKG',0.0
SXADDPAR,head, 'E_GAIN',ysas_GAIN

stop
;Determine extraction windows above and below order maximum.
  getxwd,ordim,orders,def_xwd,def_sxwd,colrange=col_range ;get extraction width


  hamflat, flat, fhead, orders, blzcoef, colrange=col_range, FXWD=def_xwd $
         , SF_SMOOTH=20., SP_SMOOTH=20., OSAMPLE=12, swath_width=400 $
         , MASK=mask, /PLOT


stop
hamflat, im, head, orders, blzcof, COLRANGE=colrange $
           , OSAMPLE=osample, FXWD=fxwd, SF_SMOOTH=sf_smooth $
           , SP_SMOOTH=sp_smooth, SWATH_WIDTH=swath_width, MASK=mask $
           , POLY=pol, PLOT=1 $
           , NOSCATTER=noscatter, THRESHOLD=100 $
           , SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , UNCERTAINTY=bunc
           
stop

oimfile='Desktop/M2FS Commissioning Run/image/b_superflatnorm_84-86.fits'
fits_write, oimfile, im, head

nflat=im
fits_read,'Desktop/M2FS Commissioning Run/image/b_superbias_93-126.fits',sbias,head

fileroot='Desktop/M2FS Commissioning Run/image/b'
im=sbias*0d
for i=77,83 do im+=readfits(fileroot+string(i,format='(i04)')+'.fits')-sbias
im/=nflat
imsave=im
stop

hamspec, im ,head, orders, dxwd, dsxwd, spec $
           , SIG=sunc, COLRANGE=colrange, SF_SMOOTH=10.0 $
           , OSAMPLE=10, SWATH_WIDTH=swath_width $
           , PLOT=1, order_range=[5,6],swath_boundaries=swath_boundaries $
           , FILENAME=filename $
           , SWATH_BOUNDARIES=swath_boundaries $
           , WING_SMOOTH_FACTOR=wing_smooth_factor $
           , SLIT_TILT_FILE=slit_tilt_file

    
end