compile_opt idl2

ysas_startup_m2fs
common ysas_common
common ysas_debug
spectra=ptrarr(2)

;Define initial RV, Vsini & temp_logg for the fits
rv=(-17595.795+7000)
vsini= 3.26d ;2010A&A...520A..79M
temp_logg='5500_4.55' ;close to G8V?


;Define fit regions for O49 & O50
o49fr=[350,3000]
o50fr=[2050,4095]


;Setup the pipline 
ysas_extract_ID='MARIO'
ysas_fit_id='ysastest'
ysasdb_verbose_fit=1000
ysasdb_fitting_progress_plots=[0,0,0]

;Define starting parameters and fitting configuration

o49psfheights=[0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0, 0.00]
o49psf=[3d-5, o49psfheights]

o49wave=[714.68180d-3,   3.6283054e-06,   3.8202046e-09,   -4.9249632e-12,$
	3.4142272e-15,   -1.3302967e-18,   2.7219493e-22,   -2.2759672e-26]	;From line center fits in linecenteralign.pro

o50wave=[714.49162d-3,   5.0549646e-06,   3.6017347e-10,  -1.1322256e-12, $
	1.3419806e-15,  -7.7950559e-19,   2.1506384e-22,  -2.2625270e-26]
	
o50psfheights=[0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0, 0.00]
o50psf=[3d-5, o50psfheights]

o50wave=[714.68180d-3,   3.6283054e-06,   3.8202046e-09,   -4.9249632e-12,$
	3.4142272e-15,   -1.3302967e-18,   2.7219493e-22,   -2.2759672e-26]

;Best output so far, solution is not right.
;Wave:      0.71467391   3.6286251e-06   3.8240148e-09  -4.9260691e-12   3.4141454e-15  -1.3304221e-18   2.7215776e-22  -2.2728853e-26
;PSF:   2.5731019e-05    -0.090318988     0.031939198    -0.024293696   0.00016124602     0.024945085    0.0083763979     0.010948209   -0.0046798999
;Airmass:      0.64255258
;Veil:      0.42095779
;Vsini:       4.2038588
;RV:      -9768.2394
;Norm:     0.076972419    -0.010488462



o49init=[o49wave,  o49psf, .63, .4, vsini, .6, rv/299792458d, 0, .01]
o50init=[o50wave,  o50psf, .63, .4, vsini, .6, rv/299792458d, 0, .01]

;Define tolerance
tol=1d-11

;Define wavelength & PSF scale
wavescale=[5d-3, 1d-7, 1d-11, 1d-13,1d-16,1d-19,1d-23,1d-27]/10 
psfsigmascale=1d-5

psfhieghtscale=replicate(.01,8)


waveub=[.72,1,1,1,1,1,1,1]
wavelb=[.71,-1,-1,-1,-1,-1,-1,-1]

psfub=replicate(2,n_elements(ysas_psf_param_ndxs))
psflb=replicate(-2,n_elements(ysas_psf_param_ndxs))

;Define different fitting configurations


psfscale=[psfsigmascale,psfhieghtscale]

fixed=intarr(ysas_MAX_NUM_FIT_PARAMS)
fixed[*]=1
fixed[ysas_wavelength_param_ndxs]=0
waveconfig={$
	fixed:fixed $
	,ub:waveub $
	,lb:wavelb $
	,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
	,scale: wavescale $
	,mode:'Final' $
	,fftenable:0b $
	,fromexisting:1b $
	,maxreps:10 $
	,tol:tol $
	,vsini:vsini}

fixed[*]=1
fixed[ysas_psf_param_ndxs]=0
psfconfig={fixed:fixed $
	,ub:psfub $
	,lb:psflb $
	,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
	,scale: psfscale $
	,mode:'Final' $
	,fftenable:0b $
	,fromexisting:1b $
	,maxreps:10 $
	,tol:tol $
	,vsini:vsini}
	
fixed[*]=1
fixed[ysas_wavelength_param_ndxs]=0
fixed[ysas_psf_param_ndxs]=0
wavepsfconfig={fixed:fixed $
	,ub:[waveub,   psfub] $
	,lb:[wavelb,   psflb] $
	,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
	,scale:[wavescale, psfscale] $
	,mode:'Final' $
	,fftenable:0b $
	,fromexisting:1b $
	,maxreps:10 $
	,tol:tol $
	,vsini:vsini}

fixed*=0
fixed[ysas_rotational_broadening_param_ndxs[1]]=1
allconfig={fixed:fixed $
	,ub:[waveub, psfub,   2, 2,   25,  .0005d,   0.5, 40000] $
	,lb:[wavelb, psflb,  .1,.1,   .1, -.0005d,   -.5,-40000] $
	,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
	,scale:[wavescale, psfscale, .1, .1, 1, 25000/299792458d, 0.01, 0.001] $
	,mode:'Final' $
	,fftenable:0b $
	,fromexisting:1b $
	,maxreps:3 $
	,tol:tol $
	,vsini:vsini}



; Load orders 49 & 50 
template={datastart:0, delimiter:32b,fieldcount:2,fieldgroups:[0,1],$
    fieldlocations:[0,4],fieldnames:['pixel','flux'],fieldtypes:[2,5],$
    missingvalue:!values.d_nan,version:1.,commentsymbol:'#'}

order=50
ysas_fit_region=o50fr
data1=read_ascii('/Users/one/Downloads/ap1.txt', template=template)
o50=buildEmpiricSpectrum(data1.flux,  order - 1, NORMALIZE=3, $
    EXTRACTID=ysas_extract_ID, instrument='M2FS_hires')
o50.synthspec=temp_logg
o50.fitparams=o50init

order=49
ysas_fit_region=o49fr
data2=read_ascii('/Users/one/Downloads/ap2.txt', template=template) ;49
o49=buildEmpiricSpectrum(data2.flux,  order-1, NORMALIZE=3, $
    EXTRACTID=ysas_extract_ID, instrument='M2FS_hires')
o49.synthspec=temp_logg
o49.WAVELENGTHDOMAIN[0]=.712
o49.fitparams=o49init


;Configure to fit O49
spectra[0]=ptr_new(o50)
ysas_fit_region=o50fr


config=allconfig

;Do the fit
stop
fitSpectra, spectra, '', temp_logg, airmass, results, $
    CONFIG=config, STATUS=STATUS, CHISQUARE_LIMIT=1000d, $
    MAX_FITSPECTRA_CYCLES=10
o50fit=*spectra[0]
stop

end

;manual_calibrate_alignment,o49,o49fit.fitparams,out
;manual_calibrate_alignment,o49,[o49wave,  3d-5, .677, .3, vsini ,.6, 13302/299792458d, 0, .01],out

;
;manual_calibrate_alignment,o50,[o50wave,  psfparams, .677, .3, vsini ,.6, 13302/299792458d, 0, .01],out


;[o50wave,  3d-5, .7, .3, vsini, .6, 13302/299792458d, 0, .01]

;first try o49:
;0.71449728   5.1479318e-06  -1.0776708e-10   6.6942156e-05      0.67692639      0.30263929       3.6783614      0.60000002   4.4342719e-05    0.0027254524    -0.015491522