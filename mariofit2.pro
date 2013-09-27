function reviewspecs, spectra

	nspec=n_elements(spectra)
	DOFIT=bytarr(nspec)

	;Goal is to come up with an initial estimate of a wavelength solution for 
	; alll the spectra, the estimates should only be a function of chip
	; location so find an initial estimate for the coadded data apertures and 
	; apply that to the other apertures
	
	;pull indices for coadded r & coadded b spectra
	coadd=[49,48,47,46,45,44, 97,96,95,94,93,92]-2
	
	;Pull fiber numbers for all the spectra
	fiberno=intarr(nspec)
	for i=0,nspec-1 do fiberno[i]=(*spectra[i]).fiberno

	;Go through the R spectra and get initial guesses
	for i=0, n_elements(coadd)-1 do begin
		plotspectrumfit,*spectra[coadd[i]],xr=ysas_fit_region,$
			star=file_basename((*spectra[coadd[i]]).specfile)
		choice=clickquad()
		case choice of
			0: dofit[coadd[i]]=0
			1: begin
				manual_calibrate_alignment,*spectra[coadd[i]],$
					(*spectra[coadd[i]]).fitparams,out, order=2
				(*spectra[coadd[i]]).fitparams=out
				i--
			end
			else: dofit[coadd[i]]=1
		endcase
	endfor
    
	;Copy the initial guesses into the non coadded frames
	for i=0, n_elements(coadd)-1 do begin
		if dofit[coadd[i]] then begin
			match=where(fiberno eq (*spectra[coadd[i]]).fiberno, nmatch)
			for j=0,nmatch-1 do begin

				(*spectra[match[j]]).fitparams=(*spectra[coadd[i]]).fitparams
	
				plotspectrumfit,*spectra[match[j]],xr=ysas_fit_region,$
					star=file_basename((*spectra[match[j]]).specfile)
					
				choice=clickquad()
				
				case choice of
					0: dofit[match[j]]=0
					1: begin
						manual_calibrate_alignment,*spectra[match[j]],$
							(*spectra[match[j]]).fitparams,out, order=2
						(*spectra[match[j]]).fitparams=out
					end
					else: dofit[match[j]]=1
				endcase
				
			endfor
		endif
	endfor

	return,dofit
	
end



compile_opt idl2

ysas_startup_m2fs
common ysas_common
common ysas_debug


;Define initial RV, Vsini & temp_logg for the fits
rv=(-17595.795+7000)
vsini= 3.26d ;2010A&A...520A..79M
temp_logg='5500_4.55' ;close to G8V?


;Setup the pipline
ysas_extract_ID='MARIO'
ysas_fit_id='ysastest'
ysas_fit_region=[350,3000]
ysasdb_verbose_fit=1000
ysasdb_fitting_progress_plots=[0,0,0]


;Define starting parameters and fitting configuration

o49psfheights=[0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0, 0.00]
o49psf=[3d-5, o49psfheights]

o49wave=[714.68180d-3,   3.6283054e-06,   3.8202046e-09,   -4.9249632e-12,$
	3.4142272e-15,   -1.3302967e-18,   2.7219493e-22,   -2.2759672e-26]	;From line center fits in linecenteralign.pro

o49waveprimer=[714.16610d-3,   5.2577146e-06,  -1.4535988e-10,replicate(0d,5)]

o49init=[o49waveprimer,  o49psf, .63, .4, vsini, .6, rv/3d8, 0, .01]

;Best output so far, solution is not right.
;Wave:      0.71467391   3.6286251e-06   3.8240148e-09  -4.9260691e-12   3.4141454e-15  -1.3304221e-18   2.7215776e-22  -2.2728853e-26
;PSF:   2.5731019e-05    -0.090318988     0.031939198    -0.024293696   0.00016124602     0.024945085    0.0083763979     0.010948209   -0.0046798999
;Airmass:      0.64255258
;Veil:      0.42095779
;Vsini:       4.2038588
;RV:      -9768.2394
;Norm:     0.076972419    -0.010488462


;Grab the spectra
snr_limit=12
spectra=loadAug13Spectra()
spectra=spectra[where(spectra ne ptr_new(),nspec)]


for i=0, nspec-1 do begin

		(*spectra[i]).fitparams=o49init
		;Do this so I can call plotspectrumfit
		(*spectra[i]).fitregion=strjoin(strtrim(ysas_fit_region,1),',')
		(*spectra[i]).fitparams[ysas_rv_param_ndx]=$
			(*spectra[i]).ysasheader.baryv/3d8

endfor


dofit=reviewspecs(spectra)

;Filter based on snr
snr=dblarr(nspec)
for i=0,nspec-1 do snr[i]=computes2n(*spectra[i],ysas_fit_region)
snrspectra=spectra[where(snr gt snr_limit,nspec)]
snrdofit=dofit[where(snr gt snr_limit,nspec)]



;Phase One
stop
fitSpectra, spectra, '', temp_logg, airmass, results, $
	CONFIG=getM2FSfitconfig('init', tol=1d-6), STATUS=STATUS, $
	CHISQUARE_LIMIT=1000d, MAX_FITSPECTRA_CYCLES=5, EXCLUDE=~DOFIT, $
	CHISQUARE_THRESH=1d

stop

fitSpectra, spectra, '', temp_logg, airmass, results, $
	CONFIG=getM2FSfitconfig('loose', tol=1d-6), STATUS=STATUS, CHISQUARE_LIMIT=1000d, $
	MAX_FITSPECTRA_CYCLES=5, EXCLUDE=~DOFIT, CHISQUARE_THRESH=1

stop
plotspectrumfit
end