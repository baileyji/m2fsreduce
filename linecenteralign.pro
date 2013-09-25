;o49dat=[[457.14069d, 0.71679313],$	;Results of run for r0090.fits
;	[643.35773, 0.71774360],$
;	[784.12265, 0.71845511],$
;	[819.41398, 0.71863334],$
;	[922.82605, 0.71915326],$
;	[964.64642, 0.71936858],$
;	[1118.9285, 0.72012327],$
;	[1180.4915, 0.72043332],$
;	[1769.6796, 0.72329326],$
;	[1808.0729, 0.72347507],$
;	[1934.3088, 0.72407463],$
;	[2448.2838, 0.72645955],$
;	[2633.9824, 0.72730036],$
;	[2731.8609, 0.72774024]]
;print,poly_fit(o49dat[0,*],o49dat[1,*],7)

function simpleWavelengthSolnGuess, spec
	common ysas_common
	
	;a*x+c
	a=(spec.wavelengthdomain[1]-spec.wavelengthdomain[0])/ $
	  (ysas_fit_region[1]-ysas_fit_region[0])
	c=spec.wavelengthdomain[0]-a*ysas_fit_region[0]
	return, [c,a,replicate(0,n_elements(ysas_wavelength_param_ndxs)-2)]
end


pro line_locations, order, line_wave, wid_wave, wid_pix

	line_wave_o49=[0.71679313d, 0.71774360, 0.71845511, 0.71863334, 0.71915326, $
		0.71936858, 0.72012327, 0.72043332, 0.72329326, 0.72347507, 0.72407463, $
		0.72645955, 0.72730036, 0.72774024]
	wid_pix_o49=replicate(15,n_elements(line_wave_o49))
	wid_wave_o49=replicate(.00010,n_elements(line_wave_o49))
	
	line_wave_o49=[0.71679313d, 0.71774360, 0.71845511, 0.71863334, 0.71915326, $
		0.71936858, 0.72012327, 0.72043332, 0.72329326, 0.72347507, 0.72407463, $
		0.72645955, 0.72730036, 0.72774024]
	wid_pix_o49=replicate(15,n_elements(line_wave_o49))
	wid_wave_o49=replicate(.00010,n_elements(line_wave_o49))
	
	if order eq 49 then begin
		line_wave=line_wave_o49
		wid_pix=wid_pix_o49
		wid_wave=wid_wave_o49
	endif else begin
		line_wave=line_wave_o50
		wid_pix=wid_pix_o50
		wid_wave=wid_wave_o50
	endelse

end

pro primelambdasoln, spec, rv, vsini, first_guess=first_guess


	common ysas_common
	common ysas_debug

	if ~keyword_set(first_guess) then begin
		;Phase One: Get a halfway decent solution
		;Define a fitting configuration
		psfguess=[3d-5, replicate(0, n_elements(ysas_psf_param_ndxs)-1)]
		init=[simpleWavelengthSolnGuess(spec), $
			  psfguess,  .63, .4, vsini, .6, rv/3d8, 0, .01]
		wavescale=[5d-3, 1d-7, 1d-11, 1d-13,1d-16,1d-19,1d-23,1d-27]
		psfsigmascale=1d-5
		psfhieghtscale=replicate(.01,8)
		psfscale=[psfsigmascale,psfhieghtscale]
		waveub=[spec.wavelengthdomain[0]*1.1,1,1,1,1,1,1,1]
		wavelb=[spec.wavelengthdomain[0]*.9,-1,-1,-1,-1,-1,-1,-1]
		psfub=replicate(2,n_elements(ysas_psf_param_ndxs))
		psflb=replicate(-2,n_elements(ysas_psf_param_ndxs))
		fixed=intarr(ysas_MAX_NUM_FIT_PARAMS)
		fixed[ysas_rotational_broadening_param_ndxs[1]]=1
		config={fixed:fixed $
			,ub:[waveub, psfub,   2, 2,   25,  .0005d,   0.5, 40000] $
			,lb:[wavelb, psflb,  .1,.1,   .1, -.0005d,   -.5,-40000] $
			,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
			,scale:[wavescale, psfscale, .1, .1, 1, 25000/3d8, 0.01, 0.001] $
			,mode:'Final' $
			,fftenable:0b $
			,fromexisting:1b $
			,maxreps:3 $
			,tol:1d-8 $
			,vsini:vsini}
	
		;Fit the spectrum
		ysas_final_fit_plot=1b
		spectra=[ptr_new(spec)]
		(*spectra[0]).fitparams=init
		fitSpectra, spectra, '', spec.synthspec, airmass, results, $
			CONFIG=config, STATUS=STATUS, CHISQUARE_LIMIT=1000d, $
			MAX_FITSPECTRA_CYCLES=3
		
		;Print the resulting fit
		printparams, (*spectra[0]).fitparams
		
		;Make sure the solution is good enough
		;Accept if the user clicks in the top half, reject in the lower half
		CURSOR,x0,y0,/down
		if  y0 lt .5 then begin
			stop
			ptr_free,spectra[0]
			return
		endif else begin
			spec.fitparams=(*spectra[0]).fitparams
			ptr_free,spectra[0]
		endelse
	endif else begin
		spec.fitparams=first_guess
	endelse

	;Get the line information	
	line_locations, order, line_wave, wid_wave, wid_pix
	
	;Compute the model
	telluricSpectrum=getTelluricSpectrum(spec)
	synthSpectrum=getSynthSpectra(spec.synthspec, spec)
	pix=indgen(spec.nPixels)
	model = computeSpectralModel(pix, synthSpectrum, $
		telluricSpectrum, params, spectrum = spec, /USE_ZSPACE, $
		pixel_wavelengths = wave, /RETURN_INTENSITIES, telluric_out=tspec)
	
	;Data Array
	obs_cent=dblarr(2, n_elements(line_wave))

	;Pixel space plot for clicking
	!p.multi=[0,2,1]
	getwindow,1,xsize=1000,ysize=500,/erase
	col=2
	
	;Fit each of the lines
	for i=0,n_elements(line_wave)-1 do begin
		
		
		;Find the nearest large minimum in the spectrum and get it's pixel
		; position, pix_center
		junk=min(wave - line_wave[i], pix_center)
		 
		;Extract a width about that pixel position & fit with a gaussian
		x_pix=indgen(wid_pix[i])-wid_pix[i]/2+pix_center
		y_pix=(*spec.intensities)[x_pix]
		pixy=mpfitpeak(x_pix,y_pix,fit_pix,/negative)
		
		;Extract a width of the telluric spectrum about the wavelength &
		; fit with a guassian
		max_wave=line_wave[i]+wid_wave[i]/2d
		min_wave=line_wave[i]-wid_wave[i]/2d
		wave_mask=where((*tspec.wavelengths) gt min_wave and $
						(*tspec.wavelengths) lt max_wave)
		x_wave=(*tspec.wavelengths)[wave_mask]
		y_wave=(*tspec.intensities)[wave_mask]
		wavey=mpfitpeak(x_wave,y_wave,fit_wave,/negative)

		;Plot the two fits along with a region of the spectrum for acceptance
		xrange=[x_pix[0]-75,x_pix[-1]+75]
		plot,pix[mask], (*spec.intensities)[mask], $
			YTITLE='Data (white) Tellurics (red)', title='Line Plot', $
			xr=xrange,/xs

		oplot,x_pix,y_pix,color=col
		col++
		oplot,x_pix,pixy,color=col
		col--

		plot,(*tspec.wavelengths),(*tspec.intensities), $
			title='Telluric Data',/xs, xr=[x_wave[0]-.0005,x_wave[-1]+.0005]
		oplot,x_wave,y_wave,color=col
		col++
		oplot,x_wave,wavey,color=col

		col++

		;Accept if the user clicks in the top half, reject in the lower half
		CURSOR,x0,y0,/down
		if  y0 gt .5 then begin
			obs_cent[*,i]=[fit_pix[1],fit_wave[1]]
			print, obs_cent[*,i]
		endif


	endfor
	!p.multi=0
	
	;Compute a wavelength solution
	soln=poly_fit(obs_cent[*,0],obs_cent[*,1],$
		n_elements(ysas_wavelength_param_ndxs)-1)
	print, soln
	
	destroySpec,[telluricSpectrum, synthSpectrum, tspec]
	
	spec.fitparams[ysas_wavelength_param_ndxs]=soln
	
end