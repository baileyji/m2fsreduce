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

pro picklineinterface, spec, params, info,fit_region,datapt,datai,WAVE=wave
    common ysas_common
    common ysas_debug
    telluricSpectrum=getTelluricSpectrum(spec)      ;+2 pointers
    synthSpectrum=getSynthSpectra(spec.synthspec, spec) ;+2 pointers
    
    pix=indgen(spec.nPixels)
    mask=indgen(ysas_fit_region[1]-ysas_fit_region[0]+1)+ysas_fit_region[0]
    
    model = computeSpectralModel(pix, synthSpectrum, telluricSpectrum, $                                                        +6 pointers
        params, spectrum = spec, bad_model = bad_model, /USE_ZSPACE, $
        pixel_wavelengths = wave,/RETURN_INTENSITIES,telluric_out=tspec,$
        synth_out=synspec)
    
    tonlyparams=params
    tonlyparams[ysas_veiling_param_ndx]=0
    tonly = computeSpectralModel(pix, synthSpectrum, telluricSpectrum, $                                                        +6 pointers
        tonlyparams, bad_model = bad_model, /USE_ZSPACE, $
        /RETURN_INTENSITIES)
    sonlyparams=params
    sonlyparams[ysas_airmass_param_ndx]=0
    sonly = computeSpectralModel(pix, synthSpectrum, telluricSpectrum, $                                                        +6 pointers
            sonlyparams, bad_model = bad_model, /USE_ZSPACE, $
            /RETURN_INTENSITIES)

    residuals=(*spec.intensities) - model
    
 
    ;Pixel space plot for clicking
    getwindow,0,xsize=2300,ysize=700,/erase
    wait,.1
    padding=60
    xrange=[pix[mask[0]]-padding,pix[mask[-1]]+padding]
    pltn=0
    plot,mask,sonly[mask],ytitle='Synth',title='Reference Plots',xr=xrange, $
        /xs,_EXTRA=gang_plot_pos(4,1,pltn++, OFFSET=[.03,.1],size=[.96,.85]), $
        yr=[0,1.2],/ys
    plot,mask,tonly[mask],ytitle='Tellu',xr=xrange, yr=[0,1.2],/ys,$
        /xs,_EXTRA=gang_plot_pos(4,1,pltn++, OFFSET=[.03,.1],size=[.96,.85])    
    plot,pix[mask], model[mask], YTITLE='Model', title='Alignment Plots', $
        xr=xrange, yr=[0,1.2],/ys,$
        /xs,_EXTRA=gang_plot_pos(4,1,pltn++, OFFSET=[.03,.1],size=[.96,.85])
        
    
    plot,pix[mask],median((*spec.intensities)[mask],5),ytitle='Spec', $
        xr=xrange,/xs, yr=[0,1.2],/ys,$
       _EXTRA=gang_plot_pos(4,1,pltn++, OFFSET=[.03,.1],size=[.96,.85])
    oplot,pix[mask], model[mask], color=2
    
    
    if keyword_set(datai) gt 0 then begin 
        oplot,[datapt[0:datai-1].x],[datapt[0:datai-1].y],psym=2,color=3
    endif

    info[5]=4*(!y.crange[1]-!y.crange[0])+!y.crange[0]
    info[4]=!y.crange[0]
    info[2]=!x.crange[0]
    info[3]=!x.crange[1]

    ptr_free, telluricSpectrum.intensities, telluricSpectrum.wavelengths, $                                                     -4 pointers
                 synthSpectrum.intensities, synthSpectrum.wavelengths

end


pro picklines, spec, line_centers, wid_wave, wid_pix

    common ysas_common
    maxpts=100
    datapt=replicate({x:0d,y:0d},maxpts)
    pointi=0
    getdata=1b
    info=[0d,3d,0d,0d,0d,0d]
    fit_region=ysas_fit_region
    params=spec.fitparams
    
    points_reqd=n_elements(ysas_wavelength_param_ndxs)+1
    
    
    print,''
    print,'Click a line center on the spectrum, '+$
        'Repeat at least n+1 times.'
    print,'Click off the image when finished.'

    picklineinterface,spec,params,info,fit_region
    repeat begin
        CURSOR,x0,y0,/data,/down

        off= (x0 lt info[2] or x0 gt info[3] or y0 lt info[4] or y0 gt info[5])

        print,x0,y0,off,getdata
        if ~off then begin
            datapt[pointi].x=x0
            datapt[pointi].y=y0
            pointi++
        endif
        off= (pointi ge 2*points_reqd and off)
        
        picklineinterface,spec,params,info,fit_region,datapt,pointi,WAVE=wave
    endrep until off
    
    
    ;return the line centers
    wid_pix=datapt[1:pointi-1:2].x-datapt[0:pointi-1:2].x
    
    cent_pix=round(wid_pix/2.0+datapt[0:pointi-1:2].x)
    wid_pix=round(wid_pix)
    wid_wave=wave[round(datapt[1:pointi-1:2].x)]-wave[round(datapt[0:pointi-1:2].x)]
    line_centers=wave[cent_pix]

end


pro line_locations, order, line_wave, wid_wave, wid_pix

	line_wave_o49=[0.71679313d, 0.71774360, 0.71845511, 0.71863334, 0.71915326, $
		0.71936858, 0.72012327, 0.72043332, 0.72329326, 0.72347507, 0.72407463, $
		0.72645955, 0.72730036, 0.72774024]
	wid_pix_o49=replicate(15,n_elements(line_wave_o49))
	wid_wave_o49=replicate(.00010,n_elements(line_wave_o49))
	
;	line_wave_o50=[0.71679313d, 0.71774360, 0.71845511, 0.71863334, 0.71915326, $
;		0.71936858, 0.72012327, 0.72043332, 0.72329326, 0.72347507, 0.72407463, $
;		0.72645955, 0.72730036, 0.72774024]
;	wid_pix_o50=replicate(15,n_elements(line_wave_o50))
;	wid_wave_o50=replicate(.00010,n_elements(line_wave_o50))
	
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

pro primelambdasoln, spec, order, rv, vsini, first_guess=first_guess, pick_manually=pick_manually


	common ysas_common
	common ysas_debug

	if keyword_set(first_guess) then begin
		;Phase One: Get a halfway decent solution
		;Define a fitting configuration
		psfguess=[3d-5, replicate(0, n_elements(ysas_psf_param_ndxs)-1)]
		init=[simpleWavelengthSolnGuess(spec), $
			  psfguess,  .63, .4, vsini, .6, rv/299792458d, 0, .01]
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
			,scale:[wavescale, psfscale, .1, .1, 1, 25000/299792458d, 0.01, 0.001] $
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
	endif

	;Get the line information
	if keyword_set(pick_manually) then begin
	    ;picklines, spec, line_wave, wid_wave, wid_pix
	    restore, 'line_locations.sav'
	endif else begin
    	line_locations, order, line_wave, wid_wave, wid_pix
	endelse
	
	
	;Compute the model
	telluricSpectrum=getTelluricSpectrum(spec)
	synthSpectrum=getSynthSpectra(spec.synthspec, spec)
	pix=indgen(spec.nPixels)
	params=spec.fitparams
	model = computeSpectralModel(pix, synthSpectrum, $
		telluricSpectrum, params, spectrum = spec, /USE_ZSPACE, $
		pixel_wavelengths = wave, /RETURN_INTENSITIES, telluric_out=tspec)
	
	;Data Array
	obs_cent=dblarr(2, n_elements(line_wave))

	;Pixel space plot for clicking
	!p.multi=[0,2,1]
	getwindow,1,xsize=1000,ysize=500,/erase
	col=2
	
	fit_region_ndxs=rnggen(ysas_fit_region[0],ysas_fit_region[1])
	
	;Fit each of the lines
	npts=0
	for i=0,n_elements(line_wave)-1 do begin
		
		
		;Find the nearest large minimum in the spectrum and get it's pixel
		; position, pix_center
		junk=min(abs(wave[fit_region_ndxs] - line_wave[i]), pix_center)
		pix_center=fit_region_ndxs[pix_center]
		 
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
		plot,pix, (*spec.intensities), $
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
		choice=clickquad()
		if  choice gt 1 then begin
			obs_cent[*,npts]=[fit_pix[1],fit_wave[1]]
			print, obs_cent[*,npts]
			npts++
		endif


	endfor
	!p.multi=0
	
	if npts gt n_elements(ysas_wavelength_param_ndxs) then begin
    	;Compute a wavelength solution
    	obs_cent=transpose(obs_cent)
    	soln=reform(poly_fit(obs_cent[0:npts-1,0],obs_cent[0:npts-1,1],$
    		n_elements(ysas_wavelength_param_ndxs)-1))
        
        spec.fitparams[ysas_wavelength_param_ndxs]=soln

    	print, soln
    	print, mean((wave-computewavelengthsoln(pix,soln))[ysas_fit_region[0]:ysas_fit_region[1]])
	endif else message,'Not enough points to update solution',/info
	
	destroySpec,telluricSpectrum
	destroyspec,synthSpectrum
	destroyspec,tspec
	
	
	
end