function getM2FSFitConfig, type, tol=tol, wavescalemult=wavescalemult, RVFINE=RVFINE

	common ysas_common


	;Define tolerance
	tol=keyword_set(tol) ? tol:5d-10
	
	wavescalemult = keyword_set(wavescalemult) ? wavescalemult:0.1d
	
	;Define wavelength, PSF, & RV scales
	wavescale=[1d-2, 1d-7, 1d-11, 1d-13, 1d-16, 1d-19, 1d-23, 1d-27]
	wavescale*=wavescalemult
	
	psfsigmascale=1d-5
	psfhieghtscale=replicate(.01,8)
	
	rvscale=25000/3d8
    rvscale = keyword_set(RVFINE) ? 10/3d8: rvscale
	
	
	waveub=[1,1,1,1,1,1,1,1]
	wavelb=[-1,-1,-1,-1,-1,-1,-1,-1]
	
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
		,vsini:0d}
		
	fixed[*]=1
	fixed[ysas_wavelength_param_ndxs[0:2]]=0
	initconfig={$
		fixed:fixed $
		,ub:waveub[0:2] $
		,lb:wavelb[0:2] $
		,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
		,scale: wavescale[0:2] $
		,mode:'Final' $
		,fftenable:0b $
		,fromexisting:1b $
		,maxreps:10 $
		,tol:tol $
		,vsini:0d}
		
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
		,vsini:0d}
		
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
		,vsini:0d}
		
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
        ,vsini:0d}
        
    fixed[ysas_wavelength_param_ndxs]=0
    fixed[ysas_psf_param_ndxs]=0
    fixed[ysas_airmass_param_ndx]=0
    fixed[ysas_norm_offset_param_ndxs]=0
    nostellarconfig={fixed:fixed $
        ,ub:[waveub,   psfub, 2,   0.5, 40000 ] $
        ,lb:[wavelb,   psflb, 0,  -0.5, -40000 ] $
        ,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
        ,scale:[wavescale, psfscale, .1, 0.01, 0.001] $
        ,mode:'Final' $
        ,fftenable:0b $
        ,fromexisting:1b $
        ,maxreps:10 $
        ,tol:tol $
        ,vsini:0d}        

		
	fixed*=0
	fixed[ysas_rotational_broadening_param_ndxs[1]]=1
	allconfig={fixed:fixed $
		,ub:[waveub, psfub,   2, 2,   25,  .0005d,   0.5, 40000] $
		,lb:[wavelb, psflb,   0, 0,   .001, -.0005d,   -.5,-40000] $
		,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
		,scale:[wavescale, psfscale, .1, .1, 1, rvscale, 0.01, 0.001] $
		,mode:'Final' $
		,fftenable:0b $
		,fromexisting:1b $
		,maxreps:3 $
		,tol:tol $
		,vsini:0d}

		
	
case type of
	
	'all':return, allconfig
	'wave':return,wavepconfig
	'psf':return,psfconfig
	'wavepsf':return,wavepsfconfig
	'nostellar':return,nostellarconfig
	'init':return,initconfig
	else: message,'Invalid type'
	
endcase	

		
end