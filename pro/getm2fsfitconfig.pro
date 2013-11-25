function getM2FSFitConfig, type, tol=tol, wavescalemult=wavescalemult, RVFINE=RVFINE

	common ysas_common


	;Define tolerance
	tol=keyword_set(tol) ? tol:5d-10
	
	wavescalemult = keyword_set(wavescalemult) ? wavescalemult:1d
	
	;Define wavelength, PSF, & RV scales
	;A rough guide is 
	;wavelength shift corresponding to a few pixels / middle pixel ^ power of term,
	; then use that power 
	wavescale=[1d-5, 1d-8, 1d-12, 1d-15, 1d-18, 1d-21, 1d-25, 1d-28]
	wavescale*=wavescalemult
	
	psfsigmascale=1d-5 ; about 2 pixels
	psfhieghtscale=replicate(.1,8)
	
	vsiniscale=5
	
	rvscale=25000/299792458d
    rvscale = keyword_set(RVFINE) ? RVFINE/299792458d: rvscale
	
	vsini_ub=25
    vsini_lb=-25
	
	waveub=[1,1,1,1,1,1,1,1]
	wavelb=[-1,-1,-1,-1,-1,-1,-1,-1]
	
	psfub=replicate(100,n_elements(ysas_psf_param_ndxs))
	psflb=replicate(-100,n_elements(ysas_psf_param_ndxs))
	
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
      
    fixed[*]=1
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

    fixed[*]=1
    fixed[ysas_rotational_broadening_param_ndxs[0]]=0
    fixed[ysas_rv_param_ndx]=0
    fixed[ysas_veiling_param_ndx]=0
    fixed[ysas_norm_offset_param_ndxs]=0
    nocalibconfig={fixed:fixed $
        ,ub:[2,   vsini_ub,  .0005d,   0.5, 40000]  $
        ,lb:[0,   vsini_lb, -.0005d,   -.5,-40000]  $
        ,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
        ,scale:[.1, vsiniscale, rvscale, 0.01, 0.001] $
        ,mode:'Final' $
        ,fftenable:0b $
        ,fromexisting:1b $
        ,maxreps:10 $
        ,tol:tol $
        ,vsini:0d}


		
	fixed[*]=0
	fixed[ysas_rotational_broadening_param_ndxs[1]]=1
	allconfig={fixed:fixed $
		,ub:[waveub, psfub,   2, 2,   vsini_ub,  .0005d,   0.5, 40000] $
		,lb:[wavelb, psflb,   0, 0,   vsini_lb, -.0005d,   -.5,-40000] $
		,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
		,scale:[wavescale, psfscale, .1, .1, vsiniscale, rvscale, 0.01, 0.001] $
		,mode:'Final' $
		,fftenable:0b $
		,fromexisting:1b $
		,maxreps:10 $
		,tol:tol $
		,vsini:0d}

    fixed[*]=0
	fixed[ysas_rv_param_ndx]=1
	fixed[ysas_rotational_broadening_param_ndxs[1]]=1
	norvconfig={fixed:fixed $
		    ,ub:[waveub, psfub,   2, 2,   vsini_ub,    0.5, 40000] $
		    ,lb:[wavelb, psflb,   0, 0,   vsini_lb,   -.5,-40000] $
		    ,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
		    ,scale:[wavescale, psfscale, .1, .1, vsiniscale, 0.01, 0.001] $
		    ,mode:'Final' $
		    ,fftenable:0b $
		    ,fromexisting:1b $
		    ,maxreps:10 $
		    ,tol:tol $
		    ,vsini:0d}

	fixed[*]=0
	fixed[ysas_wavelength_param_ndxs]=1
	fixed[ysas_rotational_broadening_param_ndxs[1]]=1
	nowaveconfig={fixed:fixed $
	    ,ub:[ psfub,   2, 2,   vsini_ub,  .0005d,   0.5, 40000] $
	    ,lb:[ psflb,   0, 0,   vsini_lb, -.0005d,   -.5,-40000] $
	    ,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
	    ,scale:[psfscale, .1, .1, vsiniscale, rvscale, 0.01, 0.001] $
	    ,mode:'Final' $
	    ,fftenable:0b $
	    ,fromexisting:1b $
	    ,maxreps:10 $
	    ,tol:tol $
	    ,vsini:0d}
	
case type of
	
	'all':return, allconfig
	'wave':return,wavepconfig
	'psf':return,psfconfig
	'wavepsf':return,wavepsfconfig
	'nostellar':return,nostellarconfig
    'nocalib':return,nocalibconfig
	'init':return,initconfig
	'norv':return,norvconfig
	'nowave':return,nowaveconfig
	else: message,'Invalid type'
	
endcase	

		
end