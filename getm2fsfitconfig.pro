function getM2FSFitConfig, type, vsini,tol=tol

	common ysas_common


	;Define tolerance
	tol=keyword_set(tol) ? tol:5d-10
	
	;Define wavelength, PSF, & RV scales
	wavescaleloose=[.008, 1d-7, 1d-11, 1d-13,1d-16,1d-19,1d-23,1d-27]
	wavescale=wavescaleloose/10
	
	psfsigmascale=1d-5
	psfhieghtscale=replicate(.01,8)
	
	rvscale=25000/3d8
	
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
		,scale: wavescaleloose[0:2] $
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
		
	fixed*=0
	fixed[ysas_rotational_broadening_param_ndxs[1]]=1
	allconfig={fixed:fixed $
		,ub:[waveub, psfub,   2, 2,   25,  .0005d,   0.5, 40000] $
		,lb:[wavelb, psflb,  .1,.1,   .1, -.0005d,   -.5,-40000] $
		,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
		,scale:[wavescale, psfscale, .1, .1, 1, rvscale, 0.01, 0.001] $
		,mode:'Final' $
		,fftenable:0b $
		,fromexisting:1b $
		,maxreps:3 $
		,tol:tol $
		,vsini:0d}

	fixed*=0
	fixed[ysas_rotational_broadening_param_ndxs[1]]=1
	alllooseconfig={fixed:fixed $
		,ub:[waveub, psfub,   2, 2,   25,  .0005d,   0.5, 40000] $
		,lb:[wavelb, psflb,  0, 0,   .001, -.0005d,   -.5,-40000] $
		,params: dblarr(ysas_MAX_NUM_FIT_PARAMS) $
		,scale:[wavescaleloose, psfscale, .1, .1, 1, rvscale, 0.01, 0.001] $
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
	'loose':return,alllooseconfig
	'init':return,initconfig
	else: message,'Invalid type'
	
endcase	

		
end