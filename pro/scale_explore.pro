;Test to explore scale & tolerance effects on pinima
;
;

conf=getM2FSfitconfig('all', tol=1d-6,/rvfine)

tols=reverse(10^(-11+dindgen(8)))

rvscales=[1, 10, 100, 1000, 10000]/299792458d


;Get extracted spectrum
extractedSpec=*spectra[2]

;Get telluric spectrum to use in fitting (+2 pointers)
telluricSpec=getTelluricSpectrum(extractedSpec, USE_FIT=conf.fromExisting)

;Get the synthetic spectra (+2 pointers)
synthSpec=getSynthSpectra(temp_logg, extractedSpec, USE_FIT=conf.fromExisting)


items=[[tols],[rvscales]]

results=dblarr(size(items,/dim)[0],size(items,/dim)[1],
            ysas_max_num_fit_params+1)

for i=0, size(items,/dim)[0]-1 do begin
    for j=0, size(items,/dim))[1]-1 do begin
    
        ;Reset the configuration
        conf.params=extractedSpec.fitparams
        
        ;Unpack Items
        conf.params[ysas_rv_param_ndx]=items[i,1]
        conf.tol=items[j,2]
        
        ;Do the fit
        optimalParameters=fit_spectrum(extractedSpec, synthSpec, telluricSpec, $
            conf.fixed, conf.params, VSCALE=conf.scale, CHISQR=chisqr, $
            NCALLS=numCalls, TOLARANCE=conf.tol, UPPERBOUNDS=conf.ub, LOWERBOUNDS=conf.lb,$
            RESIDUALS=residuals, VARIABLE_THRESH=1d, CHISQUARE_THRESH=1, maxreps=1, $
            FIT_REGION=[350,3000])


        ;Store results
        results[i,j,*]=[chisqr,optimalparameters]
        
   endfor
endfor


isurface, results[*,*,0]

end