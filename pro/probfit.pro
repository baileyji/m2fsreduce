; probfit

ysas_startup_m2fs
common ysas_common
common ysas_debug

;Setup the pipline
ysas_extract_ID='MARIO'
ysas_fit_id='ysastest'
ysas_fit_region=[350,3000]
ysasdb_verbose_fit=1000
ysasdb_fitting_progress_plots=[0,0,0]

restore, 'm2fs_data/spectra.sav'

spectra=spectra[0:2]

o49psfheights=[0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0, 0.00]
;Right1,Right2,Right3,Right4 Far right!,Left1,Left2,Left3,Left4  Far left!
o49psf=[3d-5, o49psfheights]

;for i=0,n_elements(spectra)-1 do begin
;    (*spectra[i]).nom_rv=7459d
;    (*spectra[i]).fitparams[ysas_rv_param_ndx]=((*spectra[0]).nom_rv - $
;                                                (*spectra[0]).ysasheader.baryv)/$
;                                                299792458d
;                                           
;    (*spectra[i]).fitparams[ysas_norm_offset_param_ndxs]=[0,0]
;    (*spectra[i]).fitparams[ysas_psf_param_ndxs]=o49psf
;    (*spectra[i]).synthspec='5400-4.50-0.5.PHOENIX.sav'
;
;endfor
;stop
;snrdofit[*]=0
;snrdofit[0:2]=1
;fitSpectra, snrspectra, '', '5400-4.50-0.5.PHOENIX.sav', '', results, $
;    CONFIG=getM2FSfitconfig('init', tol=1d-6), STATUS=STATUS, $
;    CHISQUARE_LIMIT=30000d+200, MAX_FITSPECTRA_CYCLES=5, EXCLUDE=~snrdofit, $
;    CHISQUARE_THRESH=1d
;spectra=snrspectra[0:2]
;save,spectra,/compress,file='./m2fs_data/probfit_exp_start1d-6_psfgt0.sav'


;Run this first, MAX_FITSPECTRA_CYCLES needs to be enough for it to converge and 200 isn't enough
; Be sure to update and name the save file

restore,'./m2fs_data/probfit_exp_start1d-6_psfgt0.sav'
fitSpectra, spectra, '', '5400-4.50-0.5.PHOENIX.sav', '', results, $
    CONFIG=getM2FSfitconfig('norv', tol=1d-10), STATUS=STATUS, $
    CHISQUARE_LIMIT=30000d+200, $
    MAX_FITSPECTRA_CYCLES=500, $
    CHISQUARE_THRESH=1d-7,$
    PARAMLOGDATA=PARAMLOGDATA, /NO_VARIABLE_CHANGE_CHECK
save,spectra,PARAMLOGDATA,/compress,file='./m2fs_data/probfit_phase1_1d-10_500_psfgt0.sav'
ptr_free,PARAMLOGDATA.paramlog
ptr_free,PARAMLOGDATA.paramndxs
;restore,'./m2fs_data/probfit_phase1_1d-8_500_psfgt0.sav'

;fitSpectra, spectra, '', '5400-4.50-0.5.PHOENIX.sav', '', results, $
;    CONFIG=getM2FSfitconfig('nowave', tol=1d-8,rvfine=100), STATUS=STATUS,$
;    CHISQUARE_LIMIT=30000d+200, $
;    MAX_FITSPECTRA_CYCLES=50,$
;    CHISQUARE_THRESH=1d-7;,SCALE_CHANGE=.5  
;save,spectra,/compress,file='./m2fs_data/probfit_exp_phase2=nowave_1d-8_150_psfgt0.sav'

;Run this second, MAX_FITSPECTRA_CYCLES needs to be enough for it to converge and 200 isn't enough
; Be sure to update and name the save file
fitSpectra, spectra, '', '5400-4.50-0.5.PHOENIX.sav', '', results, $
    CONFIG=getM2FSfitconfig('all', tol=1d-10,rvfine=10), STATUS=STATUS,$
    CHISQUARE_LIMIT=30000d+200, $
    MAX_FITSPECTRA_CYCLES=500,$
    CHISQUARE_THRESH=1d-7, PARAMLOGDATA=PARAMLOGDATA,/NO_VARIABLE_CHANGE_CHECK
save,spectra,PARAMLOGDATA,/compress,file='./m2fs_data/probfit_phase2=all_1d-10_500_psfgt0.sav'
;restore,'./m2fs_data/probfit_phase2=all_1d-8_500_psfgt0.sav'
rv=dblarr(3)
for i=0,2 do rv[i]=(*spectra[i]).fitparams[ysas_rv_param_ndx]*299792458d + (*spectra[i]).ysasheader.baryv
print,rv
        
end