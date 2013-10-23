function readMarioSpecFile, file, order, TEMP_LOGG=TEMP_LOGG, $
	IMAGE_FILE=IMAGE_FILE, FIBER=FIBER, OBJNAME=OBJNAME, NOM_VSINI=NOM_VSINI, $
	NOM_RV=NOM_RV, PTR=PTR

	template={datastart:0, delimiter:32b,fieldcount:2,fieldgroups:[0,1],$
		fieldlocations:[0,4],fieldnames:['pixel','flux'],fieldtypes:[2,5],$
		missingvalue:!values.d_nan,version:1.,commentsymbol:'#'}

	data=read_ascii(file, template=template) ;49
	nspec=buildEmpiricSpectrum(data.flux,  order-1, NORMALIZE=3, $
		EXTRACTID='MARIO', instrument='M2FS_hires')
	
	if keyword_set(TEMP_LOGG) then $
		nspec.synthspec=TEMP_LOGG

	spec={m2fsspectrum}
	COPY_STRUCT, nspec, spec
	
	spec.specfile=file
	spec.order=49
	
	if keyword_set(IMAGE_FILE) then begin
		spec.imfile=IMAGE_FILE
		spec.header=ptr_new(headfits(IMAGE_FILE))
		spec.ysasheader=ysas_modeinfo(IMAGE_FILE)
	endif
	if keyword_set(FIBER) then spec.fiber=FIBER
	if keyword_set(OBJNAME) then spec.objname=OBJNAME
	if keyword_set(NOM_VSINI) then spec.nom_vsini=NOM_VSINI
	if keyword_set(NOM_RV) then spec.nom_rv=NOM_RV

	if keyword_set(ptr) then return, ptr_new(spec,/no_copy) $
	else return, spec

end

function loadAug13Spectra

common ysas_common

datadir='/Users/one/Desktop/M2FS Commissioning Run/'

spectra=ptrarr(3+7*6*2+6*2)
ndx=0

ysas_fit_region=[350,3000]

;Grab RV standard spectra
sfile=datadir+'spectra/r0089.ec_2.txt'
spectra[ndx]=readMarioSpecFile(sfile, 49, $
	temp_logg='5500_4.55', $
	image_file=datadir+'image/r0089.fits', $
	fiber='R5-02', $
	objname='HIP 10798', $
	nom_vsini=3.26d, $
	nom_rv=7000,/PTR)
(*spectra[ndx]).fiberno=5*16+2
ndx+=1

sfile=datadir+'spectra/r0090.ec_2.txt'
spectra[ndx]=readMarioSpecFile(sfile, 49, $
	temp_logg='5500_4.55', $
	image_file=datadir+'image/r0090.fits', $
	fiber='R5-02', $
	objname='HIP 10798', $
	nom_vsini=3.26d, $
	nom_rv=7000,/PTR)
(*spectra[ndx]).fiberno=5*16+2
ndx+=1

sfile=datadir+'spectra/r0109.ec_2.txt'
spectra[ndx]=readMarioSpecFile(sfile, 49, $
   temp_logg='5500_4.55', $
   image_file=datadir+'image/r0109.fits', $
   fiber='R5-02', $
   objname='HIP 10798', $
   nom_vsini=3.26d, $
   nom_rv=7000,/PTR)
(*spectra[ndx]).fiberno=5*16+2
ndx+=1


;Load the R side data
	
	order=49
	rORb='r'
	fnums=rnggen(97,103)
	stars=replicate('',12)
	vsinis=replicate(3d,12)
	rvs=replicate(0d,12)
	temp_loggs=replicate('5500_4.55',12)
	fiberIDs=replicate('',12)
	
	;Grab individual R spectra		
	for fnum=0, n_elements(fnums)-1 do begin
		for anum=2,12,2 do begin
			
			fnum_s=string(fnums[fnum],format='(i04)')
			ifile=datadir+'image/'+rORb+fnum_s+'.fits'
			sfile=datadir+'spectra/'+rORb+fnum_s+'.ec2_'+strtrim(anum,2)+'.txt'
			
			spectra[ndx]=readMarioSpecFile(sfile, order, $
				temp_logg=temp_loggs[anum-1], $
				image_file=ifile, $
				fiber=fiberIDs[anum-1], $
				objname=stars[anum-1], $
				nom_vsini=vsinis[anum-1], $
				nom_rv=rvs[anum-1],/PTR)
			(*spectra[ndx]).fiberno=anum
			ndx+=1
				
		endfor
	endfor
	
	;Grab coadded Spectra
	for anum=2,12,2 do begin
		sfile=datadir+'spectra/blanco1_'+rORb+'.ec2_'+strtrim(anum,2)+'.txt'

		spectra[ndx]=readMarioSpecFile(sfile, order, $
			temp_logg=temp_loggs[anum-1], $
			fiber=fiberIDs[anum-1], $
			objname=stars[anum-1], $
			nom_vsini=vsinis[anum-1], $
			nom_rv=rvs[anum-1],/PTR)
		(*spectra[ndx]).fiberno=anum
		(*spectra[ndx]).ysasheader.baryv=(*spectra[anum/2+2]).ysasheader.baryv
		ndx+=1
	endfor
			
			

;Load the B side data
	order=49
	rORb='b'
	fnums=rnggen(77,83)
	stars=replicate('',12)
	vsinis=replicate(3d,12)
	rvs=replicate(0d,12)
	temp_loggs=replicate('5500_4.55',12)
	fiberIDs=replicate('',12)
	

	;Grab individual spectra
	for fnum=0, n_elements(fnums)-1 do begin
		for anum=2,12,2 do begin
		
			fnum_s=string(fnums[fnum],format='(i04)')
			ifile=datadir+'image/'+rORb+fnum_s+'.fits'
			sfile=datadir+'spectra/'+rORb+fnum_s+'.ec2_'+strtrim(anum,2)+'.txt'
			
			spectra[ndx]=readMarioSpecFile(sfile, order, $
				temp_logg=temp_loggs[anum-1], $
				image_file=ifile, $
				fiber=fiberIDs[anum-1], $
				objname=stars[anum-1], $
				nom_vsini=vsinis[anum-1], $
				nom_rv=rvs[anum-1],/PTR)
			(*spectra[ndx]).fiberno=anum+128
			ndx+=1
		endfor
	endfor

	;Grab coadded Spectra
	for anum=2,12,2 do begin
		sfile=datadir+'spectra/blanco1_'+rORb+'.ec2_'+strtrim(anum,2)+'.txt'

		spectra[ndx]=readMarioSpecFile(sfile, order, $
			temp_logg=temp_loggs[anum-1], $
			fiber=fiberIDs[anum-1], $
			objname=stars[anum-1], $
			nom_vsini=vsinis[anum-1], $
			nom_rv=rvs[anum-1],/PTR)
		(*spectra[ndx]).fiberno=anum+128
		(*spectra[ndx]).ysasheader.baryv=(*spectra[anum/2+2]).ysasheader.baryv
		ndx+=1
	endfor

	return, spectra

end