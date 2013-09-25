function get_fits, datapath, flat=flat, bias=bias, object=object, quadrant=quadrant
  fits=file_search(datapath+'*.fits')
  
  BIAS_TYPE_STR='Bias'
  FLAT_TYPE_STR='Flat'
  OBJECT_TYPE_STR='Object'
  bias_map=bytarr(n_elements(fits))
  flat_map=bytarr(n_elements(fits))
  object_map=bytarr(n_elements(fits))
  
  for i=0, n_elements(fits)-1L do begin
    bias_map[i]=fxpar(headfits(f),'EXPTYPE') eq BIAS_TYPE_STR
    flat_map[i]=fxpar(headfits(f),'EXPTYPE') eq FLAT_TYPE_STR
    object_map[i]=fxpar(headfits(f),'EXPTYPE') eq OBJECT_TYPE_STR
  endfor
  
  if keyword_set(flat) then ret=fits[where(flat_map)]
  if keyword_set(bias) then ret=fits[where(bias_map)]
  if keyword_set(object) then ret=fits[where(object_map)]
  
  if keyword_set(quadrant) then begin
    quad=where(strposa(ret,quadrant+'.fits') ne -1,count)
    if count eq 0 then return, !NULL
    return, ret[quad]    
  endif
end


  cd,'~'
  datapath='IDLWorkspace82/ut20130826/'
  
  frames=file_basename(file_search(datapath+'*c1.fits'),'c1.fits')
  
  quadrants=['c1','c2','c3','c4']
  nx=4096
  ny=4112
  
  overscan=[[2048,2175],[2048,2175],[2048,2175],[2048,2175]]
  crop=[[0,2047,0,2055], $
        [0,2047,0,2055], $
        [0,2047,0,2055], $
        [0,2047,0,2055]]


  
  quad_pos=[[0, 2047, 0, 2055, 0],$
            [2048, 4095, 0, 2055, 5],$
            [2048, 4095, 2056, 4111, 2],$
            [0, 2047, 2056, 4111, 7]]

  for i=0l, n_elements(frames)-1 do begin
      i=(where(frames eq 'r0090'))[0]
      frame=frames[i]
      
      outfile=datapath+frame+'.fits'
      
      im=intarr(nx,ny)
      
      ;Process each quadrant 
      for quad=0, 3 do begin
        ;load, overscan correct, & trim each quadrant
        
        file=datapath+frame+quadrants[quad]+'.fits'
        
        im_quad=readfits(file, hdr, /silent)
        
        ;Correct column bias and crop
        bias_corr=colbias(im_quad, $
          overscan[0,quad], overscan[1,quad], $
          crop[0,quad], crop[1,quad], crop[2,quad], crop[3,quad], $
          biasval=biasval, noiseval=noiseval)
        
;        x=0 & w=100 & h=1500
;        bias_corr[x:x+w,x+10:x+10+w]=h
;        x=1000 & w=100
;        if quad ge 1 then bias_corr[x:x+w,x:x+w]=h
;        x=1200 & w=100
;        if quad ge 2 then bias_corr[x:x+w,x:x+w]=h
;        x=1400 & w=100
;        if quad ge 3 then bias_corr[x:x+w,x:x+w]=h

        im[ quad_pos[0,quad]:quad_pos[1,quad], $
            quad_pos[2,quad]:quad_pos[3,quad]] = rotate(bias_corr, $
                                                        quad_pos[4,quad])
        
        print, string(9b)+quadrants[quad]+" biasv: "+strtrim(biasval,2)+"   noisev: "+strtrim(noiseval,2)
        
     endfor


     
     ;Update header with new dims
     FXADDPAR,hdr, 'NAXIS1', nx
     FXADDPAR,hdr, 'NAXIS1', ny
     FXADDPAR,hdr, 'CHOFFX',0.0d
     FXADDPAR,hdr, 'CHOFFY',0.0d
     FXADDPAR,hdr, 'DATASEC','[1:4096,1:4112]'
     FXADDPAR,hdr, 'TRIMSEC','[1:4096,1:4112]'
     FXADDPAR,hdr, 'BIASSEC',''
     ;Write out the file
     writefits, outfile, im, hdr
     stop
     
  endfor

end
