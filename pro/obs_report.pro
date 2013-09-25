datapath='IDLWorkspace82/ut20130826/'
fits=file_search(datapath+'*c1.fits')
for i=0l, n_elements(fits) -1 do begin
  f=fits[i]
  print, file_basename(f)
  print, string(9b), fxpar(headfits(f),'OBJECT')
  print, string(9b), fxpar(headfits(f),'COMMENT') 
  print, string(9b), fxpar(headfits(f),'EXPTYPE')
  print, string(9b), fxpar(headfits(f),'BINNING')
endfor


end