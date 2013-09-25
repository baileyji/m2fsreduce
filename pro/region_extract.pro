

function reg_extract, type,unit, low, high

  if strmatch(type, 'lte5500', /FOLD_CASE) then begin
    if strmatch(unit,'u') then begin
       low*=10000d
       high*=10000d
    endif
    file='/Users/One/Dropbox/research/telluricspec/4rus/lte5500_4.55-0.0.MainSeq.hires1.7.gz'
    rdsyngf,file,lamb,flux
    n=where((lamb gt low) and (lamb lt high))
    return, {wave:lamb[n], flux:flux[n]}
  endif
  if strmatch(type, 'telluric', /FOLD_CASE) then begin
    file='/Users/One/Dropbox/research/telluricspec/transdata_0.5_1_mic'
    template={commentsymbol:'',datastart:0,delimiter:' ',fieldcount:2l,$
      fieldgroups:[0,1l],fieldlocations:[2,12l],fieldnames:['k','f'],$
      fieldtypes:replicate(5l,2), missingvalue:!values.f_nan, version:1.0}
    dat=read_ascii(file,template=template)
    tlamb=reverse(1d3/dat.k)*10 & tfluxorig=reverse(dat.f)>0
    n=where((tlamb gt low) and (tlamb lt high))
    return, {wave:tlamb[n], flux:tflux[n]}
  endif


end 
  
  
  
d=reg_extract('lte5500','u',.710517d,716517d)

plot,d.wave,d.flux


end