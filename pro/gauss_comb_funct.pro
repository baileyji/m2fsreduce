
pro gauss_comb_funct,x,b,f,pder
    compile_opt idl2
    
  zthresh=1.0d-100
  
  f=make_array(n_elements(x),value=0d,/double)
  
  ;if n_params(0) gt 3 then pder = dblarr(n_elements(x), n_elements(b))
  
  for i=0l, n_elements(b)/3-1 do begin
      a=b[3*i:3*i+2]
      f+=gauss1(x, a , /peak)
      
;      z = (x-a[0])/(a[1] > zthresh)   ;get z
;      ez = exp(-z^2/2d)
;      f+=a[2]*ez
;      
;      if n_params(0) gt 3 then begin
;          pder[*,3*i] = ez      ;compute partials
;          
;          pder[*,3*i+1] = a[0] * ez * z/(a[2] > zthresh)
;          pder[*,3*i+2] = pder[*,1] * z
;        
;      endif
      
  endfor

end

function gauss_comb_funct, x, p, dp
  gauss_comb_funct,x, p, f, dp
  return,f
end