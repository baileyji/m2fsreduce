function clickquad
	;clockwise from upperleft starting at 0
	CURSOR,x,y,/down,/normal
	if x le .5 then begin
		return, (y ge .5 ? 0:3)
	endif else begin
		return, (y ge .5 ? 1:2)
	endelse

end