pro make_telluric_file, telluricascii, outfile, domain

    template={commentsymbol:'',datastart:0,delimiter:' ',fieldcount:2l,$
        fieldgroups:[0,1l],fieldlocations:[2,12l],fieldnames:['k','f'],$
        fieldtypes:replicate(5l,2), missingvalue:!values.f_nan, version:1.0}
        
    dat=read_ascii(telluricascii,template=template)
    
    tlamb=reverse(1d3/dat.k)*10 & tfluxorig=reverse(dat.f)>0
    n=where((tlamb gt domain[0]) and (tlamb lt domain[1]))

    wavelengths=tlamb[n]
    flux=tfluxorig[n]
    save,flux,wavelengths,/compress,file=outfile

end


pro make_templogg_file, tloggfile, outfile, domain

    rdsyngf,tloggfile,lamb,flux
    lamb=temporary(lamb)/10000d
    n=where((lamb gt domain[0]) and (lamb lt domain[1]))

    wavelengths=lamb[n]
    flux=flux[n]
    save,flux,wavelengths,/compress,file=outfile
    
end



make_telluric_file, $
    '/Users/One/Dropbox/research/telluricspec/transdata_0.5_1_mic', $
    '/Users/one/m2fs_data/modeldata/telluric_spectrum_0.65-0.75.sav', $
    [0.65,0.75]

make_templogg_file, $
        '/Users/One/Dropbox/research/telluricspec/4rus/lte5500_4.55-0.0.MainSeq.hires1.7.gz', $
        '/Users/one/m2fs_data/modeldata/lte5500_4.55_hires.sav', $
        [0.65,0.75]

end