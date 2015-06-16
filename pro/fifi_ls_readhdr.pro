function fifi_ls_readhdr, header, CHECKREQ=checkreq, CANCEL=cancel

  cancel = 0
  if n_elements(checkreq) eq 0 then checkreq = 1 $
  else checkreq = keyword_set(checkreq)
  hdr = hash()

  pkgpath = file_dirname(file_dirname(file_which('fifi_ls_readhdr.pro')),/MARK)
  datapath = filepath('data', ROOT_DIR=pkgpath) + path_sep()
  
  ; read header definition table
  readcol, datapath + 'headerdef.dat', $
           keyword, req, default, type, min, max, enum, $
           COMMENT='#', FORMAT='A,A,A,A,A,A,A', /SILENT
  type = strlowcase(strtrim(type,2))

  ; read header comment table
  readcol, datapath + 'headercomment.dat', $
           keys, default_comment, DELIMITER=',', $
           COMMENT='#', FORMAT='A,A', /SILENT
  keys = strupcase(strtrim(keys,2))
  default_comment = ' '+strtrim(default_comment,2)
  comments = hash(keys, default_comment)

  ; first just populate structure with values from header
  for i=0,n_elements(keyword)-1 do begin
     key = strupcase(strtrim(keyword[i],2))
     val = sxpar(header, key, COMMENT=comment, COUNT=found)
     if n_elements(comment) eq 0 then begin
        comment = ' ' 
     endif else begin
        comment = ' '+strtrim(comment,2)
     endelse

     default_val = default[i]
     if type[i] eq 'float' then default_val = float(default_val)
     if type[i] eq 'int' then default_val = fix(default_val)

     if found gt 0 then begin
        if checkreq && strupcase(req[i]) eq 'Y' then begin
           if type[i] eq 'float' || type[i] eq 'int' then begin

              ; check for valid numbers
              if ~valid_num(val) then begin
                 print, 'Required keyword '+key+$
                        ' has wrong type (value: '+strtrim(val,2)+$
                        ').  Should be '+type[i]
                 cancel = 1 
                 hdr[key] = default_val
                 continue
              endif

              ; convert to appropriate value
              if type[i] eq 'float' then val = float(val)
              if type[i] eq 'int' then val = fix(val)

              ; check for min/max
              if min[i] ne '.' && val le float(min[i]) then begin
                 print, 'Required keyword '+key+$
                        ' has wrong value '+strtrim(val,2)+$
                        '.  Should be >'+min[i]
                 cancel = 1 
                 hdr[key] = default_val
                 continue
              endif
              if max[i] ne '.' && val ge float(max[i]) then begin
                 print, 'Required keyword '+key+$
                        ' has wrong value '+strtrim(val,2)+$
                        '.  Should be <'+max[i]
                 cancel = 1 
                 hdr[key] = default_val
                 continue
              endif

           endif else begin
              if type[i] eq 'bool' then begin
                 if size(val,/TYPE) ne 1 then begin
                    print, 'Required keyword '+key+$
                           ' has wrong type (value: '+strtrim(val,2)+$
                           ').  Should be '+type[i]
                    cancel = 1 
                    hdr[key] = default_val
                    continue
                 endif

              endif else begin

                 val = strtrim(val,2)

                 ; check for enum
                 if enum[i] ne '.' then begin
                    opt = strupcase(strsplit(enum[i],'|',/EXTRACT))
                    idx = where(strupcase(val) eq opt)
                    if idx[0] eq -1 then begin
                       print, 'Required keyword '+key+$
                              ' has wrong value '+val+$
                              '.  Should be within '+enum[i]
                       cancel = 1 
                       hdr[key] = default_val
                       continue
                    endif
                 endif
              endelse
           endelse
        endif

        if type[i] eq 'bool' then begin
           ; convert byte to 'T' or 'F'
           if val eq 1B then val = 'T' else val = 'F'
        endif
        if type[i] eq 'string' then val = strtrim(val,2)

        hdr[key] = val
        if comment ne ' ' then begin
           comments[key] = comment
        endif
     endif else begin
        if checkreq && strupcase(req[i]) eq 'Y' then begin
           print, 'Required keyword '+key+' not found'
           cancel = 1
        endif
        hdr[key] = default_val
     endelse
  endfor

  ; add comments hash to header hash 
  hdr['COMMENTS'] = comments

  ; now calculate values needed from header

  ; paths
  hdr['PKGPATH'] = pkgpath
  hdr['DATAPATH'] = datapath

  ; decimal date
  bad_date = 1
  if hdr['DATE-OBS'] ne 'UNKNOWN' then begin
     date_arr = strsplit(hdr['DATE-OBS'],'-T',/EXTRACT)
     if n_elements(date_arr) ge 3 then begin
        y = date_arr[0]
        m = date_arr[1]
        d = date_arr[2]
        if valid_num(y) && valid_num(m) && valid_num(d) then begin
           y = fix(y)
           m = fix(m)
           d = fix(d)
           if y gt 100 then y-=2000
           date = y+m/100.+d/1.0e+04
           bad_date = 0
           hdr['FDATE'] = date
        endif
     endif
  endif
  if bad_date then begin
     print,'fifi_ls_readhdr: DATE-OBS '+hdr['DATE-OBS']+$
           ' not understood, using date='+$
           strtrim(hdr['FDATE'],2)
  endif
  
  ; cardmode
  ; standardize values to upper case
  hdr['CARDMODE'] = strupcase(hdr['CARDMODE'])

  ; instcfg
  ; standardize values to upper case
  hdr['INSTCFG'] = strupcase(hdr['INSTCFG'])
  case hdr['INSTCFG'] of
     'HI-MED': hdr['INSTCFG'] = 'HIGH_MED'
     'HIMED': hdr['INSTCFG'] = 'HIGH_MED'
     'HI-LO': hdr['INSTCFG'] = 'HIGH_LOW'
     'HILO': hdr['INSTCFG'] = 'HIGH_LOW'
     'MED': hdr['INSTCFG'] = 'MEDIUM'
     'LO': hdr['INSTCFG'] = 'LOW'
     'CAM': hdr['INSTCFG'] = 'CAMERA'
     'PUP': hdr['INSTCFG'] = 'CAMERA'
     else: ; leave alone
  endcase

  ; obsmode
  ; standardize values to upper case
  ; also track whether is marked as fowler
  hdr['OBSMODE'] = strupcase(hdr['OBSMODE'])
  case hdr['OBSMODE'] of
     'FOWLER': begin
        hdr['FOWLER'] = 1
        hdr['OBSMODE'] = 'NOD'
     end
     'CHOP-NOD': hdr['OBSMODE'] = 'CHOPNOD'
     'DARK': hdr['OBSMODE'] = 'STARE'
     else: ; leave alone
  endcase

  ; get tort parameters from config file
  ; columns are:
  ;   date hrfl0 xdfl0 slitrot krot brl x0brl y0brl
  readcol, datapath + 'tortparm.dat', $
           dt,hr,xdf,sr,kr,br,x0,y0, $
           COMMENT='#', /SILENT
  i = (where(hdr['FDATE'] lt dt))[0]
  hdr['HRFL0'] = hr[i]
  hdr['XDFL0'] = xdf[i]
  hdr['SLITROT'] = sr[i]
  hdr['BRL'] = br[i]
  hdr['X0BRL'] = x0[i]
  hdr['Y0BRL'] = y0[i]
  ; only use default for krot if it is not already set
  if hdr['KROT'] eq -9999. then hdr['KROT'] = kr[i]

  ; get grating angle by mode
  instcfg = hdr['INSTCFG']
  if instcfg eq 'MEDIUM' || instcfg eq 'HIGH_MED' then begin
     hdr['GRATANGL'] = hdr['ECHELLE']
     hdr['XDDGR'] = hdr['XDMRDGR']
  endif
  if instcfg eq 'LOW' || instcfg eq 'HIGH_LOW' then begin
     hdr['GRATANGL'] = hdr['LORES']
     hdr['XDDGR'] = hdr['XDLRDGR']
  endif

  ; xdr from wave number, iorder, xddgr
  
  if instcfg ne 'CAMERA' then begin
     waveno0 = hdr['WAVENO0']
     xddgr = hdr['XDDGR']
;;;;here - sanity check for iorder, waveno
     iorder = round(2.0 * xddgr * sin(hdr['GRATANGL']/!RADEG) * waveno0)
     sinang = iorder/(2.0 * xddgr * waveno0)

     theta = asin(sinang)
     hdr['XDR'] = tan(theta+0.025)

  endif

  ; set addtime and cumtime to exptime to begin with
  if hdr['ADDTIME'] eq -9999. then hdr['ADDTIME'] = hdr['EXPTIME']
  if hdr['CUMTIME'] eq -9999. then hdr['CUMTIME'] = hdr['EXPTIME']

  ; slitval
  ; convert slit angle from degrees to cm
  slit = hdr['SDEG']
  readcol, datapath + 'slitval.dat', dt, lim, d1, d2, c1, c2,$
           v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,$
           v13,v14,v15,v16,v17,v18,v19,v20,v21,v22,v23,v24,$
           COMMENT='#', /SILENT

  i = (where(hdr['FDATE'] lt dt))[0]
  val = [v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],v7[i],v8[i],$
         v9[i],v10[i],v11[i],v12[i],v13[i],v14[i],v15[i],v16[i],$
         v17[i],v18[i],v19[i],v20[i],v21[i],v22[i],v23[i],v24[i]]
  
  limit = lim[i]
  div1 = d1[i]
  div2 = d2[i]
  const1 = c1[i]
  const2 = c2[i]
     
  if slit lt limit then begin
     islit = fix(slit/div1) + const1 - 1
  endif else begin
     islit = fix(slit/div2) + const2 - 1
  endelse
  
  if islit ge 0 && islit lt 24 then begin
     slitval = val[islit]
     if slitval eq 0 then begin
        print, 'Warning: slit angle out of range'
        slitval = 0.01
     endif
  endif else begin
     print, 'Warning: slit angle out of range'
     slitval = 0.01
  endelse
  hdr['SLITVAL'] = slitval

  ; set efl0 from telescop if necessary
  tel = hdr['TELESCOP']
  case strupcase(tel) of
     else: efl0 = 10. * 300.
  endcase
  hdr['EFL0'] = efl0

  ; set some tort params from focal reducer / wave number
  ; Note: no focal reducer in fifi_ls
  waveno0 = hdr['WAVENO0']
  fred = 1.0
  hdr['HRFL'] = hdr['HRFL0']/fred
  hdr['XDFL'] = hdr['XDFL0']/fred
  hdr['EFL'] = hdr['EFL0']/fred
  hdr['PLTSCALE'] = hdr['PIXELWD'] / (hdr['EFL']*4.848e-06)
  hdr['SLITWID'] = hdr['SLITVAL'] / (hdr['EFL0']*4.848e-06)
  hdr['OMEGAP'] = hdr['PLTSCALE'] * hdr['SLITWID'] * (4.848e-06)^2
  hdr['HRDGR'] = 0.3*2.54*0.996*cos(hdr['HRG'])

  wno0 = hdr['WNO0']
  if ~valid_num(wno0) || wno0 eq -9999 || wno0 le 0 then begin
     wno0 = waveno0
  endif
  hrdgr = hdr['HRDGR']
  hrr = hdr['HRR']
  hrorder = round(2.0*hrdgr*wno0*sin(atan(hrr)))
  hrr0 = tan(asin(hrorder/(2.0*hrdgr*wno0)))
  if abs(hrr-hrr0) gt 0.01 then begin
     hdr['HRR'] = hrr0
  endif
     
  ; bad pixel mask
  if file_test(datapath+'bpm.fits') then $
     hdr['BPM'] = datapath+'bpm.fits'

  ; linearity file
  if file_test(datapath+'lincoeff.fits') then $
     hdr['LINFILE'] = datapath+'lincoeff.fits'

  ; check for pinhole mode
  ; not used for FIFI_LS yet
  hdr['PINHOLE'] = 0

  ; add the current date/time
  get_date, dt, /TIMETAG
  hdr['DATE'] = dt

  return, hdr
end
