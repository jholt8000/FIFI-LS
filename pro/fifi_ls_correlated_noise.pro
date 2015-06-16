;+
; NAME: fifi_ls_correlated_noise 
;
; PURPOSE: Find correlated noise before fitting ramps
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=correlated_noise(FILENAME, OUTDIR, s2n=s2n,
; DYN_BADPIX=DYN_BADPIX, verbose=verbose) 
;
; INPUTS: FILENAME - Name of the raw data file (including the path if not in
;                    in current working directory) 
;
; OPTIONAL INPUTS:
;
; KEYWORDS: S2N - signal-to-noise value below which data will be 
;                 considered questionable and flagged (wavelength will 
;                 be set to NaN).  If this keyword is not set, a default
;                 S/N = 10 (value subject to change) will be used.  Note that
;                 S2N = 0 means keyword is NOT set, and therefore
;                 S/N = 10 will be used.  To disable S/N filtering,  
;                 set S2N = -1.
;
;           DYN_BADPIX - factor applied to dynamic bad pixel detection.
;                      Dynamic bad pixel detection flags data points that
;                      are larger (either positive or negative) than
;                      f * median value in the same spectral channel. To
;                      disable dynamic bad pixel detection, must explicitly
;                      set DYN_BADPIX = -1.  If this keyword it not set,
;                      f = 2.  To change the value of f, set this keyword
;                      to the value, example DYN_BADPIX = 5
;
; OUTPUTS: RESULT: 1. -1 if failure (or return from displaying help file), 1 if
;                  successful
;                  2. fits file containing linear fit (a+bx) result and sigma
;                     File name: *.chop#.CLEAN.rampfits.fits
;                     fits content:
;                      prime header
;                      binary table - (header + [3,450] data array) x 
;                                     (# of ramps)
; 
; OPTIONAL OUTPUTS:
;
; SIDE EFFECTS: 1. FIFILS fits keyword DATAPROD gets updated to FITTED
;               2. If the 26th "sky" pixel and partial ramps haven't 
;                  already been removed (a procedure performed in 
;                  remove_noise.pro), they will be removed. 
;               3. Data points above an emperically determined cutoff value
;                  (search for the variable "cut" in the code) will be 
;                  not used in line fitting.
;               4. If choplength is longer than ramplength, the ramps on the
;                  same chop will be averaged before line fitting. 
;               5. If ramplength is longer than choplength, the ramps will
;                  be divided into subramps equal to the choplength, and each
;                  subramp will be fitted individually.
;		6. Garbage ramps (as defined as averaged ramps or sub-ramps having
;          less than 3 data points to fit) are assigned NaN for their fit
;          parameters and sigma.  Valid ramps use all valid data points for 
;          fitting (this is a change from "First data point of each ramp is 
;          omitted from line fitting." statement used prior to 2014-10-22.
; 
; RESTRICTIONS: the raw data file must consiste of complete frames
; (full frames are not necessary)
;
; PROCEDURE: 1. Read input file.
;            2. Flag saturated data points; fit a straight line through
;               remaining points on each (sub)ramp; calculate error.
;            3. Create fits file and write results to disk.
;
; EXAMPLE: result=correlated_noise('00006_rawlw021.L1.chop0.CLEAN.fits') 
;
; MODIFICATION HISTORY:
; 4Mar08  k.n.  based on fitrawramps.pro
;               accepts a fits file as input, not the raw binary file (.fif)
;
; 31Mar08 k.n.  modified to remove the 26th sky pixel (chopper info) and the 
;		first frame (partial or full) if they already hadn't been by
;		remove_noise.pro
;
; 17Sep08 k.n.  updated documentation
;
; 23oct08 k.n.  modified to average ramps before fitting if necessary
;
; 16jan09 k.n.  added line to update newly introduced keyword DATAPROD
;
; 20feb09 k.n.  introduced new cutoff for saturated ramps (~1.4e4), omit last
;               data point in each non-saturated ramps from fitting
;                 note: the cutoff and last data point exclusion might change
;                 in the future.
;
; 18jun09 k.n.  square brackets ([]) now allowed in file name
;
; 17may10 k.n.  omit first data point of each ramp from fitting
;		  note: see 20feb09 entry above
;
; 17feb15 J.H.  Fix typo in check_scn call and update mrdfits to use
;               filename instead of file unit to work with newer astrolib.
; 
; 24oct14 k.n.  1. revert parameter "cut" back to 0.8V
;        2. use all valid data points for fitting ramps 
;        3. updated fitting into 2 cases:
;           a. non-saturated ramps and saturated ramps with more than 
;              3 data points remaining:
;                - see item 2 above 
;           b. saturated ramps with 3 or less data points remaining:
;                - don't fit and set fit parameters and sigma to NaN
;        4. deals with missing ramps when co-adding multiple frames
;           (due to removing 1st ramp, or last ramp missing). 
;
; 14nov14 k.n.  1. Added S2N keyword to enable S2N filtering.  S2N 
;                  filtering was previously part of calibrate.pro. 4POINT
;                  data not tested.
;               2. Added DYN_BADPIX keyword to enable dynamic bad pixel
;                  flagging.  Dynamic bad pixel detection was previously
;                  part of calibrate.pro.
;    
; $Id: correlated_noise.pro,v 1.26 2014-03-24 23:55:09 klein Exp $ 
;-
function fifi_ls_correlated_noise, filename, outdir, parent=parent, s2n=s2n, DYN_BADPIX=DYN_BADPIX, verbose=verbose

;--------------
;  help file
;--------------
;if n_params() ne 2 then begin
;  doc_library,'correlated_noise'
;  return, -1
;endif

;-------------------------------------
;  make sure the input file exists
;-------------------------------------


;------------------------
;  define parameters
;------------------------ 
HEADER_WORDS=8
ARRAY_SIZE_X=5
ARRAY_SIZE_Y=5
ARRAY_SIZE_LAM=18
ARRAY_SIZE=(ARRAY_SIZE_X*ARRAY_SIZE_Y+1)*ARRAY_SIZE_LAM
FRAME_SIZE=HEADER_WORDS+ARRAY_SIZE

;---- header words --
P_START=0L
START='8000'X
P_FRM_CNT_L=1L
P_FRM_CNT_H=2L
P_FLAG=3L
P_SMPL_CNT=4L
P_RMP_CNT=5L
P_SCN_IDX=6L
P_SPARE=7L
SPARE=not START

;--------------------------
;  read input data file 
;--------------------------
openr,unit, filename,/get_lun
raw_data=mrdfits(unit,0,primehead,/unsigned)
raw_data=mrdfits(unit,0,exthead,/unsigned)
close,unit & free_lun,unit
scn=scnstrFromFits(primehead)

; to work with new astrolib
raw_data2 = mrdfits(filename, 1, /unsigned)
raw_data = raw_data2

if check_scn(scn,errors) then message,/info,errors
C_SCHEME=scn.ch_scheme               ; chopper scheme
CHOP_LN=scn.choplength               ; choplength
redBlue=((raw_data[0].header)[P_FLAG]  and 2) / 2   ; yields 0 for red and 
                                                    ;        1 for blue
ramplength=scn.ramplength[redBlue]
subramp_length=scn.subramp_length[redBlue]

;; -- define output structure and file name --
fitfilename=strmid(filename,0,strlen(filename)-4)+'rampfits.fits'

ff=0L                           ;fitfile pointer
i=0L                            ;current record (file pointer)

;-----------------------------------------------------------
;  "trim" data if it already hasn't been trimmed (i.e. 26th
;  pixel, partial ramps at the beginning or end of data file)
;-----------------------------------------------------------
s=size(raw_data.data)
if s[1] eq ARRAY_SIZE_X*ARRAY_SIZE_Y+1 then begin	; if the 26th pixel is
        ; there, remove it (chopper info)
 foo=(raw_data.data)[0:ARRAY_SIZE_X*ARRAY_SIZE_Y-1,*,*]
 temp=create_struct('header',uintarr(8),'data',fltarr(25,18))
 data0=replicate(temp,s[3])     ; define smaller data.data array
 data0.header=raw_data.header
 data0.data=foo
 raw_data=temporary(data0)          ; rename it to the original structure name
 foo=0 & temp=0                 ; free up some memory
endif

;RK 2014-03-09
; old code too fancy
;i=0L
;first_ramp=(raw_data[0].header)[P_RMP_CNT]  
;while (raw_data[i]).header[P_SMPL_CNT] ne 0 do i=i+1
;head=i
;anz_records=n_elements(raw_data)
;last_ramp=(raw_data[anz_records-1]).header[P_RMP_CNT]
;last_frame=(raw_data[anz_records-1]).header[P_SMPL_CNT]
;if last_frame eq ramplength-1 then raw_data=raw_data[head:*] $	; remove 1st ramp only 
;else begin
;while (raw_data[i]).header[P_RMP_CNT] lt last_ramp do i=i+1
;  tail=i
;  raw_data=raw_data[head:tail]		; remove 1st and last partial ramps
;endelse
;just drop first ramp after long syncs, i.e. ramp_count == 0
head = 0 
if raw_data[0].header[P_RMP_CNT] eq 0 then begin
   if raw_data[ramplength].header[P_SMPL_CNT] eq 0 then begin
      head=ramplength
   endif else begin
      i=1L
      while (raw_data[i]).header[P_SMPL_CNT] ne 0 do i=i+1
      head = i
      print,"correlated_noise: data didn't start at sampl count 0 but at ",raw_data[ramplength].header[P_SMPL_CNT]
      print,"correlated_noise: ramp starts with ramp count",raw_data[head].header[P_RMP_CNT]
   endelse
   raw_data=temporary(raw_data[head:*])
endif
;
; to be done later if needed: 
; if last ramp is partial, discard last ramp 
;

anz_records=n_elements(raw_data)        ; number of readouts
nramps=anz_records/subramp_length           ; number of ramps

;-------------------------------
;  define effective choplength
;-------------------------------
case strtrim(C_SCHEME,2) of
  '2POINT': C=CHOP_LN*1.
  '4POINT': C=CHOP_LN/2.
      else: begin
            print, 'Invalid chopper scheme. Exiting correlated_noise.pro.'
            return, -1
            end
endcase

;------------
;  fitting
;------------
Volts2Counts=65536/3.63d0
cut=0.8*Volts2Counts	; cut = 14443.196

if C gt ramplength then begin	; if choplength > ramplength average ramps 
  n=C/ramplength		; before fitting

  message, /info, 'Ramps will be averaged before fitting.'
  if verbose then begin
      print, 'C / ramplength = ', n
      print, 'C=',C,' ramplength=',ramplength
      print, 'anz_records=',anz_records
      print, 'subramp_length=',subramp_length
      print, 'nramps=',nramps
   endif
  
  fullframe=fltarr(ramplength,ARRAY_SIZE_X*ARRAY_SIZE_Y,ARRAY_SIZE_LAM)
  fit=fltarr(3,ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM)
  x=findgen(ramplength)

  foo=create_struct('header',uintarr(8), $
      'data',fltarr(3,ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM))
  ;if head gt 0 then nramps += 1     ; if the first frame was dropped
  if nramps mod 2 ne 0 then nramps += 1  ; if there are odd #s of ramps
  fits=replicate(foo,nramps/n)		; anz_records/C = nramps/n

  i=0L
  while i+ramplength le anz_records do begin
    ;print, i, ff
    ;-- check if we are at a start of a full frame
    if (raw_data[i]).header[P_SMPL_CNT] eq 0 then begin
      nr=0L & p=0L
      ;while p lt n do begin
      while p lt n and i lt anz_records do begin
        header=raw_data[i].header
        rampn=header[P_RMP_CNT]
        if (rampn mod n) eq p then begin
          ;-- stack frames 
          for r=0,ramplength-1 do begin
             fullframe[r,*,*]=fullframe[r,*,*]+raw_data[i].data

            i=i+1
          endfor
          r=0
          nr=nr+1
        endif
        p=p+1
      endwhile
      fullframe=fullframe/nr	; average ramps
      fullframe=reform(fullframe, ramplength, ARRAY_SIZE_X* $
                       ARRAY_SIZE_Y*ARRAY_SIZE_LAM, /overwrite)
    ;-- linear fit of ramps 
      for pix=0, ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM-1 do begin
        ; crude attempt at removing "bent" (saturated) ramps
        q=where(fullframe[*,pix] gt cut, nx)
        if nx ge 0 and nx lt subramp_length-3 then begin
          if nx eq 0 then begin  ; non-saturated ramps
            q_start = 0
            q_end = subramp_length-1    ; should this be ramplength ?
            ; use the full ramp because detector noise is MUCH smaller
            ; than other sources of noises IN FLIGHT.  For lab data,
            ; you may want to discard the 1st and last points by setting
            ; q_start = 1 and q_end = subramp_length-2.   2014-10-22
          endif else begin   ; saturated ramps with more than 3 points to fit
            q_start = 0
            q_end = min(q)-1
            ; change q_start (and q_end) if using Lab data (maybe).
          endelse

          fit[0:1,pix]=linfit(x[q_start:q_end],fullframe[q_start:q_end, pix], $
                        /double,chisq=chisq,sigma=sigma)  ;y=A+Bx
          fit[2,pix]=sigma[1] ;/sqrt(chisq/(ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM-2)) ; sigma for slope
          line_fit = fit[1,pix]*x[q_start:q_end] + fit[0,pix]

          ;plot,line_fit, title='pix='+string(pix)
          ;oplot, fullframe[q_start:q_end,pix],psym=2
          ;wait,0.02
          
        endif else fit[*, pix] = !VALUES.F_NAN    ; garbage ramps
        ;PRINT, fit[*, pix]
      endfor
      fits[ff].header=header
      fits[ff].data=fit
      ff=ff+1
      fullframe[*]=0. & fit[*]=0.
    endif else begin
      message,/info,"Full frame lost, searching..."
      while (raw_data[i].header)[P_SMPL_CNT] ne 0 and i lt anz_records-1 do i=i+1
      if (raw_data[i].header)[P_SMPL_CNT] ne 0 then message,/info,"...found at record "+string(i)
    endelse


 endwhile
endif else begin		; other cases - divide ramps into subramps 
;  if ramplength eq C then n=1
;  if ramplength eq 2*C and strtrim(C_SCHEME,2) eq '4POINT' then n=1
;  if ramplength gt C and strtrim(C_SCHEME,2) eq '2POINT' then $
;     n=ramplength/(2*C) else n=1
  fullframe=fltarr(subramp_length,ARRAY_SIZE_X*ARRAY_SIZE_Y,ARRAY_SIZE_LAM)
  fit=fltarr(3,ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM)
  x=findgen(subramp_length-1)	; x=findgen(C-1)

  message, /info, 'Ramps will be sub-divided before fitting.'
if verbose then  print, 'ramplength / C = ', n

  foo=create_struct('header',uintarr(8), $
      'data',fltarr(3,ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM))
  fits=replicate(foo,nramps)

  i=0L
  if (raw_data[i]).header[P_SMPL_CNT] eq 0 then begin
    while i+subramp_length le anz_records do begin
   ;  print, i, ff
      ;-- check if we are at a start of a subramp
      if ((raw_data[i]).header[P_SMPL_CNT] mod subramp_length eq 0) then begin 
        header=raw_data[i].header			; save header
        for r=0,subramp_length-1 do begin			; copy frames
          fullframe[r,*,*]=fullframe[r,*,*]+raw_data[i].data
          i=i+1
        endfor
        fullframe=reform(fullframe,C,ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM,/overwrite)
      ;-- linear fit of ramps 
        for pix=0,ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM-1 do begin
          ; a crude attempt at removing "bent" ramps
          q=where(fullframe[*,pix] gt cut, nx)
          if nx ge 0 and nx lt subramp_length-3 then begin
            if nx eq 0 then begin  ; non-saturated ramps
              q_start = 0
              q_end = subramp_length-1
              ; use the full ramp because detector noise is MUCH smaller
              ; than other sources of noises IN FLIGHT.  For lab data, 
              ; you may want to discard the 1st and last points by setting
              ; q_start = 1 and q_end = subramp_length-2.   2014-10-22 
            endif else begin    ; saturated ramps with more than 3 points to fit
              q_start = 0
              q_end = min(q)-1
              ; change q_start (and q_end) if using Lab data (maybe).
            endelse
            fit[0:1,pix]=linfit(x[q_start:q_end],fullframe[q_start:q_end, pix], $
                       /double,chisq=chisq,sigma=sigma)  ;y=A+Bx
            fit[2,pix]=sigma[1] ;/sqrt(chisq/(ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM-2))	; sigma for slope
          endif else fit[*, pix] = !VALUES.F_NAN   ; garbage ramps
        endfor

        fits[ff].header=header
        fits[ff].data=fit
        ff=ff+1
        fullframe[*]=0. & fit[*]=0.
      endif else begin
        message,/info,"Full frame lost, searching..."
        while (raw_data[i].header)[P_SMPL_CNT] ne 0 and i lt anz_records-1 do i=i+1
        if (raw_data[i].header)[P_SMPL_CNT] ne 0 then message,/info,"...found at record "+string(i)
      endelse
    endwhile
  endif
endelse

;------------------------
; S/N filtering   nov2014
;------------------------
if ~keyword_set(s2n) then s2n=10.   ; default subject to change
s2n= float(s2n)
if s2n gt 0 then begin
  message, /info, 'S/N filter applied with S/N =' + string(s2n)
  n_ramps = n_elements(fits)
  for i = 0, n_ramps - 1 do begin
    fit = fits[i].data
    for pix = 0, ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM - 1 do $;begin
      if fit[1,pix] / fit[2,pix] lt s2n then $ 
       fit[1,pix] = !VALUES.F_NAN    ; set flux to NaN
    ;endfor
    fits[i].data = fit
  endfor
endif else message, /info, 'S/N filter NOT applied.'

;--------------------------
; median filtering  nov2014
;--------------------------
; implementation of Rainer's dynamic bad pixel detection    oct2014

if ~keyword_set(dyn_badpix) then f = 2 else f = dyn_badpix

if f gt 0 then begin
  message, /info, 'Performing median filtering using f=' + string(f)
;  n_ramps = n_elements(fits)
;  ; med = median((fits[*].data)[1,*])     ; median over the entire file  
;  for i = 0, n_ramps - 1 do begin        ; loop over each frame 
;    flux = fits[i].data[1,*]   ; flux: 1 x 450
;    flux_r= reform(flux, 25, 18)  ; 25 x 18
;    for j = 0, 24 do begin       ; loop over each spatial pixel
;      tempf = flux_r[j, 1:16]    ; pix 1-16 are the "real" channels
;      med = median(tempf, /even)
;      q = where(tempf gt f * med or tempf lt med / f, nx)
;      ;PRINT, med / f, med, f * med, nx, q
;      if nx gt 0 then begin
;        tempf[q] = !VALUES.F_NAN    ; set flux to NaN
;        flux_r[j, 1:16] = tempf
;      endif
;    endfor
;    something = reform(flux_r, 1, 450)
;    fits[i].data[0,*] = something
;  endfor
  ; calculate mean slope for all pixels
  pix_mean = fltarr(ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM)
  for i = 0, ARRAY_SIZE_X*ARRAY_SIZE_Y*ARRAY_SIZE_LAM - 1 do $ 
    pix_mean[i] = mean(fits.data[1, i, *])
  ; reshape pix_mean into 25 x 18
  flux_r = reform(pix_mean, 25, 18)
  for j = 0, 24 do begin       ; loop over each spatial pixel
    tempf = flux_r[j, 1:16]    ; pix 1-16 are the "real" channels
    med = median(tempf, /even)
    q = where(tempf gt f * med or tempf lt med / f, nx)
    if nx gt 0 then begin
      tempf[q] = !VALUES.F_NAN    ; set flux to NaN
      flux_r[j, 1:16] = tempf
    endif
  endfor
  
endif else message, /info, 'Skipping dynamic bad pixel detection.'

;--------------------------
;  write results to disk 
;--------------------------
;---- update prime header 
;version='$Id: correlated_noise.pro,v 1.26 2014-03-24 23:55:09 klein Exp $'
fifisoftware=getenv('FIFISOFTWARE') 
version = gitdescribe(fifisoftware)
sxaddhist,'Processed by correlated_noise.pro '+version,primehead

;--- is there already a file by the same name?  if so, remove it.
foo=file_search(fitfilename)
if foo[0] ne '' then file_delete, fitfilename,/quiet

; fluxer keywords
sxaddpar,primehead,'CRPIX1',2.5
sxaddpar,primehead,'CRPIX2',2.5
sxaddpar,primehead,'CRPIX3',0
sxaddpar,primehead,'CRVAL1',sxpar(primehead,'OBSLAM')
sxaddpar,primehead,'CRVAL2',sxpar(primehead,'OBSBET')
sxaddpar,primehead,'CRVAL3',0
sxaddpar,primehead,'CDELT1',6./3600. ; deg/pix
sxaddpar,primehead,'CDELT2',6./3600. ;
sxaddpar,primehead,'CDELT3',1        ;constant delta lambda in this case

fxaddpar, primehead, 'EQUINOX',2000   ; right now set to FK5
fxaddpar, primehead, 'DATAPROD','FITTED'	; update DATAPROD keyword
fxwrite, fitfilename,primehead,/noupdate
mwrfits,fits, fitfilename

;write output file in fluxer readable format
data=fits.data
data=reform(data[1,*],25,18)

fluxerfilename=strmid(filename,0,strlen(filename)-4)+'rampfits_25by18.fits'
fxwrite,fluxerfilename,primehead,/noupdate
mwrfits,data,fluxerfilename

;fluxer,fluxerfilename
close,/all

return, fits
end
