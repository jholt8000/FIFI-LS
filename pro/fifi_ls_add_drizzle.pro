;+
; NAME: fifi_ls_add_drizzle 
;
; PURPOSE: Combine multiple drizzled data into one final data cube. 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE:
; RESULT=fifi_ls_add_drizzle(FILELIST,FILEROOT,CONFIGFILENAME=CONFIGFILENAME,
; verbose=verbose)
;
; INPUTS: FILELIST - name of the file containing the names (including 
;                    directory path) of drizzled data to be combined
;         FILEROOT - directory path and root name for output file.
;
; KEYWORDS: CONFIGFILENAME - name of the file defining drizzle parameters,
;                            including the directory path to it; if not used,
;                            a fallback file will be accessed. 
;
; OPTIONAL INPUTS: 
;
; OUTPUTS: 1. RESULT - 1 if successful, -1 if failure
;          2.  FITS file containing a combined drizzled data 
;              [flux,weight,sigma] as a primary array, and the data header 
;              of the first data readout of each input drizzled data file
;              as an extension binary table.
;              File name: FILEROOT.cube.FITS
;              File content:
;              prime header and prime array 
;                [# of x pixels,# of y pixels,# of lambda pixels,3]
;              binary table = (data header of the first data readout) x
;                             (# of combined files) 
;           
; OPTIONAL OUTPUTS: 
;
; SIDE EFFECTS: 1. Uses a lot of DO-loops in the code to save memory.
;               2. Astrometric FITS keywords may be updated.
; 
; RESTRICTIONS: 1. All data files are assumed to be of the same
;                  spectral channel (Red or Blue). 
;
; PROCEDURE: 1. Inspect drizzle configuration file (CONFIGFILENAME or 
;               fall-back file) to determine method of output cube creation.
;            2. Create output cube.
;            3. Combine the files specified in FILELIST.
;            4. Create FITS file and write results to disk.
;
; EXAMPLE: result=fifi_ls_add_drizzle(filelist,fileroot)
;
; MODIFICATION HISTORY:
; 26aug09 k.n.  original version 
;
; 02sep09 k.n.  bug fixes
;
; 10sep09 k.n.  made CONFIGFILENAME a keyword (formally a required input)
;
; 18sep09 k.n.  corrected an error in calculating errors
;
; 12oct09 k.n.  exit if list of input files (FILELIST) is empty
;
; 04jan10 k.n.  bug fix: normalize added weighted flux by total weight
;
; 11may10 k.n.  bug fix: initialize temporary array used to paste input
;               cube in to for the case where grid boundaries are not
;               specified in the drizzle output grid config file.
; 
; 13may10 k.n.  all "for i=0" statements are now "for i=0L"
;
; 24feb14 k.n.  added "return, 1" at end of program (successful execution of
;               program). 
; 
; $Id: fifi_ls_add_drizzle.pro,v 1.20 2014-03-24 23:48:06 klein Exp $
;-
function fifi_ls_add_drizzle, filelist, fileroot, configfilename=configfilename, verbose=verbose

  for i=0,n_elements(filelist)-1 do begin
     
;-------------
;  help file 
;-------------
if n_params() lt 2 then begin
  doc_library,'fifi_ls_add_drizzle'
  return, -1
endif

;----------------
; sanity checks
;----------------
; check to see if the file list exists
if file_test(filelist[i]) eq 0 then begin
  print, 'fifi_ls_add_drizzle.pro: List of input files does not exist.  Exiting.'
  return, -1
endif
;if keyword_set(configfilename) and file_test(configfilename) eq 0 then begin 
;  print, 'fifi_ls_add_drizzle.pro: Drizzle config file does not exist.  Exiting.'
;  return, -1
;endif
; check to see if the input files all exist
if (file_info(filelist[i])).size eq 0 then begin
  print, 'fifi_ls_add_drizzle.pro: List of input files is empty.  Exiting.'
  return, -1
endif

filename=filelist

;readcol,filelist[i],filename,format='a',delim=' '
for i=0L,n_elements(filelist)-1 do begin
  if file_test(filename[i]) eq 0 then begin
    print, 'fifi_ls_add_drizzle.pro: At least one of the input data files does not exist. Exiting.'
    return, -1
  endif
endfor
if verbose then print,n_elements(filename)
  ; read first input cube
  openr,unit,filename[0],/get_lun
  data=mrdfits(unit,0,primehead,/unsigned)	; drizzled data
  header=mrdfits(unit,0,exthead,/unsigned)	; data header
  
  data2 = mrdfits(filename[0],0,/unsigned)
  data = data2
  if verbose then begin
     print,filename[0]
     print,size(data)
  endif
  
  close,unit & free_lun,unit

;-----------------------------------------------------------
; read the drizzle config file: user-specified or fallback
;-----------------------------------------------------------
if strmatch(filename[0], '*lw*') eq 1 then begin
      lambda0 = fxpar(primehead,'G_WAVE_R')
      dx = 12
      dy = 12
endif else begin
      lambda0 = fxpar(primehead,'G_WAVE_B')
      dx = 6
      dy = 6
endelse
   
if strtrim(fxpar(primehead,'OBJ_NAME'),2) ne 'TELSIM' OR strtrim(fxpar(primehead,'OBJECT'),2) ne 'TELSIM' then begin
   
    x0=fxpar(primehead,'OBSLAM')
    y0=fxpar(primehead,'OBSBET')  
;    rot = fxpar(primehead,'DET_ANGL')  
;    dlambda = 0.1  ; this is the thing that needs a change! let's make this an input param and make it smaller by default -JH

    ;these arent used -jh
    ;coordsys0=fxpar(primehead,'COORDSYS')  
    ;crdsysmp0=fxpar(primehead,'CRDSYSMP')  
    ;primaray0=fxpar(primehead,'PRIMARAY')

    ; I cannot find anywhere that lmin is being set in the config file -JH
    lmin = !values.f_nan
    lmax = !values.f_nan
     
  ; I cannot find xmin and xmax
  ; either in pipeline.script or config file
    xmin = !values.f_nan
    xmax = !values.f_nan
    ymin = !values.f_nan
    ymax = !values.f_nan
    
endif

;-------------------------------------------------------------------------
; if lmin, xmin, ymin (and max's too) are defined, all input
; data cubes can be added up
;-------------------------------------------------------------------------

; combine the input files: open one file at a time so 
; we won't have memory shortage

if finite(lmin) and finite(xmin) and finite(ymin) and $
   finite(lmax) and finite(xmax) and finite(ymax) then begin
  print,'fifi_ls_add_drizzle.pro: output cube has edges defined; input cubes will be added.'

  ; check if crpix1 is a long integer (i.e. new drizzle.pro has been used
  ; grid boundaries are x0 + n * dx and so on - it not, exit
;RK 2014-03-05 Now CRPIX is double, disabled checking
;  if size(fxpar(primehead,'CRPIX1'),/type) ne 3 then begin
;    print, 'fifi_ls_add_drizzle.pro: rerun drizzle.pro version 03may10 or newer. Exiting.'
;    return, -1
;  endif

  s0=size(data)
  if s0(0) eq 4 then begin
     sum=data                   ; array[#lambda,#x,#y,3]
    ; data[*,*,*,0] is flux for all x,y,lambda(0.01)
    ; 
    sum[*,*,*,0]=double(data[*,*,*,0]*data[*,*,*,1])	; weighted average
    sum[*,*,*,2]=double(data[*,*,*,2]^2*data[*,*,*,1]^2)	; sigma^2
    head=header	; data header array 
  endif else begin
    print, 'fifi_ls_add_drizzle.pro: data not 4 dimensional.  Skipping file ' $
            +strtrim(filename[0],2)
    sum=[!values.f_nan]	; set output cube to NaN
  endelse

  ; add the rest 
  for i=1L,n_elements(filename)-1 do begin
    openr,unit, filename[i],/get_lun
    data=mrdfits(unit,0,exthead,/unsigned)	; drizzled data
    header=mrdfits(unit,0,exthead,/unsigned)    ; data header
    data=mrdfits(filename[i], 1, /unsigned)
    close,unit & free_lun,unit

    s=size(data)
    if total(s-s0) eq 0. then begin
       ; why isn't this done in drizzle? -jh
      if n_elements(sum) eq 1 then begin	; if previous output cube is 
        sum=data    ; array[#x,#y,#lambda,3]	; NaN 
        sum[*,*,*,0]=double(data[*,*,*,0]*data[*,*,*,1])      ; weighted average
        sum[*,*,*,2]=double(data[*,*,*,2]^2*data[*,*,*,1]^2)  ; sigma^2
        head=header ; data header array 
      endif else begin	; "normal" case 
        sum[*,*,*,0]=double(temporary(sum[*,*,*,0])+data[*,*,*,0]*data[*,*,*,1])
		; weighted average
        sum[*,*,*,1]=double(temporary(sum[*,*,*,1])+data[*,*,*,1])	
		; accumulated weight
        sum[*,*,*,2]=double(temporary(sum[*,*,*,2])+data[*,*,*,2]^2*data[*,*,*,1]^2)			; sigma^2
        head=[[head],[header]]	; data header array
      endelse
    endif else $ 
      print, 'fifi_ls_add_drizzle.pro: data not same dimention as first data cube.  Skipping file '+strtrim(filename[i],2)
  endfor

  sum[*,*,*,2]=temporary(sqrt(sum[*,*,*,2]))            ; sigma

  q=where(sum[*, *, *, 1] gt 0)  ; select non-zero weights
  if q[0] ne -1 then begin
    q=array_indices(sum,q)
    for i=0L,n_elements(q)/4-1 do begin
      sum[q[0,i],q[1,i],q[2,i],0]=sum[q[0,i],q[1,i],q[2,i],0]/sum[q[0,i],q[1,i],q[2,i],1]  ; flux
      sum[q[0,i],q[1,i],q[2,i],2]=sum[q[0,i],q[1,i],q[2,i],2]/sum[q[0,i],q[1,i],q[2,i],1]  ; sigma 
    endfor
  endif else $
    print,'fifi_ls_add_drizzle.pro: all weights were zero... something is very wrong here.'

endif else begin
  
;------------------------------------------------------------------------ 
; if lmin etc are not defined, output cube must be created first.
;------------------------------------------------------------------------ 

  ; check if crpix1 is a long integer (i.e. new drizzle.pro has been used
  ; grid boundaries are x0 + n * dx and so on - it not, exit
;RK 2014-03-05 Now CRPIX is double, disabled checking
;  if size(fxpar(primehead,'CRPIX1'),/type) ne 3 then begin
;    print, 'fifi_ls_add_drizzle.pro: rerun drizzle.pro version 03may10 or newer on ' $
;          +filename[0]+'. Exiting.'
;    return, -1
;  endif

  s=size(data)  ; s(0)=number of axes, s(1)=# of X elements, s(2)=# of Y...
if verbose then  print,'s=',s
  if s(0) eq 4 then begin
; check [x0,y0,lambda0] against [x0,y0,lambda0] - should be the same
    crval1=fxpar(primehead,'CRVAL1') & crpix1=fxpar(primehead,'CRPIX1')
    crval2=fxpar(primehead,'CRVAL2') & crpix2=fxpar(primehead,'CRPIX2')
    crval3=fxpar(primehead,'CRVAL3') & crpix3=fxpar(primehead,'CRPIX3')
    naxis1=fxpar(primehead,'NAXIS1') & naxis2=fxpar(primehead,'NAXIS2')
    naxis3=fxpar(primehead,'NAXIS3')

    xdiff=abs(crval1-double(x0)) & ydiff=abs(crval2-double(y0))
    ldiff=abs(crval3-double(lambda0)) & gap=1d-4

    if (xdiff lt gap) and (ydiff lt gap) and (ldiff lt gap) then begin 
      xmin=1-crpix1 & xmax=xmin+naxis1
      xminmax=[xmin,xmax]
      ymin=1-crpix2 & ymax=ymin+naxis2
      yminmax=[ymin,ymax]
      lmin=1-crpix3 & lmax=lmin+naxis3
      lminmax=[lmin,lmax]
    endif else begin
      print,'fifi_ls_add_drizzle.pro: user-specified reference value(s) do not match actual value(s).  Skipping file '+strtrim(filename[0],2)
      xminmax=[-!values.f_nan,+!values.f_nan]
      yminmax=[-!values.f_nan,+!values.f_nan]
      lminmax=[-!values.f_nan,+!values.f_nan]
    endelse 
    good=[0]	; keep track of which files are valid
  endif else begin
    print,'fifi_ls_add_drizzle.pro: data not 4-dimentional. Skipping file.'
    good=[!values.f_nan]
  endelse

  ; read the rest
  for i=1L,n_elements(filename)-1 do begin
    openr,unit, filename[i],/get_lun
    data=mrdfits(unit,0,primehead,/unsigned)
    header=mrdfits(unit,0,exthead,/unsigned) ; data header
    data=mrdfits(filename[i],0,/unsigned)
    close,unit & free_lun,unit

    s=size(data)  ; s(0)=number of axes, s(1)=# of X elements, s(2)=# of Y...
    if s(0) eq 4 then begin 
      crval1=fxpar(primehead,'CRVAL1') & crpix1=fxpar(primehead,'CRPIX1')
      crval2=fxpar(primehead,'CRVAL2') & crpix2=fxpar(primehead,'CRPIX2')
      crval3=fxpar(primehead,'CRVAL3') & crpix3=fxpar(primehead,'CRPIX3')
      naxis1=fxpar(primehead,'NAXIS1') & naxis2=fxpar(primehead,'NAXIS2')
      naxis3=fxpar(primehead,'NAXIS3')

      xdiff=abs(crval1-double(x0)) & ydiff=abs(crval2-double(y0))
      ldiff=abs(crval3-double(lambda0)) & gap=1d-4

      if (xdiff lt gap) and (ydiff lt gap) and (ldiff lt gap) then begin 
        xmin=1-crpix1 & xmax=xmin+naxis1
        xtemp=[xmin,xmax]
        ymin=1-crpix2 & ymax=ymin+naxis2
        ytemp=[ymin,ymax]
        lmin=1-crpix3 & lmax=lmin+naxis3
        ltemp=[lmin,lmax]

        xminmax=[xminmax,xtemp]
        yminmax=[yminmax,ytemp]
        lminmax=[lminmax,ltemp] 

        xminmax=minmax([xminmax[where(finite(xminmax))]]) ; update minmax values
        yminmax=minmax([yminmax[where(finite(yminmax))]])
        lminmax=minmax([lminmax[where(finite(lminmax))]])

        good=[good,i]
      endif else $ 
        print,'fifi_ls_add_drizzle.pro: user-specified reference value(s) do not match actual value(s).  Skipping file '+strtrim(filename[i],2)
    endif else print,'fifi_ls_add_drizzle.pro: data not 4 dimentional.  Skipping file'$
               +strtrim(filename[i],2)+'.'

  endfor

  good=good[where(finite(good))]	; get rid of the NaN in case it's there 

  ; define the dimentions of the output grid
  mm=(xminmax[1]-xminmax[0]) + 1
  nn=(yminmax[1]-yminmax[0]) + 1
  ll=(lminmax[1]-lminmax[0]) + 1

  if finite(mm*nn*ll) and (mm*nn*ll gt 0) then $ 
   sum=fltarr(mm,nn,ll,3) $	; [x,y,l,[flux,weight,sigma]] 
   else begin 
     print, 'fifi_ls_add_drizzle.pro: invalid output grid dimensions.  Exiting.'
     return, -1
   endelse

  temp=sum
  s0=size(sum)

; now add each cube to the output cube
  for i=0L,n_elements(filename[good])-1 do begin
    openr,unit,filename[good[i]],/get_lun
    data=mrdfits(unit,0,primehead,/unsigned)	  ; data 
    header=mrdfits(unit,0,exthead,/unsigned)      ; data header
    data = mrdfits(filename[good[i]],0,/unsigned)
    close,unit & free_lun,unit

    crval1=fxpar(primehead,'CRVAL1') & crpix1=fxpar(primehead,'CRPIX1')
    crval2=fxpar(primehead,'CRVAL2') & crpix2=fxpar(primehead,'CRPIX2')
    crval3=fxpar(primehead,'CRVAL3') & crpix3=fxpar(primehead,'CRPIX3')
    naxis1=fxpar(primehead,'NAXIS1') & naxis2=fxpar(primehead,'NAXIS2')
    naxis3=fxpar(primehead,'NAXIS3')

    xmin=1-crpix1 & xmax=xmin+naxis1
    ymin=1-crpix2 & ymax=xmin+naxis2
    lmin=1-crpix3 & lmax=xmin+naxis3

    m=xmin-xminmax[0]	; indices of cube corner
    n=ymin-yminmax[0]
    l=lmin-lminmax[0]

    temp[*]=0.
    temp[m,n,l,0]=temporary(data[*,*,*,*])
    sum[*,*,*,0]=double(temporary(sum[*,*,*,0])+temp[*,*,*,0]*temp[*,*,*,1]) 
	; weighted average
    sum[*,*,*,1]=double(temporary(sum[*,*,*,1])+temp[*,*,*,1])  
	; accumulated weight
    sum[*,*,*,2]=double(temporary(sum[*,*,*,2])+temp[*,*,*,2]^2*temp[*,*,*,1]^2)        ; sigma^2
    if i eq 0 then head=header else head=[[head],[header]]  ; data header array
  endfor

;; normalize flux and sigma 
  sum[*,*,*,2]=temporary(sqrt(sum[*,*,*,2]))            ; sigma

  q=where(sum[*, *, *, 1] gt 0)  ; select non-zero weights
  if q[0] ne -1 then begin
    q=array_indices(sum,q)
    for i=0L,n_elements(q)/4-1 do begin
      sum[q[0,i],q[1,i],q[2,i],0]=sum[q[0,i],q[1,i],q[2,i],0]/sum[q[0,i],q[1,i],q[2,i],1]  ; flux
      sum[q[0,i],q[1,i],q[2,i],2]=sum[q[0,i],q[1,i],q[2,i],2]/sum[q[0,i],q[1,i],q[2,i],1]  ; sigma 
    endfor
  endif else $
    print,'drizzle.pro: all weights were zero... something is very wrong here.'


; update CRPIX1,CRPIX2,CRPIX3 (pix1,pix2,pix3) 
  pix1=1-xminmax[0] ; the 'first' pixel in FITS convention is 1
  pix2=1-yminmax[0] ; 
  pix3=1-lminmax[0] ; is this correct now? RK 2012-06-28
  fxaddpar,primehead,'CRPIX1',pix1,'Reference pixel number for X'
  fxaddpar,primehead,'CRPIX2',pix2,'Reference pixel number for Y'
  fxaddpar,primehead,'CRPIX3',pix3,'Reference pixel number for LAMBDA'

endelse

;----------------------------------------
; update FITS header, write output file
;----------------------------------------
;fifisoftware=getenv('FIFISOFTWARE') 
;version = gitdescribe(fifisoftware)


;sxaddhist,'Processed by fifi_ls_add_drizzle.pro '+version,primehead
;sxaddhist,'fifi_ls_add_drizzle file list: '+filelist[i],primehead

fxhmake,primehead,sum,/extend  ; this will update the header to change it 
        ; to a prime-array only FITS header without erasing previous entries. 
        ; /extend adds extension tables 
outfilename=fileroot+'.cube.fits'
fxwrite, outfilename,primehead,sum,/noupdate
mwrfits,head, outfilename	; save data header as
        ; extention table - adds file size but seems the easiest way to carry
        ; over information 




return, sum

endfor
  
end
