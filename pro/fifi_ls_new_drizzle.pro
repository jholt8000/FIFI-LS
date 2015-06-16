;+
; NAME: drizzle 
;
; PURPOSE: drizzle FIFI-LS raw data 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=drizzle(FILEROOT,CONFIGFILENAME=CONFIGFILENAME,verbose=verbose)
;
; INPUTS: FILEROOT - root of the file to be drizzled, including the path to it
;
; OPTIONAL INPUTS:
;
; OUTPUTS: 1. result - 1 if successful, -1 if failure
;          2. FITS file containing drizzled data [flux, weight, sigma] as a
;             a primary array, and the data header of the first data readout
;             an extension binary table.
;             File name: *.drizzle.L3.FITS
;             File content: 
;              prime header and prime array 
;                [# of x pixels,# of y pixels,# of lambda pixels,3]
;              binary table = data header of the first data readout
;
; OPTIONAL OUTPUTS: 
;
; KEYWORDS: CONFIGFILENAME - name of the file defining drizzle parameters,
;                            including the path to it.  If unused, a fallback
;                            file will be used.
;
; SIDE EFFECTS: 1. Adds astrometric keywords/values to FITS header.
;               2. SOFIA FITS keyword PROCSTAT gets updated to LEVEL_3.
;               3. Output FITS file now becomes prime array + binary table
;
; RESTRICTIONS: 1. x,y values are absolute mm when FITS keyword OBJ_NAME is
;                  TELSIM, but x,y offset in mm for other cases. 
;
; PROCEDURE: 1. Read input file.
;            2. Define output grid (4-D cube) based on values in either
;               CONFIGFILENAME or the fallback file.
;            3. Drizzle data.
;            4. Create FITS file and write results to disk.
;
; EXAMPLE: result=drizzle('00121_221031_FullScan-G_STRT_R-752500_lw','config.txt')
;
; MODIFICATION HISTORY:
; 21apr09  k.n.  started coding
;
; 26jun09  k.n.  square brackets ([]) in file name now allowed
;
; 02jul09  k.n.  opps, wrote input data as output data.  Corrected.
;
; 16jul09  k.n.  drizzled data array saved as primary array, data headers are
;                saved as a binary table in an extension.
;
; 21jul09  k.n.  changed output to contain only flux, weight, and sigma
;
; 27aug09  k.n.  1. changed output filename suffix to L3.drizzle.FITS from
;                drizzle.FITS
;                2. added line to update SOFIA keyword PROCSTAT to indicate that
;                output file is Level 3 
;
; 01sep09  k.n.  made CONFIGFILENAME a keyword, not a required input; if
;                CONFIGFILENAME is not used, a fallback drizzle configuration
;                file is used.
;
; 12sep09  k.n.  modified for cases with multiple grating positions
;
; 18sep09  k.n.  ignore any flux falling outside the output grid
;
; 09oct09  k.n.  exit if all wavelengths are NaN (i.e. flagged as bad data)
;
; 12oct09  k.n.  changed output filename suffix to drizzle.L3.FITS (to be
;                consistent with L1 and L2 files)
;
; 10dec09  k.n.  corrected error in calculating standard deviation
;
; 03may10  k.n.  change/corrected way the output grid boundaries are
;                calculated; grid boundaries are x0 + n * dx (and so on)
;
; 14may10  k.n.  changed "for i=0" statements to "for i=0L" (long integer)
; 
; 23oct14  k.n.  if lmin is not supplied in drizzle config file and needs
;                to be determined from input data, use NaN keyword with 
;                funcion min so if one grating position has all NaN as
;                its wavelength, it will be ignored but other position
;                will be processed.
;  
;-
function fifi_ls_new_drizzle, infile

  device, decomposed=0
  loadct, 13
  
  fits_info, infile, n_ext=n_ext, /silent
  ncol=255/(n_ext)
  colors = ncol*indgen(n_ext+1)+ncol

  primehead = headfits(infile)
  channel = strtrim(fxpar(primehead, 'CHANNEL'),2)

  if (channel eq 'RED') then begin
     dx = 12
     dy = 12
  endif else begin
     dx = 6
     dy = 6
  endelse

  psym_arr = indgen(7) + 1
  psym_arr[2] = 6
  psyms=[psym_arr, psym_arr, psym_arr, psym_arr, psym_arr, psym_arr]
  k=1
 ; set_plot,'PS'
 ; device,filename='test1.ps',/color, bits=8
  while (k lt n_ext+1) do begin
     
     fits_struct = mrdfits(infile, k, /silent, /unsigned)
     data = [fits_struct.data]
     lambda = fits_struct.lambda
     stddev = [fits_struct.stddev]
     xs = reform(fits_struct.xs, 5, 5, 16)
     ys = reform(fits_struct.ys, 5, 5, 16)

     ;exthdr = fits_struct.exthdr
                                ;ind_pos = fxpar(exthdr,'INDPOS')
     
     i=0
     while (i lt 5) do begin
        j=0
        while (j lt 5) do begin
           if (k eq 1) and (i eq 0) and (j eq 0) then  begin
              lambda_pix = lambda[i,j,*]            
              flux_pix = data[i,j,*]
              xs_pix = xs[i,j,*]
              ys_pix = ys[i,j,*]
              
           endif else begin
              lambda_pix = [lambda_pix, lambda[i,j,*]]
              flux_pix = [flux_pix, data[i,j,*]]
              xs_pix = [xs_pix, xs[i,j,*]]
              ys_pix = [ys_pix, ys[i,j,*]]


              if (i eq 2) and (j eq 2) then begin
                 if k eq 1 then begin
                  plot, lambda[2,2,*], data[2,2,*], psym=1, $
                       xrange=[min(lambda[2,2,*]), max(lambda[2,2,*])],$
                        yrange=[min(data[2,2,*])-1, max(data[2,2,*])-0.1], color=colors[1]
                  oplot, lambda[2,2,*], data[2,2,*],  color=colors[1]
                 endif else begin
                   oplot, lambda[2,2,*], data[2,2,*], psym=psyms[k], color=colors[k]
                   oplot, lambda[2,2,*], data[2,2,*],  color=colors[k]
                        
                 endelse
               
                 ;window, 1
                 ;plot, lambda[2,2,*], psym=1, yrange=[63.02, 63.35],$
                 ;      xrange=[0,17],yticklen=1.0
                 ;oplot,dblarr(22)+63.15
                 ;oplot,dblarr(22)+63.14
                 ;oplot,dblarr(22)+63.13
                 ;oplot,dblarr(22)+63.12
                 ;oplot,dblarr(22)+63.11                 
                 ;oplot,dblarr(22)+63.16
                 ;oplot,dblarr(22)+63.17
                 ;oplot,dblarr(22)+63.18
                 ;oplot,dblarr(22)+63.19
                 ;oplot,dblarr(22)+63.21                 
                 ;oplot,dblarr(22)+63.25
                 ;oplot,dblarr(22)+63.24
                 ;oplot,dblarr(22)+63.23
                 ;oplot,dblarr(22)+63.22
                 ;oplot,dblarr(22)+63.21                 
                 ;oplot,dblarr(22)+63.26
                 ;oplot,dblarr(22)+63.27
                 ;oplot,dblarr(22)+63.28
                 ;oplot,dblarr(22)+63.29                 
              endif else if i eq 2 and j eq 2 then begin
                 
                 oplot, lambda[2,2,*], psym=psyms[k], color = colors[k]
              endif
           endelse
           
           j+=1
           
        endwhile
        i+=1
        
     endwhile

     k+=1
  endwhile
  return, lambda_pix
  rlp = reform(lambda_pix, 5, 5, 16*n_ext)
  rfp = reform(flux_pix, 5, 5, 16*n_ext)
  rxs = reform(xs_pix, 5, 5, 16*n_ext)
  rys = reform(ys_pix, 5, 5, 16*n_ext)
  
  ; want a big hypercube with sorted lambdas, and fluxes for each pixel
  ; (5, 5, 2, 16*n_ext)
  hypercube = dblarr(5,5,2,16*n_ext)
  
  i=0
  while (i lt 5) do begin
     j=0
     while (j lt 5) do begin
        lp = rlp[i,j,*]
        fp = rfp[i,j,*]
        
        hypercube[i,j,0,*] = lp[sort(lp)]
        hypercube[i,j,1,*] = fp[sort(lp)]
        
        j+=1
     endwhile
     i+=1
  endwhile

  lambda_cube = reform(hypercube[*,*,0,*], 5, 5, 16*n_ext)
  flux_cube = reform(hypercube[*,*,1,*], 5, 5, 16*n_ext)
  x_cube = reform(xs_pix, 5, 5, 16*n_ext)
  y_cube = reform(ys_pix, 5, 5, 16*n_ext)
    
  ;now have 150 spaxels with interleaved lambdas
  deltalambda=0.02
  minlambda=min(lambda_pix)
  maxlambda=max(lambda_pix)

  nlambdas=(maxlambda - minlambda) / deltalambda
  lambda_grid = findgen(nlambdas+1)*deltalambda + min(lambda_pix)
;  print,'nlambdas, maxlambda, minlambda=', nlambdas, maxlambda, minlambda
;  print,'lambda grid=', lambda_grid
  ll = nlambdas + 2
;  print,'ll = ', ll

  minx=min(xs_pix)

  print,'minx here=',minx
  maxx=max(xs_pix)
  pix1=long(-floor(minx/dx))    ; using now offsets from central pix
  minx=double(-pix1*dx)
  mm=floor((maxx-minx)/dx) + 2
;  print,'minx,maxx,pix1,dx=',minx,maxx,pix1,dx
  print,'mm = ', mm
;  mm=5
  miny=min(ys_pix)
;  print, 'min y here=',miny
  maxy=max(ys_pix)
  pix2=long(-floor(miny/dy)) ; using now offsets from central pix RK 2014-03-07
  miny=double(-pix2*dy)
  nn=floor((maxy-miny)/dy) + 2
;  print,'miny,maxy,pix2,dy=',miny,maxy,pix2,dy
  print,'nn=',nn
;  nn=5
  outside=0

  psfa=dblarr(ll,mm,nn,6)

  print, 'll = ',ll
  
  lambdas = reform(lambda_pix, 400*n_ext)
  fluxes = reform(flux_pix, 400 * n_ext)
  xs_pix = reform(xs_pix, 400*n_ext)
  ys_pix = reform(ys_pix, 400*n_ext)
  
  plot,lambdas[4,*,*],fluxes[4,*,*]
  ; distribute and add up flux and weight in cube cells 
  lines=400*n_ext

  for spax=0, lines-1 do begin

     if finite(lambdas[spax]) then begin ; if wavelength is not NaN
        lambdap = lambdas[spax]	; get wavelength
        xp = xs_pix[spax]  	; get x position
        yp = ys_pix[spax]	; get y position
        flux = fluxes[spax]      ; get flux
        ;print, 'spax,lambdap, xp, yp,flux=',spax,lambdap, xp, yp,flux
        ;sd = data[spax,5,i]		; get standard dev. (sigma)
        sd=0
        left = Floor((xp - minx)/ dx) ; + 1 	; define output x spaxex
        right = left + 1              ; define output x spaxex
        down = Floor((yp - miny) / dy) ; + 1	; define output y spaxex
        up = down + 1                  ; define output y index
        blue = Floor((lambdap - minlambda)/deltalambda) ; + 1	; define output 
        red = blue + 1                                  ; wavelength index
        xr = minx + right * dx	; calculate x values @ output pix boundaries
        xl = xr - dx 
        yu = miny + up * dy     ; calculate y values @ output pix boundaries
        yd = yu - dy
        lambdar = minlambda + red * deltalambda	; calculate wavelengths @ output
        lambdab = lambdar - deltalambda         ; pix boundaries

        print,'xp, yp, lambdap, flux=',xp,yp,lambdap,flux
        print,'blue, red, left, right, up, down, lambdar, lambdab=', blue, red, left, right, up, down, lambdar, lambdab
        print,'xr, xl, yu, yd=',xr,xl,yu,yd
        print,'lambdar, lambdab=',lambdar, lambdab

        ; drizzle only flux inside the output array
        if (blue ge 0) and (blue lt ll) and (left ge 0) and (left lt mm) and $
           (up ge 0) and (up lt nn) then begin
           weight=(xr-xp)*(yp-yd)*(lambdar-lambdap)
           print,'weight1=',weight
           ; flux
           psfa[blue,left,up,3]=psfa[blue,left,up,3]+flux*weight
           print, 'flux*weight1=', flux*weight
           ; weight
           psfa[blue,left,up,4]=psfa[blue,left,up,4]+weight
           ;psfa[blue,left,up,5]=psfa[blue,left,up,5]+weight^2*sd^2
        endif else outside=1 

        if (blue ge 0) and (blue lt ll) and (left ge 0) and (left lt mm) and $
           (down ge 0) and (down lt nn) then begin
           weight=(xr-xp)*(yu-yp)*(lambdar-lambdap)
           print,'weight2=',weight
           psfa[blue,left,down,3]=psfa[blue,left,down,3]+flux*weight
           psfa[blue,left,down,4]=psfa[blue,left,down,4]+weight
           ;psfa[blue,left,down,5]=psfa[blue,left,down,5]+weight^2*sd^2
        endif else outside=1 

        if (blue ge 0) and (blue lt ll) and (right ge 0) and (right lt mm) and $
           (up ge 0) and (up lt nn) then begin
           weight=(xp-xl)*(yp-yd)*(lambdar-lambdap)
           print,'weight3=',weight
           psfa[blue,right,up,3]=psfa[blue,right,up,3]+flux*weight
           psfa[blue,right,up,4]=psfa[blue,right,up,4]+weight
           ;psfa[blue,right,up,5]=psfa[blue,right,up,5]+weight^2*sd^2
        endif else outside=1 

        if (blue ge 0) and (blue lt ll) and (right ge 0) and (right lt mm) and $
           (down ge 0) and (down lt nn) then begin
           weight=(xp-xl)*(yu-yp)*(lambdar-lambdap)
           print,'weight4=',weight
           psfa[blue,right,down,3]=psfa[blue,right,down,3]+flux*weight
           psfa[blue,right,down,4]=psfa[blue,right,down,4]+weight
           ;psfa[blue,right,down,5]=psfa[blue,right,down,5]+weight^2*sd^2
        endif else outside=1 

        if (red ge 0) and (red lt ll) and (left ge 0) and (left lt mm) and $
           (up ge 0) and (up lt nn) then begin
           weight=(xr-xp)*(yp-yd)*(lambdap-lambdab)
           print,'weight5=',weight
           psfa[red,left,up,3]=psfa[red,left,up,3]+flux*weight
           psfa[red,left,up,4]=psfa[red,left,up,4]+weight
           ;psfa[red,left,up,5]=psfa[red,left,up,5]+weight^2*sd^2
        endif else outside=1 

        if (red ge 0) and (red lt ll) and (right ge 0) and (right lt mm) and $ 
           (up ge 0) and (up lt nn) then begin
           weight=(xp-xl)*(yp-yd)*(lambdap-lambdab)
           print,'weight6=',weight
           psfa[red,right,up,3]=psfa[red,right,up,3]+flux*weight
           psfa[red,right,up,4]=psfa[red,right,up,4]+weight
           ;psfa[red,right,up,5]=psfa[red,right,up,5]+weight^2*sd^2
        endif else outside=1 

        if (red ge 0) and (red lt ll) and (left ge 0) and (left lt mm) and $ 
           (down ge 0) and (down lt nn) then begin
           weight=(xr-xp)*(yu-yp)*(lambdap-lambdab)
           print,'weight7=',weight
           psfa[red,left,down,3]=psfa[red,left,down,3]+flux*weight
           psfa[red,left,down,4]=psfa[red,left,down,4]+weight
           ;psfa[red,left,down,5]=psfa[red,left,down,5]+weight^2*sd^2
        endif else outside=1 

        if (red ge 0) and (red lt ll) and (right ge 0) and (right lt mm) and $
           (down ge 0) and (down lt nn) then begin
           weight=(xp-xl)*(yu-yp)*(lambdap-lambdab)
           print,'weight8=',weight
           psfa[red,right,down,3]=psfa[red,right,down,3]+flux*weight 
           psfa[red,right,down,4]=psfa[red,right,down,4]+weight
           ;psfa[red,right,down,5]=psfa[red,right,down,5]+weight^2*sd^2
        endif else outside=1 

     endif
  endfor

  ;return, psfa
  if outside eq 1 then $
     print, 'fifi_ls_drizzle.pro: some flux fell outside the output grid.'

  ;psfa[*,*,*,5]=sqrt(psfa[*,*,*,5])   ; simga

  ; normalize accumulated flux by accumulated weight
  q=where(psfa[*, *, *, 4] gt 0) ; select non-zero weights
  if q[0] ne -1 then begin
     q=array_indices(psfa,q)
     for i=0L,n_elements(q)/4-1 do begin 
        psfa[q[0,i],q[1,i],q[2,i],3]=psfa[q[0,i],q[1,i],q[2,i],3]/psfa[q[0,i],q[1,i],q[2,i],4] ; flux
        ;    psfa[q[0,i],q[1,i],q[2,i],5]=psfa[q[0,i],q[1,i],q[2,i],5]/psfa[q[0,i],q[1,i],q[2,i],4]  ; sigma 
     endfor
  endif else $ 
     print,'drizzle.pro: all weights were zero... something is very wrong here.'

  ; this is what got written
  cube=float(transpose(psfa[*,*,*,3:5],[1,2,0,3])) ; [x,y,lambda,results]

  
  x0=fxpar(primehead,'OBSLAM')
  y0=fxpar(primehead,'OBSBET')  
  rot = fxpar(primehead,'DET_ANGL')
;want ra in deg = hour * 15 + min / 4 + sec /240
  DA      = fxpar(primehead,'DET_ANGL') * !DTOR ; ccw array rotation
  make_astr,ast,crpix=[pix1+1,pix2+1],crval=[x0,y0],$
            cd=[[-cos(DA),-sin(DA)],[sin(DA),cos(DA)]]*dx/3600d0,$
            ctype=['RA---TAN','DEC--TAN']
  putast,primehead,ast,EQUINOX=2000

  fxaddpar, primehead, 'CRVAL3', string(minlambda+(maxlambda-minlambda)/2), 'LAMBDA at ref. pixel'

  fxaddpar, primehead, 'CRPIX3', string(nlambdas/2), 'Reference pixel number for LAMBDA' 

  fxaddpar,primehead,'CTYPE3','LAMBDA' 
  fxaddpar,primehead,'CDELT3',string(deltalambda),'LAMBDA grid spacing'
  fxaddpar,primehead,'PROCSTAT','LEVEL_3' ; update SOFIA keyword
 
  outfilename=infile+'.drizzle.L3.fits'

  fxwrite, outfilename, primehead, /noupdate
  mwrfits, cube, outfilename
                                

  flx_data_cube = reform(cube(*,*,*,0))
  fluxerfilename=strmid(outfilename,0,strlen(outfilename)-4)+'drizzle_flux.fits'
  fxwrite,fluxerfilename, primehead, /noupdate
  mwrfits,flx_data_cube,fluxerfilename

  return, flx_data_cube
end

