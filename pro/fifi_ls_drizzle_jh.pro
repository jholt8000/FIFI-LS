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
; OPTIONAL INPUTS:new_n_pix_x = new_n_pix_x, $
;                             new_n_pix_y = new_n_pix_y, $
;                             delta_x_out = delta_x_out, $
;                             delta_y_out = delta_y_out, y_shift = y_shift, $
;                             x_shift = x_shift, delta_l_out=delta_l_out
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

;-
function fifi_ls_drizzle_jh, infile, new_n_pix_x = new_n_pix_x, $
                             new_n_pix_y = new_n_pix_y, $
                             delta_x_out = delta_x_out, $
                             delta_y_out = delta_y_out, y_shift = y_shift, $
                             x_shift = x_shift, delta_l_out=delta_l_out

  if ~keyword_set(infile) then print, 'Cannot reduce without filename'
  if ~keyword_set(new_n_pix_x) then new_n_pix_x = 5
  if ~keyword_set(new_n_pix_y) then new_n_pix_y = 5

  if ~keyword_set(x_shift) then x_shift = 0
  if ~keyword_set(y_shift) then y_shift = 0
  if ~keyword_set(delta_l_out) then delta_l_out = 0.02
  
  primehead = headfits(infile)
  channel = strtrim(fxpar(primehead, 'CHANNEL'),2)

  print,'channel=',channel
  if (channel eq 'RED') then begin
     deltax = 12L
     deltay = 12L
     if ~keyword_set(delta_x_out) then delta_x_out = 12L
     if ~keyword_set(delta_y_out) then delta_y_out = 12L
  endif else begin
     deltax = 6L
     deltay = 6L
     if ~keyword_set(delta_x_out) then delta_x_out = 6L
     if ~keyword_set(delta_y_out) then delta_y_out = 6L
  endelse

  dlam_map = fxpar(primehead, 'DLAM_MAP')
  dbet_map = fxpar(primehead, 'DBET_MAP')

  shift_map = 1
  if shift_map then begin
      if abs(dlam_map) gt 0 then begin
         x_shift = -dlam_map / (delta_x_out)
      endif
  
      if abs(dbet_map) gt 0 then begin
         y_shift = dbet_map / float((delta_y_out))
      endif
  endif
  
  in_struct = mrdfits(infile, 1, /silent)
  fluxes = in_struct.data
  lambdas = in_struct.lambda
  xes = in_struct.xs
  yes = in_struct.ys
  stddev = in_struct.stddev

  lambdas = reform(lambdas, n_elements(lambdas))
  fluxes = reform(fluxes, n_elements(fluxes))
  xes = reform(xes, n_elements(xes))
  yes = reform(yes, n_elements(yes))

  new_grid_size_x = new_n_pix_x * delta_x_out
  new_grid_size_y = new_n_pix_y * delta_y_out
  
  old_grid_size_x = 5 * deltax
  old_grid_size_y = 5 * deltay
 
  extra_pix_x = (new_grid_size_x - old_grid_size_x ) / delta_x_out
  extra_pix_y = (new_grid_size_y - old_grid_size_y ) / delta_y_out

  print, 'extra_pix_x=', extra_pix_x
  print, 'extra_pix_y=',extra_pix_y
  ;final grid size, at minimum need to add a deltax/2. to each end
  ; findgen doesn't include last point so need to add one more for 
  ; that and two at each end for a total of 3
  ;xes = xes - x_shift
  ;print, 'shifted xes by ', x_shift
  xmin = min(xes)
  xmax = max(xes)
  x_grid = (findgen(new_n_pix_x) * delta_x_out) + xmin - $
           delta_x_out * extra_pix_x/2.  ;- (deltax/2.)
  ymin=min(yes)
  ymax=max(yes)
  y_grid = (findgen(new_n_pix_y) * delta_y_out) + ymin - $
           delta_y_out * extra_pix_y/2. ;- (deltay/2.)

  
  ;if (channel eq 'RED') then $
 ; x_grid = reverse(x_grid) + delta_x_out / 2.

  y_grid = reverse(y_grid) + delta_y_out / 2.
  minlambda=min(lambdas)
  maxlambda=max(lambdas)
  nlambdas=floor((maxlambda - minlambda) / delta_l_out)
  lambda_grid = findgen(nlambdas+3)*delta_l_out + minlambda - (delta_l_out / 2.)
  ll = n_elements(lambda_grid)
  
  big_flux_grid = dblarr(new_n_pix_x, new_n_pix_y, ll)
  n_times = dblarr(new_n_pix_x, new_n_pix_y, ll)
  weights = dblarr(new_n_pix_x, new_n_pix_y, ll)
  fracs = dblarr(new_n_pix_x, new_n_pix_y, ll)
  
  ; do the drizzle
  for spax=0, (n_elements(fluxes))-1 do begin

     flux = fluxes[spax]
     xstart=0
     xend=0
     ystart=0
     yend=0
     lambdastart=0
     lambdaend=0
     lhs_x=0
     rhs_x=0
     lhs_y=0
     rhs_y=0
     lhs_lambda=0
     rhs_lambda=0
     
     ; first go through x
     x0 = xes[spax] - deltax / 2.
     x1 = xes[spax] + deltax / 2.
     ; this is causing a lot of overlap, is that what we want?
     xcenter = where(x_grid ge x0 and x_grid le x1, count)
     
     if count ne 0 then begin
        if xcenter[0] lt 1 then begin
           xstart = 0
        endif else begin
           xstart = xcenter[0] - 1
        endelse
        if xcenter[-1] ge n_elements(x_grid)-1 then begin
           xend = xcenter[-1]
        endif else begin
           xend = xcenter[-1] + 1
        endelse  

        wtd_x = x_grid[xstart:xend]
        lhs_x = dblarr(n_elements(wtd_x))
        rhs_x = dblarr(n_elements(wtd_x))
        for j = 0, n_elements(wtd_x) - 1 do begin
           if (x1 lt (wtd_x[j] - delta_x_out / 2.)) then begin
              lhs_x[j] = 0
              rhs_x[j] = 0
           endif else if ((x0 lt (wtd_x[j] - delta_x_out / 2.)) and $
                           ((wtd_x[j] - delta_x_out / 2.) lt x1) and $
                           (x1 lt wtd_x[j])) then begin
              lhs_x[j] = (delta_x_out / 2.) - (wtd_x[j] - x1)
              rhs_x[j] = 0
            endif else if ((x0 lt (wtd_x[j] - delta_x_out / 2.) and $
                           (wtd_x[j] lt x1) and $
                           (x1 lt (wtd_x[j] + delta_x_out / 2.)))) then begin
              lhs_x[j] = delta_x_out / 2.
              rhs_x[j] = x1 - wtd_x[j]
           endif else if ((x0 lt (wtd_x[j] - delta_x_out / 2.) and $
                           ((wtd_x[j] + delta_x_out / 2.) lt x1))) then begin
              lhs_x[j] = delta_x_out / 2.
              rhs_x[j] = delta_x_out / 2.
           ;case c
           endif else if ((x0 lt wtd_x[j]) and $
                          ((wtd_x[j] - delta_x_out / 2.) lt x0) and $
                          (wtd_x[j] lt x1) and $
                          (x1 lt (wtd_x[j] + (delta_x_out / 2.)))) then begin
              lhs_x[j] = wtd_x[j] - x0
              rhs_x[j] = x1 - wtd_x[j]
           endif else if ((x0 lt wtd_x[j]) and $
                          ((wtd_x[j] - delta_x_out / 2.) lt x0) and $
                          ((wtd_x[j] + (delta_x_out / 2.)) lt x1)) then begin
              lhs_x[j] = wtd_x[j] - x0
              rhs_x[j] = delta_x_out / 2.
           endif else if ((wtd_x[j] lt x0) and $
                          (x0 lt (wtd_x[j] + (delta_x_out / 2.)) and $
                           ((wtd_x[j] + (delta_x_out / 2.)) lt x1))) then begin
              lhs_x[j] = 0
              rhs_x[j] = wtd_x[j] + (delta_x_out / 2.) - x0
           endif else if ((wtd_x[j] + delta_x_out / 2.) lt x0) then begin
              lhs_x[j] = 0
              rhs_x[j] = 0
           endif
        endfor
     endif
     ;;
     ;; Go through y, do not need to do this if lhs_x, lhs_y both eq 0
     y0 = yes[spax] - deltay / 2.
     y1 = yes[spax] + deltay / 2.
     ycenter = where(y_grid ge y0 and y_grid le y1, ycount)
     if ycount ne 0 then begin
        if ycenter[0] lt 1 then begin
           ystart = 0
        endif else begin
           ystart = ycenter[0]-1
        endelse
        if ycenter[-1] ge n_elements(y_grid)-1 then begin
           yend = ycenter[-1]
        endif else begin
           yend = ycenter[-1] + 1
        endelse  

        wtd_y = y_grid[ystart:yend]
        lhs_y = dblarr(n_elements(wtd_y))
        rhs_y = dblarr(n_elements(wtd_y))
        for j = 0, n_elements(wtd_y) - 1 do begin
           ; case A
           if (y1 lt (wtd_y[j] - delta_y_out / 2.)) then begin
              lhs_y[j] = 0
              rhs_y[j] = 0
           endif else if ((y0 lt (wtd_y[j] - delta_y_out / 2.)) and $
                           ((wtd_y[j] - delta_y_out / 2.) lt y1) and $
                           (y1 lt wtd_y[j])) then begin
              lhs_y[j] = (delta_y_out / 2.) - (wtd_y[j] - y1)
              rhs_y[j] = 0
            endif else if ((y0 lt (wtd_y[j] - delta_y_out / 2.) and $
                           (wtd_y[j] lt y1) and $
                           (y1 lt (wtd_y[j] + delta_y_out / 2.)))) then begin
              lhs_y[j] = delta_y_out / 2.
              rhs_y[j] = y1 - wtd_y[j]
           endif else if ((y0 lt (wtd_y[j] - delta_y_out / 2.) and $
                           ((wtd_y[j] + delta_y_out / 2.) lt y1))) then begin
              lhs_y[j] = delta_y_out / 2.
              rhs_y[j] = delta_y_out / 2.
           ;case c
           endif else if ((y0 lt wtd_y[j]) and $
                          ((wtd_y[j] - delta_y_out / 2.) lt y0) and $
                          (wtd_y[j] lt y1) and $
                          (y1 lt (wtd_y[j] + (delta_y_out / 2.)))) then begin
              lhs_y[j] = wtd_y[j] - y0
              rhs_y[j] = y1 - wtd_y[j]
           endif else if ((y0 lt wtd_y[j]) and $
                          ((wtd_y[j] - delta_y_out / 2.) lt y0) and $
                          ((wtd_y[j] + (delta_y_out / 2.)) lt y1)) then begin
              lhs_y[j] = wtd_y[j] - y0
              rhs_y[j] = delta_y_out / 2.
           endif else if ((wtd_y[j] lt y0) and $
                          (y0 lt (wtd_y[j] + (delta_y_out / 2.)) and $
                           ((wtd_y[j] + (delta_y_out / 2.)) lt y1))) then begin
              lhs_y[j] = 0
              rhs_y[j] = wtd_y[j] + (delta_y_out / 2.) - y0
           endif else if ((wtd_y[j] + delta_y_out / 2.) lt y0) then begin
              lhs_y[j] = 0
              rhs_y[j] = 0
           endif
        endfor
     endif

     ;; lambda, do not need to do this if lhs_y, rhs_y eq 0
     lambda0 = lambdas[spax] - delta_l_out / 2.
     lambda1 = lambdas[spax] + delta_l_out / 2.
     lambdacenter = where(lambda_grid ge lambda0 and lambda_grid le lambda1, lcount)
    
     if lcount ne 0 then begin
        if lambdacenter[0] lt 1 then begin
           lambdastart = 0
        endif else begin
           lambdastart = lambdacenter[0]-1
        endelse
        if lambdacenter[-1]+1 gt n_elements(lambda_grid) then begin
           lambdaend = lambdacenter[-1]
        endif else begin
           lambdaend = lambdacenter[-1] + 1
        endelse  

        wtd_lambda = lambda_grid[lambdastart:lambdaend]
        lhs_lambda = dblarr(n_elements(wtd_lambda))
        rhs_lambda = dblarr(n_elements(wtd_lambda))

        for j=0, n_elements(wtd_lambda)-1 do begin
           ; case A
           if (lambda1 lt (wtd_lambda[j] - delta_l_out / 2.)) then begin
              lhs_lambda[j] = 0
              rhs_lambda[j] = 0
           endif else if ((lambda0 lt (wtd_lambda[j] - delta_l_out / 2.)) and $
                           ((wtd_lambda[j] - delta_l_out / 2.) lt lambda1) and $
                           (lambda1 lt wtd_lambda[j])) then begin
              lhs_lambda[j] = (delta_l_out / 2.) - (wtd_lambda[j] - lambda1)
              rhs_lambda[j] = 0
            endif else if ((lambda0 lt (wtd_lambda[j] - delta_l_out / 2.) and $
                           (wtd_lambda[j] lt lambda1) and $
                           (lambda1 lt (wtd_lambda[j] + delta_l_out / 2.)))) then begin
              lhs_lambda[j] = delta_l_out / 2.
              rhs_lambda[j] = lambda1 - wtd_lambda[j]
           endif else if ((lambda0 lt (wtd_lambda[j] - delta_l_out / 2.) and $
                           ((wtd_lambda[j] + delta_l_out / 2.) lt lambda1))) then begin
              lhs_lambda[j] = delta_l_out / 2.
              rhs_lambda[j] = delta_l_out / 2.
           ;case c
           endif else if ((lambda0 lt wtd_lambda[j]) and $
                          ((wtd_lambda[j] - delta_l_out / 2.) lt lambda0) and $
                          (wtd_lambda[j] lt lambda1) and $
                          (lambda1 lt (wtd_lambda[j] + (delta_l_out / 2.)))) then begin
              lhs_lambda[j] = wtd_lambda[j] - lambda0
              rhs_lambda[j] = lambda1 - wtd_lambda[j]
           endif else if ((lambda0 lt wtd_lambda[j]) and $
                          ((wtd_lambda[j] - delta_l_out / 2.) lt lambda0) and $
                          ((wtd_lambda[j] + (delta_l_out / 2.)) lt lambda1)) then begin
              lhs_lambda[j] = wtd_lambda[j] - lambda0
              rhs_lambda[j] = delta_l_out / 2.
           endif else if ((wtd_lambda[j] lt lambda0) and $
                          (lambda0 lt (wtd_lambda[j] + (delta_l_out / 2.)) and $
                           ((wtd_lambda[j] + (delta_l_out / 2.)) lt lambda1))) then begin
              lhs_lambda[j] = 0
              rhs_lambda[j] = wtd_lambda[j] + (delta_l_out / 2.) - lambda0
           endif else if ((wtd_lambda[j] + delta_l_out / 2.) lt lambda0) then begin
              lhs_lambda[j] = 0
              rhs_lambda[j] = 0
           endif
        endfor
     endif

     debug = 1
     if debug then begin
        print, 'x, y, lambda at spax=', xes[spax], yes[spax], lambdas[spax], spax
        print, 'x0, x1 =', x0, x1,' xcenter=', xcenter, 'xstart=',xstart,' xend=',xend
        print,' wtd_x = ',wtd_x
        print,'x grid=',x_grid
        print, ' lhs_x, rhs_x=', lhs_x, rhs_x
        ;print, 'y0, y1, lhs_y, rhs_y=', y0, y1, lhs_y, rhs_y
        ;print, 'lambda0, lambda1, lhs_lambda, rhs_lambda=', lambda0, lambda1, $
        ;       lhs_lambda, rhs_lambda
        
     endif
     
     ;dont need to go through full grid, know lambdastart, lambdaend, etc
                                ; now drizzle the flux down on the new grid
     ii=0
     print,'xstart=',xstart, 'xend=',xend
     ;print,'ystart=',ystart, 'yend=',yend
    
     i = xstart + x_shift
     while i lt xend + x_shift-1 do begin
     ;for i = xstart, xend do begin
        jj=0
        j = ystart + y_shift
        while j lt yend + y_shift-1 do begin
        ;for j = ystart, yend do begin
           kk=0
           for k=lambdastart, lambdaend do begin

              orig_volume = (x1-x0) * (y1-y0) * (lambda1 - lambda0)
              volume = ((lhs_x[ii] + rhs_x[ii]) * (lhs_y[jj] + rhs_y[jj]) * $
                                          (lhs_lambda[kk] + rhs_lambda[kk]))
             
              if volume gt 0 and finite(flux) and (flux gt 0) then begin
                 
                          big_flux_grid[i, j, k] = big_flux_grid[i, j, k] + (flux * volume / orig_volume)
                          if (big_flux_grid[i, j, k] gt 10) and (debug eq 1) then begin
                             print, 'volume=',volume,'  i=',i,' j=',j,' k=',k,' ii=', $
                                    ii,' jj=',jj,' kk=',kk, ' big_flux_grid[i,j,k]=', $
                                    big_flux_grid[i, j, k]
                             print,'flux=',flux,' volume=',volume, '  orig_volume=', $
                                   orig_volume
                             
                             print,'x0,x1,y0,y1,lambda0,lambda1=', x0, x1, y0, y1, $$$$
                                   lambda0,lambda1
                              
                              print,'xgrid=',x_grid
                              print,'xmin=',xmin,' xmax=',xmax
                           endif

                          ; need to record the fraction of a single
                          ; pixel included in each final pixel,
                          ; that will be the weighting factor when added
                          
                          n_times[i, j, k] = n_times[i, j, k] + 1
                          weights[i, j, k] = volume / orig_volume
                          ;shifted_grid[i-dlam_map,j,k] = big_flux_grid[i,j,k]
               endif

               kk+=1
                       
           endfor
           jj+=1

           j+=1
        endwhile
        
        ;endfor
        ii+=1
        i+=1
      endwhile
     ;endfor
     ;print,'---spax = ', spax,' done -------------------------'

  endfor

  window,1
  !X.MARGIN = !X.MARGIN/2.
  !Y.MARGIN = !Y.MARGIN/2.
  !P.MULTI=[0, new_n_pix_x, new_n_pix_y]

  ; have to switch i, j -> j, i because multi does col first
  loadct, 0, /silent
  print,'size ntimes=',size(n_times)
  
  for j=new_n_pix_y-1,0,-1 do begin
     for i=0, new_n_pix_x-1 do  begin
        for k=0, nlambdas+2 do begin
           if n_times[i, j, k] gt 0 then begin
              big_flux_grid[i, j, k] = (big_flux_grid[i, j, k] / n_times[i, j, k]) ; / weights[i, j, k]
              ; fraction of a single pixel in new pixel

              fracs[i, j, k] = weights[i, j, k] * n_times[i, j, k]
           endif
           

        endfor
        plot, lambda_grid, big_flux_grid[i, j, *], psym=2, yrange=[0,40]
        endfor
  endfor

  device, decomposed=0
  loadct, 13, /silent
  !P.MULTI = 0
  window, 3
  for i=0, n_elements(x_grid) - 1 do begin
     for j=0, n_elements(y_grid) - 1 do begin
       
     x0 = x_grid[i] - delta_x_out/2. 
     y0 = y_grid[j] - delta_y_out/2. 
     xlength = delta_x_out
     ylength = delta_y_out
     if (i eq 0) and (j eq 0) then begin
         plot, [x_grid[i]], [y_grid[j]], psym=6, linestyle=3, color=100, xrange=[-30, 30], yrange=[-30, 30]
         oplot, [x0,x0+xlength,x0+xlength,x0,x0],[y0,y0,y0+ylength,y0+ylength,y0],color=100
      endif else begin
         oplot, [x_grid[i]], [y_grid[j]], psym=6,color=100
         oplot, [x0,x0+xlength,x0+xlength,x0,x0],[y0,y0,y0+ylength,y0+ylength,y0],color=100
      endelse
   endfor     
  endfor

  rxs = reform(xes, 5, 5, n_elements(xes)/25)
  rys = reform(yes, 5, 5, n_elements(yes)/25)
  for i=0, n_elements(rxs[*,*,1,*])-1 do begin
     xs=rxs[*,*,1,*]
     ys=rys[*,*,1,*]
     x0=xs[i]- deltax / 2.
     y0=ys[i]- deltay / 2.
     xlength = deltax
     ylength = deltay
     if (i eq 0)  then begin
         oplot, [xs[i]], [ys[i]], psym=4, color=255
         oplot, [x0,x0+xlength,x0+xlength,x0,x0],[y0,y0,y0+ylength,y0+ylength,y0],color=255
      endif else begin
         oplot, [xs[i]], [ys[i]], psym=4,color=255
         oplot, [x0,x0+xlength,x0+xlength,x0,x0],[y0,y0,y0+ylength,y0+ylength,y0],color=255
      endelse
   endfor         

  
  x0=fxpar(primehead,'OBSLAM')
  y0=fxpar(primehead,'OBSBET')  
  rot = fxpar(primehead,'DET_ANGL')
  ;want ra in deg = hour * 15 + min / 4 + sec /240
  DA      = fxpar(primehead,'DET_ANGL') * !DTOR ; ccw array rotation
  ;make_astr,ast,crpix=[pix1+1,pix2+1],crval=[x0,y0],$
  ;          cd=[[-cos(DA),-sin(DA)],[sin(DA),cos(DA)]]*dx/3600d0,$
  ;          ctype=['RA---TAN','DEC--TAN']
  ;putast,primehead,ast,EQUINOX=2000

  fxaddpar, primehead, 'CRVAL3', string(minlambda+(maxlambda-minlambda)/2), 'LAMBDA at ref. pixel'
  fxaddpar, primehead, 'CRPIX3', string(nlambdas/2), 'Reference pixel number for LAMBDA' 

  fxaddpar,primehead,'CTYPE3','LAMBDA' 
  fxaddpar,primehead,'CDELT3',string(delta_l_out),'LAMBDA grid spacing'
  fxaddpar,primehead,'PROCSTAT','LEVEL_3' ; update SOFIA keyword


  outfilename=str_replace(infile, '.scancomb.fits', '.drizzle.fits')

  ; if the output file already exists, delete it
  if file_test(outfilename) gt 0 then file_delete, outfilename

  fxwrite, outfilename, primehead, /noupdate, /append
  mwrfits, big_flux_grid, outfilename, /silent
  mwrfits, fracs, outfilename
           
  ;flx_data_cube = reform(cube(*,*,*,0))
  ;fluxerfilename=strmid(outfilename,0,strlen(outfilename)-4)+'drizzle_flux.fits'
  ;fxwrite,fluxerfilename, primehead, /noupdate
  ;mwrfits,flx_data_cube,fluxerfilename

  ;return, n_times
end


pro add_drizzle, files

  added = 0
  big_cube = dblarr(100, 100, 20)
  primehead = headfits(files[0])
  tfracs=0
  ;hc = dblarr(100, 100, 20, n_elements(files))
  
  for i = 0, n_elements(files) - 1 do begin
     d = mrdfits(files[i], 1, /silent)
     fracs= mrdfits(files[i], 2, /silent)
     
     big_cube = big_cube + d * fracs
     ;tfracs = tfracs + 1
     ; need to go through i,j,k and find the zeros and the ones with <1contributor
     added = added + 1
     ;print,'sized=',size(d)
     ;print,'sizehc=',size(hc[*,*,*,i]),i
     ;hc[*,*,*,i] = d
     ;print,'sd=',size(d)
     
  end

  big_cube = big_cube / added

  ;big_cube = hc
  fxwrite, 'add_driz.fits', primehead, /noupdate
  mwrfits, big_cube / added, 'add_driz.fits'
end

     
