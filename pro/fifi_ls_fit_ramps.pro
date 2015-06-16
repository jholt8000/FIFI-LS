;+
; NAME: fifi_ls_fit_ramps 
;
; PURPOSE: Fit straight lines to ramps 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=fit_ramps(FILENAME, OUTDIR, s2n=s2n,
;                                    verbose=verbose) 
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
; OUTPUTS: RESULT: 1. Returns the last extension fits structure
;                  2. fits file containing linear fit (a+bx) result and sigma
;                     File name: *.chop#.CLEAN.rampfits.fits
;                     fits content:
;                      prime header
;                      extensions for each grating with the following:
;                      ext.header = extension header
;                      ext.data = flux cube, shape=(5,5,18)
;                      ext.stddev = err cube, shape=(5,5,18)
; 
; OPTIONAL OUTPUTS:
;
; SIDE EFFECTS: 
;               1. Removes the 26th "sky" pixel if it has not 
;                  already been removed 
;               2. Data points in bent ramps will be not used in line fitting.
;               3. If choplength is longer than ramplength, the ramps on the
;                  same chop will be averaged before line fitting. 
;		4. Garbage ramps (as defined as averaged ramps or
;		sub-ramps having less than 3 data points to fit) are
;		assigned NaN for their fit parameters and sigma.  
;
; PROCEDURE: 1. Average ramps on single inductosyn and chop position
;            2. Flag saturated data points; fit a straight line through
;               remaining points on each (sub)ramp; calculate error.
;            3. Create fits file and write results to disk.
;
; EXAMPLE: IDL> result=fit_ramps('00006_rawlw021.L1.chop0.CLEAN.fits')
;          IDL> fluxer, result.data
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
; 2015feb17 jholt:  Fix typo in check_scn call and update mrdfits to use
;               filename instead of file unit to work with newer
;               astrolib.
;
; 2015Mar22 jholt: Separate indpos, use max peak index instead of
;                  hardcoded number to find bent ramps

function fifi_ls_fit_ramps, filename, outdir, s2n = s2n, verbose = verbose

  comb_before_fit = 'no'
  
  if ~keyword_set(outdir) then outdir = '.'
  if ~keyword_set(s2n) then s2n = 10
  if ~keyword_set(verbose) then verbose = 0

  ;  make sure the input file exists
  if (file_search(filename, /full))[0] eq '' then begin
     print, 'fifi_ls_fit_ramps.pro: input file does not exist. Exiting.'
     return, -1
  endif

  ;  define parameters
  ARRAY_SIZE_X = 5
  ARRAY_SIZE_Y = 5
  ARRAY_SIZE_LAM = 18
  ARRAY_SIZE = (ARRAY_SIZE_X*ARRAY_SIZE_Y+1)*ARRAY_SIZE_LAM

  ;  read input data file header
  primehead = headfits(filename)

  ; -- define output structure and file name --
  fitfilename = strmid(filename, 0, strlen(filename) - 4)+'rampfits.fits'

  C_SCHEME = fxpar(primehead, 'C_SCHEME')
  CHOP_LN = fxpar(primehead, 'C_CHOPLN')
  CHANNEL = strtrim(fxpar(primehead, 'CHANNEL'))

  if (CHANNEL eq 'RED') then begin
     ramplength = fxpar(primehead, 'RAMPLN_R')
     n_ramps_to_remove = 1
  endif else begin
     ramplength = fxpar(primehead, 'RAMPLN_B')
     n_ramps_to_remove = 2
  endelse
  
  ; get number of inductosyn extensions
  fits_info, filename, N_ext = n_ext, /silent

  ; if the output file already exists, delete it
  if file_test(fitfilename) gt 0 then file_delete, fitfilename

  ; write the file with the primary header and make it appendable 
  fxwrite, fitfilename, primehead, /noupdate, /append

  flux = fltarr(ARRAY_SIZE_X*ARRAY_SIZE_Y, ARRAY_SIZE_LAM)
  stddev_arr = fltarr(ARRAY_SIZE_X*ARRAY_SIZE_Y, ARRAY_SIZE_LAM)

  k = 1
  while (k lt n_ext+1) do begin

     ; read the data extension
     raw_data = mrdfits(filename, k, exthdr, /unsigned, /silent)

     ; remove the 26th pixel
     raw_data = raw_data[0:24, *, *, *]

     total_readouts = n_elements(raw_data) ; number of readouts
     nramps = total_readouts/ramplength    ; number of ramps, total
     nramps_per_spaxel = nramps / 450.

     ;reshape the data into separate ramp dimensions (takes 18, 25, 1, nramps)
     ; into (18, 25, ramplength, nramps_per_spaxel)

     raw_data2 = reform(raw_data, 25, 18, ramplength, nramps_per_spaxel)

     ; remove first ramp from each spaxel and remove first and last readout
     ; in each ramp
     raw_data2 = raw_data2[*, *, 1:ramplength - 2,$
                           n_ramps_to_remove:nramps_per_spaxel-1]
     
     x = findgen(ramplength - 2)
     
     fit = fltarr(3, ARRAY_SIZE_X*ARRAY_SIZE_Y, ARRAY_SIZE_LAM)

     if comb_before_fit eq 'yes' then begin
        for i = 0, 24 do begin
           for j = 0, 17 do begin
              
              ; average each set of ramps
              avg_ramps = mean(raw_data2, dimension = 4)
              stddev_ramps = stddev(raw_data2, dimension = 4)
              ; stddev_ramps will be used for the correlated noise 

              ;uncomment to see the ramp plotting
              ;plot,raw_data2[i,j,*,*],psym = 4
              ;oplot,avg_ramps[i,j,*]
              
              rmax = max(avg_ramps[i, j, *], locmax)
              
              ; if the ramp is bent, then the maximum will not be at the
              ; last position
              ; do not use the readout right before the bend
              ; we cannot know when it went nonlinear -- maybe remove more?
              if locmax lt ramplength-n_ramps_to_remove - 3 then begin
                 q_end = locmax-1
              endif else q_end = ramplength-3
                 ; use the ramplength-3 is the first and last discarded readout
                 ; along with the maximum of ramp that we will not use for
                 ; bent ramps
              
              if q_end gt 0 then begin
                 ;y = A+Bx
                 fit[0:1, i, j] = linfit(x[0:q_end], avg_ramps[i, j, 0:q_end], $
                                         /double,chisq = chisq, sigma = sigma) 
                 ; sigma for slope
                 fit[2, i, j] = sigma[1] 

                 if s2n gt 0 then begin
                    if fit[1, i, j] / fit[2, i, j] lt s2n then $
                       fit[1, i, j] = !VALUES.F_NAN
                 endif
              endif else fit[1, i, j] = !VALUES.F_NAN

              flux[i, j] = fit[1, i, j]
              stddev_arr[i, j] = fit[2, i, j]
           endfor
        endfor

     endif else begin ; end comb before fit, begin fit before comb
        
        fit = fltarr(3, ARRAY_SIZE_X*ARRAY_SIZE_Y, ARRAY_SIZE_LAM, $
                     n_elements(raw_data2[0, 0, 0, *]))

        for i = 0, 24 do begin
           for j = 0, 17 do begin
              for kk = 0, n_elements(raw_data2[0, 0, 0, *])-1 do begin
                 ;if kk eq 0  and j ne 0 and j ne 17 then begin
                     ;   plot,  raw_data2[i, j,  *, kk], psym = 4, $
                     ;title = 'file  = '+filename+' pix = '+ $
                     ; string(i)+string(j)+' gratpos = '+string(k), color = 200
                 ;endif else if ( j ne 0 and j ne 17) then  begin
                 ;   oplot, raw_data2[i, j, *, kk], psym = 4, color = 200
                 ;endif
                                  
                 rmax = max(raw_data2[i, j, *, kk], locmax)
                 
                 if locmax lt ramplength-3 then begin
                    q_end = locmax-1
                                      
                 endif else q_end = ramplength-3

                 if q_end gt 0 then begin
                    ;y = A+Bx
                    fit[0:1, i, j, kk] = linfit(x[0:q_end], $
                                                raw_data2[i, j, 0:q_end, kk], $
                                                /double, chisq = chisq, $
                                                sigma = sigma)
                    
                    ; sigma for slope
                    fit[2, i, j, kk] = sigma[1] 

                    ;if (j ne 0 and j ne 17) then begin
                       ; oplot, fit[0, i, j, kk]+fit[1, i, j, kk]*$
                       ; indgen(ramplength), color = 200
                    ; endif
  
                    
                    if (s2n gt 0) then begin
                       if fit[1, i, j, kk] / fit[2, i, j, kk] lt s2n then begin
                          fit[1, i, j, kk] = !VALUES.F_NAN

                       endif  
                    endif
                    
                 endif else begin
                    fit[1, i, j, kk] = !VALUES.F_NAN
                 endelse
              endfor
              
              if verbose then print, 'slope = ', fit[1, i, j, *]
              
              robustmeancomb, reform(fit[1, i, j, *], n_elements(fit[1, i, j, *])), 2, mn, mvar, DATAVAR = reform(fit[2, i, j, *], n_elements(fit[2, i, j, *]))
              
              flux[i, j] = mn
              stddev_arr[i, j] = mvar
              
              
           endfor
        endfor
     endelse


     flux_cube = reform(flux, 5, 5, 18)
     stddev_cube = reform(stddev_arr, 5, 5, 18)
     fits_struct = create_struct('header', exthdr, 'data', flux_cube, $
                                 'stddev', stddev_cube)
     mwrfits, fits_struct, fitfilename, /silent

     k += 1
  endwhile

  return, fits_struct
end
