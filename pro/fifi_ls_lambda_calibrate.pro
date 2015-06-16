;+
; NAME: fifi_ls_lambda_calibrate 
;
; PURPOSE: Apply spectral calibration 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=fifi_ls_lambda_calibrate(FILEROOT,
;                          OBSDATE=OBSDATE, verbose=verbose)
;
; INPUTS: FILEROOT - root filename to apply spectral calibration
;
; OPTIONAL INPUTS:
;
; KEYWORDS: OBSDATE - Date of observation.  Intended for files that do not 
;                     have the DATE-OBS keyword (and value) in the fits 
;                     primary header (earlier files do not).  If used, input
;                     is [YYYY,MM,DD]
;           verbose - Not used
;           
;
; OUTPUTS: RESULT - 1. Returns last IDL struct in list of files and
;                   inductosyn postions. 
;                   
;                   2. FITS file fileroot+'lambda.fits'
;                   containing an IDL structure for
;                   each grating postion in its own FITS extension:
;                   flux data in x.data componenent in a size 5,5,16 array
;                   noise data in x.sttdev in a size 5,5,16 array
;                   lambda data in x.lambda in a size 5,5,16 array
;                   extension header in x.header
;
; OPTIONAL OUTPUTS: 
;
;
; RESTRICTIONS: Input file must be a single chop
;
; PROCEDURE: 1. Read input file; discard dummy/open spectral channels.
;            2. Convert spectral channel (ind pos) to wavelength.
;            3. Create fits file and write results to disk.
;
; EXAMPLE:
;          result=fifi_ls_lambda_calibrate('00121_221031_FullScan-G_STRT_R-752500_lw')
;          print, result.lambda
;
;
; MODIFICATION HISTORY:
; 2015Mar22 J.Holt: Split original calibrate script by K. Nishikida and
;                    R. Klein into separate lambda and xy calibration
;                    scripts. Updated to handle single inductosyn
;                    position extensions.
;-
function fifi_ls_lambda_calibrate, filename, obsdate = obsdate, $
                                   verbose = verbose

   ;  help file
   if n_params() eq 0 then begin
     doc_library, 'fifi_ls_lambda_calibrate'
     return, -1
   endif

   if file_test(filename) eq 0 then begin
      print, 'cannot find file ', filename
      return, -1
   endif
   
   primehead = headfits(filename, ext = 0)

   channel = strtrim(fxpar(primehead, 'CHANNEL'))
   b_order = fxpar(primehead, 'G_ORD_B')
   
   fits_info, filename, n_ext = n_ext, /silent

   outfile = str_replace(filename, '.fits', '.lambda.fits')

   if (file_test(outfile) gt 0) then file_delete, outfile
   
   k = 1
   while (k lt n_ext+1) do begin

      fits_struct = mrdfits(filename, k, /unsigned, /silent) 

      data = fits_struct.data
      stddev = fits_struct.stddev
      exthdr = fits_struct.header
      
      ; remove the first and last spectral channels (0 and 17)
      data = data[*, *, 1:16]
      stddev = stddev[*, *, 1:16]
      ;-----------------------------------
      ; spectral  calibration
      ;-----------------------------------
      ; get date of observation
      if ~keyword_set(OBSDATE) then begin
         
         temp = strtrim(fxpar(primehead, 'DATE-OBS'), 2)
         ; in YYYY-MM-DDThh:mm:ss[.sss] format
         
         if temp eq '0' then begin
            ; this is a hack to get the test data working -jh
            temp = '2015-03-02T12:12:12.12'
            print, 'fifi_ls_lambda_calibrate: OBSDATE not in file header'
            ;return, -1
         endif 
         obsdate = fix(strsplit(strmid(temp, 0, 10), '-', /extract)) 
      endif

      ; get inductosyn position
      ind = fxpar(exthdr, 'INDPOS')

      if (channel eq 'BLUE') then begin
         if b_order eq 1 then blue = 'B1' else $
            if b_order eq 2 then blue = 'B2' else begin
            print, 'Invalid Blue grating order.'
            return, -1
         endelse
         
         lambda = fifi_ls_new_wave(ind, obsdate, blue = blue)    
      endif else begin
         
         lambda = fifi_ls_new_wave(ind, obsdate)
      endelse

      if lambda[0] eq -1 then begin
         print, 'fifi_ls_lambda_calibrate: wavelength calibration failed.'
         return, -1
      endif
      lambda = reform(lambda, 5, 5, 16)
      
      ;----------------------------
      ; flag bad pixels, bad data -this needs to be elsewhere JH
      ;----------------------------
      ; read static dead pixel file	[module#, channel#]
      if (channel eq 'RED') then  badpix = fifi_ls_bad_pixels(obsdate) else $
         badpix = fifi_ls_bad_pixels(obsdate, /blue)

      ; convert [module#, channel#] to pixel # (in the 0 though 399 convention) 
      flagpix = (badpix[0, *]-1)+(badpix[1, *]-1) * 25 ; assuming modules 1-25 
      
      ; set wavelengths of bad pixels to NaN.
      data[flagpix] = !values.f_nan
      
      ; should these go in primehead? -jh
      sxaddpar, exthdr, 'CRPIX1', 2.5
      sxaddpar, exthdr, 'CRPIX2', 2.5
      sxaddpar, exthdr, 'CRPIX3', 0
      sxaddpar, exthdr, 'CRVAL1', sxpar(exthdr, 'OBSLAM')
      sxaddpar, exthdr, 'CRVAL2', sxpar(exthdr, 'OBSBET')
      sxaddpar, exthdr, 'CRVAL3', lambda[0]
      sxaddpar, exthdr, 'CDELT3', abs(lambda[25]-lambda[0])
      ;constant delta lambda in this case
      sxaddpar, exthdr, 'CDELT1', 6./3600.   ; deg/pix
      sxaddpar, exthdr, 'CDELT2', 6./3600. ; deg/pix
      
      ;----------------------------------------------------------
      ; update fits header, write grating position to output file
      ;----------------------------------------------------------
      fits_struct = create_struct('header', exthdr, 'data', data, $
                                  'lambda', lambda, 'stddev', stddev)

      if ~file_test(outfile) then fxwrite, outfile, $
                                           primehead, /noupdate, /append
      mwrfits, fits_struct, outfile, /silent
      k += 1

   endwhile

   ;ds9file = str_replace(outfile, '.fits', '.fluxonly.fits')
   ;fxwrite, ds9file, primehead
   ;mwrfits, data, ds9file, /silent
   
   return, fits_struct
   

end
