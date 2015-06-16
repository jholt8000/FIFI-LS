;+
; NAME: fifi_ls_spatial_calibrate 
;
; PURPOSE: Apply spatial calibration (offsets)
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=new_calibrate(FILEROOT,OBSDATE=OBSDATE,verbose=verbose)
;
; INPUTS: FILEROOT - root name to apply spectral/spatial calibration
;
; OPTIONAL INPUTS:
;
; KEYWORDS: OBSDATE - Date of observation.  Intended for files that do not 
;                     have the DATE-OBS keyword (and value) in the fits 
;                     primary header (earlier files do not).  If used, input
;                     is [YYYY,MM,DD]
;
; OUTPUTS: RESULT - 1. 1 if successful, -1 if failure
;                   2. fits file containing [lambda,x,y,module#,flux,sigma] 
;                      all 400 pixels (a 6x400 array) as a binary table. 
;                      File name: *.cal.fits
;                      fits content:
;                       prime header
;                       binary table - (header + [6,400] data array) x
;                                      number of unique grating positions 
;
; OPTIONAL OUTPUTS: 
;
; SIDE EFFECTS: 1. The two open/dummy "spectral" channels are removed.
;               2. Bad data (bad/dead pixel, questionable data quality) are
;                  flagged and the wavelength (and only wavelength) is set
;                  to NaN. 
;
; RESTRICTIONS: Input file must be a demodulated file (FIFILS fits keyword
;               DEMOD) but not necessarily has to be normalized.  Demodulated
;               files may contain multiple grating positions, but the are all
;               unique (i.e. no repeat Inductosyn Positions).
;
; PROCEDURE: 1. Read input file; discard dummy/open spectral channels.
;            2. Convert spectral channel to wavelength.
;            3. Convert spatial modules to position (coord. sys TBD) 
;            4. flag bad/dead pixels, bad data by setting wavelength (and 
;               only wavelength) to NaN 
;            5. Create fits file and write results to disk.
;
; EXAMPLE: result=new_calibrate('00121_221031_FullScan-G_STRT_R-752500_lw')
;
; MODIFICATION HISTORY:
;      J.Holt 2015Mar22 Split original calibrate script by K. and
;      R. Klein into separate lambda and xy calibration
;      scripts. Updated to handle only single inductosyn position
;      files.
 
; $Id: fifi_ls_spatial_calibrate,v 1.20 2014-03-24 23:51:18 klein Exp $
;-
function fifi_ls_spatial_calibrate, filename, obsdate=obsdate, verbose=verbose

   ;  help file
   if n_params() eq 0 then begin
     doc_library, 'fifi_ls_spatial_calibrate'
     return, -1
   endif
   
   if file_test(filename) eq 0 then begin
      print, 'cannot find file ', filename
      return, -1
   endif

   primehead = headfits(filename, ext=0)

   channel=strtrim(fxpar(primehead, 'CHANNEL'))
   b_order=fxpar(primehead, 'G_ORD_B')
   
   fits_info, filename, n_ext=n_ext, /silent

   outfile = str_replace(filename, '.fits', '')
   outfile = outfile + '_xy.fits'

   if (file_test(outfile) gt 0) then file_delete, outfile
   
   k=1
   while (k lt n_ext+1) do begin
      
      fits_struct = mrdfits(filename, k, /unsigned, /silent) 
      
      data = fits_struct.data
      stddev = fits_struct.stddev
      exthdr = fits_struct.header
      lambda = fits_struct.lambda

      ; get date of observation to determine which lookup table to use
      if ~keyword_set(OBSDATE) then begin
         ; in YYYY-MM-DDThh:mm:ss[.sss] format
         temp=strtrim(fxpar(primehead, 'DATE-OBS'), 2) 
         if temp eq '0' then begin
            print, 'fifi_ls_spatial_calibrate: OBSDATE not found in file header'
            return, -1
         endif
         ; [yyyy, mm, dd]
         obsdate=fix(strsplit(strmid(temp, 0, 10), '-', /extract))	
      endif

      if (channel eq 'BLUE') then begin
         b_order=fxpar(primehead, 'G_ORD_B')
         if b_order eq 1 then blue = 'B1' else $
            if b_order eq 2 then blue = 'B2' else begin
            print, 'Invalid Blue grating order.'
            return, -1
         endelse
         xy = fifi_ls_offset_xy(obsdate, /blue)
         G_STRT = fxpar(primehead, 'G_STRT_B') 

      endif else begin
         xy=fifi_ls_offset_xy(obsdate)
         G_STRT= fxpar(primehead, 'G_STRT_R')
         
      endelse

      temp=fltarr(2, 400)
      ; repeat xy values every 25th index
      for i=0, 400/25-1 do temp[0, i*25]=xy

      ; (x, y) positions 
      if (strtrim(fxpar(primehead, 'OBJ_NAME'), 2) eq 'TELSIM') or $
         (strtrim(fxpar(primehead, 'OBJECT'), 2) eq 'TELSIM')then begin
         ; x-y stage, mm
         delta_x = fxpar(primehead, 'DLAM_MAP')
         delta_y = fxpar(primehead, 'DBET_MAP')
         if !err eq -1 then begin
            ;  no DLAM/BET_MAP, ie. old header 
            delta_x = fxpar(primehead, 'OBSDEC')/10. ; OBSRA  
            delta_y = fxpar(primehead, 'OBSRA')/10.  ;OBSDEC
         endif
         xs = -temp[1, *]+ delta_x
         ys = -temp[0, *]+ delta_y
         
      endif else begin

         DA      = fxpar(primehead, 'DET_ANGL') * !DTOR ; ccw array rotation
         DLAM_MAP= fxpar(primehead, 'DLAM_MAP') ; Map offset in arcsec           
         DBET_MAP= fxpar(primehead, 'DBET_MAP') ; Map offset in arcsec    
         PLATSCAL= fxpar(primehead, 'PLATSCAL') ; Plate scale in "/mm
         PRIMARAY= strtrim(fxpar(primehead, 'PRIMARAY'))
         DICHROIC= fxpar(primehead, 'DICHROIC')
         INDPOS = fxpar(exthdr, 'INDPOS')
         
         if ((channel eq 'RED') && (PRIMARAY eq 'RED')) || $
            ((channel eq 'BLUE')  && (PRIMARAY eq 'BLUE')) then begin
            Dx_RB = 0
            Dy_RB = 0
         endif else begin
            case DICHROIC of
               105: begin 
                  ;Dx_RB = -5.9365e-7 * G_STRT - 1.5010
                  ;Dy_RB = -4.5569e-8 * G_STRT + 0.2429
                  Dx_RB = -5.9365e-7 * INDPOS - 1.5010
                  Dy_RB = -4.5569e-8 * INDPOS + 0.2429
               end
               130: begin
                  ;Dx_RB = -8.4225e-7 * G_STRT - 0.7083
                  ;Dy_RB = -9.0334e-8 * G_STRT - 0.0807
                  Dx_RB = -8.4225e-7 * INDPOS - 0.7083
                  Dy_RB = -9.0334e-8 * INDPOS - 0.0807                 
               end
            endcase
            if (channel eq 'BLUE') then begin
               Dx_RB = -Dx_RB
               Dy_RB = -Dy_RB
            endif
         endelse
         
         xs = PLATSCAL* (temp[0, *] + Dx_RB) + cos(DA)*DLAM_MAP - $
              sin(DA)*DBET_MAP
         ys = PLATSCAL* (temp[1, *] + Dy_RB) + sin(DA)*DLAM_MAP + $
              cos(DA)*DBET_MAP

         ; recall that temp is just the poscal xy repeated every 25
         
      endelse  ; real target, not TELSIL
               ; module #;output.data[3, *]=indgen(400) mod 25 +1; module number (1-25)

;--------------------------------------------------------
; Flip sign for off-source (C_BEAM!=1) flux   20feb14
;what is this? ; -jh 2015
;--------------------------------------------------------
; off-source beam has the "source" position reversed
      C_BEAM=strtrim(fxpar(primehead, 'C_BEAM'))
      C_TIP =strtrim(fxpar(primehead, 'C_TIP'))
      if C_TIP ne 0 then message, /cont, 'Only symmetric chopping supported'
; C_BEAM = 1: nod A
; C_BEAM = 0: C_TIP must be 1 (asymmetric chopping), nod B
; C_BEAM = -1: nod B
      ;if C_BEAM ne 1 then print, ' was multiplying by -1'
;if C_BEAM ne 1 then output.data[4, *] = -1. * output.data[4, *]

      fits_struct = create_struct('header', exthdr, 'data', data, 'lambda', $
                                  lambda, 'stddev', stddev, 'xs', xs, 'ys', ys)

      if ~file_test(outfile) then fxwrite, outfile, primehead, /noupdate, $
                                           /append
      mwrfits, fits_struct, outfile, /silent
      k+=1

   endwhile
   
   return, fits_struct

end
