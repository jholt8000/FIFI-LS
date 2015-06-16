;+
; NAME: fifi_ls_subtract_chops 
;
; PURPOSE: Subtract_Chops chops of fitted data (written fifi_ls_fit_ramps.pro)
 
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: result=subtract_chops(FILEROOT, verbose=verbose) 
;
; INPUTS: FILEROOT - root name of the input file (*.chop#.rampfits.fits)
;
; OPTIONAL INPUTS:
;
; OUTPUTS: 1. 
;          2. fits file containing averaged, subtract_chopsd ramp slope 
;             (i.e. flux in arbitrary units) and sigma.
;             File name: *.demod2pt.L2.fits or *.demod4pt1.L2.fits
;             fits file contains:
;                      prime header
;                      extensions for each grating with the following:
;                      ext.header = extension header
;                      ext.data = chop subtracted flux cube, shape=(5,5,18)
;                      ext.stddev = err cube, shape=(5,5,18)
;              
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS: 
; 
; RESTRICTIONS: CAUTION on and off might be switched; results may need to be
;               multiplied by -1 for them to make sense.
;
; PROCEDURE: 1. Read input file; use the ramp count in the header to identify
;               the chop cycle (file does not need to start with a complete 
;               chop cycle).
;            2. Subtract OFF from ON.
;            3. Sum of sigma squares
;            4.  write results to FITS extension.
; 
; EXAMPLE: result=subtract_chops('00003_rawlw031')
;          fluxer, result.data
;
; MODIFICATION HISTORY:
; 2015Mar22  jholt
;-
function fifi_ls_subtract_chops, fileroot, verbose = verbose

;--------------
;  help file
;--------------
if n_params() eq 0 then begin
  doc_library, 'fifi_ls_subtract_chops'
  return, -1
endif

filenames = file_search(fileroot+'*rampfits.fits')

; there should be a chop0 and a chop1 file to subtract
if (n_elements(filenames) ne 2) then begin
   print, ' I cannot find only a chop0 and a chop1 file'
   return, -1
endif

fits_info, filenames[0], N_ext = n_ext0, /silent
fits_info, filenames[1], N_ext = n_ext1, /silent

; Only handle same number of grating steps in each chop (should make
; this better later) 
if n_ext0 ne n_ext1 then begin
   print, ' the # of inductosyn postions is not the same between chop files'
   return, -1
endif

; open the first file primary header
primehead  =  headfits(filenames[0])

outfilename1 = strmid(filenames[0], 0, $
                      strlen(filenames[0]) - 4) + 'chopsub.fits'
outfilename1 = str_replace(outfilename1, 'chop0.', '')
outfilename1 = str_replace(outfilename1, 'chop1.', '')
outfilename1 = str_replace(outfilename1, 'rampfits.', '')

; if the output file already exists, delete it
if file_test(outfilename1) gt 0 then file_delete, outfilename1

fxwrite, outfilename1, primehead, /noupdate, /append

for i = 1, n_ext0 do begin
   
   f0 = mrdfits(filenames[0], i, /unsigned, /silent)
   f1 = mrdfits(filenames[1], i, /unsigned, /silent)

   hdr0 = f0.header
   hdr1 = f1.header

   data0 = f0.data
   data1 = f1.data

   stddev0 = f0.stddev
   stddev1 = f1.stddev

   nodpos0 = fxpar(primehead, 'NODBEAM')
   nodstyle = strlowcase(strtrim(fxpar(primehead, 'NODSTYLE')))
   multiplier = 1.
   if (nodstyle eq 'symmetric') and (strtrim(nodpos0, 2) eq 'B') then $
      multiplier = -1. 

   ; subtract chop0-1 for A and chop1-0 for B
   
   if (fxpar(hdr0, 'INDPOS')) eq (fxpar(hdr1, 'INDPOS')) then begin
      if (fxpar(hdr0, 'CHOPNUM') lt 1) then begin
        data = multiplier * (data0 - data1)
     endif else begin
        data = multiplier * (data1 - data0)
     endelse
     stddev = sqrt(stddev0^2 + stddev1^2) 
  endif else begin
     print, 'inductosyn positions do not line up'

  endelse

  fits_struct = create_struct('header', hdr0, 'data', data, 'stddev', stddev ) 
  mwrfits, fits_struct, outfilename1, /silent

endfor

return, fits_struct

end
