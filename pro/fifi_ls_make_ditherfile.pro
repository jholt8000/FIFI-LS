;+
; NAME: fifi_ls_make_ditherfile
;
; PURPOSE: Combine nods of ramp-fitted, chop-subtracted data 
         
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: result=fifi_ls_make_ditherfile(FILEROOT, verbose=verbose) 
;
; INPUTS: FILEROOT - root name of the input files to be combined 
;
; OPTIONAL INPUTS:
;
; OUTPUTS: 1. returns IDL struct of last nod + grating postion 
;          2. FITS file fileroot+'nodcomb.fits' containing an IDL structure for
;             each grating postion in its own FITS extension:
;             flux data in x.data componenent in a size 5,5,16 array
;             noise data in x.sttdev in a size 5,5,16 array
;             extension header in x.header
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
; PROCEDURE: 1. Read header information from each of the input files found
;            2. Find unique combinations of channel+dither+inductosyn postions
;            3. Add the 'A' nods and subtract the 'B' nods
;            4. Create fits file and write results to disk.
; 
; EXAMPLE: foo=fifi_ls_make_ditherfile('00003_rawlw031')
;          fluxer, foo.data
;
; MODIFICATION HISTORY:
; 2015Mar22  jholt:  Created new routine to handle separate inductosyn
;                    extensions
; 
;-
function fifi_ls_make_ditherfile, fileinfo, verbose
;--------------
;  help file
;--------------
  not_comb = []
  
if n_params() eq 0 then begin
  doc_library,'fifi_ls_subtract_chops'
  return, -1
endif

if ~keyword_set(verbose) then verbose=0

if (n_elements(fileinfo) eq 1) then begin
   fileroot=fileinfo
   filenames=file_search('*'+fileroot+'*.fits')
endif else begin
   filenames = fileinfo
   fileroot = strmid(fileinfo[0], 19, 28)
endelse

if file_test('dither_info.txt') gt 0 then file_delete, 'dither_info.txt'


openw, lun, 'dither_info.txt', /get_lun, /append
printf, lun, 'filename', 'chan','nodbea','nodsty', 'telra', 'teldec', 'obslam',  'obsbet', 'dlammap', 'dbetmap',$
        format='(a49,2x,a3,2x,a7,2x,a7,3x,a10,3x,a10,3x,a10,3x,a10, 2x,a7, 2x, a7)'
i=0
while (i lt n_elements(filenames)) do begin
   
   primehead = headfits(filenames[i])

   channel = strtrim(fxpar(primehead, 'CHANNEL'))
   if channel eq 'BLUE' then channel=1 else channel=0

   telra = fxpar(primehead, 'TELRA')
   teldec = fxpar(primehead, 'TELDEC')
   ; read dither positions
   obslam =  fxpar(primehead, 'OBSLAM')
   obsbet = fxpar(primehead, 'OBSBET')
   
   dlam_map = fxpar(primehead, 'DLAM_MAP')
   dbet_map = fxpar(primehead, 'DBET_MAP')

   ; read nod position
   nod_beam = strtrim(fxpar(primehead, 'NODBEAM'))

   nod_style = strtrim(fxpar(primehead, 'NODSTYLE'))
    
   printf, lun, filenames[i], channel, nod_beam, nod_style, telra, teldec, obslam,  obsbet, dlam_map, dbet_map, $
            format='(a49,2x,i3,2x,a7,2x,a7,3x,f10.6,3x,f10.6,3x,f10.5,3x,f10.5,3x, i7,i7)'
   
  
    i+=1
    
endwhile

free_lun, lun

return, obsbet

end
