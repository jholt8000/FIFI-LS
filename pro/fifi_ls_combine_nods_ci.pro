;+
; NAME: fifi_ls_combine_nods_ci
;
; PURPOSE: Combine nods of ramp-fitted, chop-subtracted data 
         
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: result=fifi_ls_combine_nods_ci(files)
;
; INPUTS: files - input files, must have matching A's and B's positions
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
;            2. ## removed, do all your checks by hand and only put
;            files to be combined!!!  Find unique combinations of
;            channel+dither+inductosyn postions
;            3. Add the 'A' nods and subtract the 'B' nods
;            4. Create fits file and write results to disk.
; 
; EXAMPLE: foo=fifi_ls_combine_nods_for_chris(infiles)
;          fluxer, foo.data
;
; MODIFICATION HISTORY:
; 2015Jun09  jholt:  Created new routine to work under IDL 7.0,
; WARNING! this is toilet code. ;)
; 
;-
function fifi_ls_combine_nods_ci, fileinfo, verbose

  ; if symmetric chopping then use dlam_map and dbet_map ->
  ; they should match (along with TELRA/TELDEC/INDPOS/CHANNEL)
  ; if asymmetric chopping then use dlam_map and dbet_map only for As
  ; use closest in time B nods, what about TELRA and TELDEC?

  
if n_params() eq 0 then begin
  doc_library,'fifi_ls_subtract_chops'
  return, -1
endif

 not_comb = []
if ~keyword_set(verbose) then verbose = 0

if (n_elements(fileinfo) eq 1) then begin
   fileroot=fileinfo
   if fileroot eq 'all' then begin
      filenames = file_search('*chopsub*.fits')
      chp_pos = strpos(filenames[0], 'chopsub')
      fileroot = strmid(filenames[0], chp_pos-30, chp_pos-20)
      fileroot = str_replace(fileroot, 'chopsub.fits','nodcomb.fits')
   endif else begin
      filenames=file_search('*'+fileroot+'*chopsub*.fits')
   endelse
   
endif else begin
   filenames = fileinfo
   fileroot = strmid(fileinfo[0], 19, 28)
endelse

if verbose then print,'filenames=',filenames

; create empty list, to be added to
;datalist = list()

;if (write_same_position_to_one) then begin
;   outfile = str_replace(filenames[0],'.chopsub.fits', '.nodcomb.fits')
;   if file_test(outfile) then file_delete, outfile
;endif

;------------------------------------------------------
; set up two master lists of all the information needed
; one array containing only items that need to match, this will have
; duplicates removed later
; the other list contains those match items and all associated metadata and data
; to combine ## 
;------------------------------------------------------

i = 0
while (i lt n_elements(filenames)) do begin
   
   primehead = headfits(filenames[i], ext=0)

   nodstyle = strlowcase(strtrim(fxpar(primehead, 'NODSTYLE')))
   if (nodstyle eq 'asymmetric') then begin
      result = fifi_ls_combine_nods_asym(filenames, verbose)
      return, result
   endif
   
   nodpatt = strtrim(fxpar(primehead, 'NODPATT'))
   if (nodpatt ne 'ABBA') then begin
      print,'WARNING: nodpattern is ', nodpatt,$
            ' I might not be combining nods correctly!'
   endif
   
   channel = strtrim(fxpar(primehead, 'CHANNEL'))
   if channel eq 'BLUE' then channel = 1 else channel = 0

   ; read dither positions
   dlam_map = fxpar(primehead, 'DLAM_MAP')
   dbet_map = fxpar(primehead, 'DBET_MAP')

   if verbose then print,'fn, dlam, dbet = ',filenames[i], dlam_map, dbet_map
   ; read nod position
   nod_beam = strtrim(fxpar(primehead, 'NODBEAM'))

   ; determine the number of extensions
   ; (and hence the number of grating postions)
   fits_info,filenames[i], n_ext = n_ext, /silent
 
   k = 1
   while (k lt n_ext+1) do begin
      fits_struct = mrdfits(filenames[i], k, /silent)
      data = fits_struct.data
      exthdr = fits_struct.header
      sigma = fits_struct.stddev
      ; read inductosyn postion
      ind_pos = fxpar(exthdr,'INDPOS')

      ; this is the array of keywords that need to match 
      ; in each combined nod pos
      need2match = [channel, dlam_map, dbet_map, ind_pos]
      
      ; create master match array, if it isn't already there
      ; otherwise append to it
      if (i eq 0) and (k eq 1) then need2match_big = [need2match] else $
         need2match_big = [[need2match_big],[need2match]]

      ; this is the list that has the need-to-match information as well as all
      ; other relevant metadata and data needed to combine nods
      ;datalist.Add, list(filenames[i], channel, dlam_map, dbet_map, ind_pos,$
      ;                   nod_beam, data, sigma, primehead, exthdr)
      k += 1
   endwhile
 
   i += 1
endwhile

; find only unique (grating postions + dither position + channel) combinations
; need2match_big is an array with rows: channel, dlam, dbet, indpos
if verbose then print,'need2match_big = ',need2match_big

size_n2m = size(need2match_big)

;return, need2match_big
; if there is more than one row
if (size_n2m[0] gt 1) then begin
    ; first sort on inductosyn number (grating position)
    n2m = need2match_big[*,sort(need2match_big[3,*])]

    if verbose then print,'n2m = ', n2m
    rn = 10.7*n2m[0,*] + 14.2*n2m[1,*] + 16.3 * n2m[2,*] + n2m[3,*] / 23.

    uniq_indices = uniq(rn, sort(rn))

    uniq_n2m = n2m[*,uniq_indices]
    
    s = size(uniq_n2m)
    ; number of unique rows
    if s[0] gt 1 then n_uniq_rows = s[2] $
    else n_uniq_rows = 1
 
 endif else begin
    uniq_n2m = need2match_big
    n_uniq_rows = 1

 endelse

if verbose then print, 'need2match_uniq = ',uniq_n2m

ii = 0
; ii is looping over uniq channel+dither+grat combos
while ii lt n_uniq_rows do begin
   print,'----------------'
   uniq_row = uniq_n2m[*, ii]
   u_channel = uniq_row[0]
   u_dlam = uniq_row[1]
   u_dbet = uniq_row[2]
   u_indpos = uniq_row[3]
   
   comb_data = dblarr(5, 5, 18)
   stddev_data_sqd = dblarr(5, 5, 18)

   combined = 0
   files_combed = ''
   i = 0
  
   while i lt n_elements(filenames) do begin
      dl_primehead = headfits(filenames[i])
      
      fits_info,filenames[i], n_ext = n_ext, /silent
           
      dl_filename = filenames[i]
      file_for_this = dl_filename

      ; read dither positions
      dl_dlam = fxpar(dl_primehead, 'DLAM_MAP')
      dl_dbet = fxpar(dl_primehead, 'DBET_MAP')

      ; read nod position
      dl_nodbeam = strtrim(fxpar(dl_primehead, 'NODBEAM'))

      dl_channel = strtrim(fxpar(dl_primehead, 'CHANNEL'))
      if dl_channel eq 'BLUE' then dl_channel = 1 else dl_channel = 0

      k = 1
      while (k lt n_ext+1) do begin
         fits_struct = mrdfits(filenames[i], k, /silent)
         dl_data = fits_struct.data
         dl_exthdr = fits_struct.header
         dl_sigma = fits_struct.stddev
         ; read inductosyn postion
         dl_indpos = fxpar(dl_exthdr,'INDPOS')


         print, 'ii=',ii,' i=',i,' k=',k, dl_indpos
         print,'working on uniq row=', uniq_row
         print,'file info=',dl_dlam, dl_dbet, dl_indpos
        
         if ((u_channel eq dl_channel) and (u_dlam eq dl_dlam) and $
            (u_dbet eq dl_dbet) and (u_indpos eq dl_indpos)) then begin

            print,'here'
            exthdr = dl_exthdr
            primehead = dl_primehead
            file_dlam = strtrim(fix(dl_dlam),2)
            file_dbet = strtrim(fix(dl_dbet),2)

            files_combed = files_combed+' '+dl_filename
            ;math handled in subtract chops -jh
            ;if dl_nodbeam eq 'B' then begin
            ;   comb_data = comb_data - dl_data
            ;endif else if dl_nodbeam eq 'A' then begin
            comb_data = comb_data + dl_data
         
            stddev_data_sqd = stddev_data_sqd + dl_sigma^2 
            combined = combined + 1
            print,'combined=',combined
         endif
         k += 1
      endwhile
      
      i += 1
   endwhile

   ; if more than one file were actually combined, write output
   if (combined gt 1) then begin
      if u_channel eq 1 then file_chan = '_sw' $
      else file_chan = '_lw'

      ; take sqrt of the sum of the squares
      stddev_data = sqrt(stddev_data_sqd)

      fxaddpar, exthdr,  'NODFILES', files_combed, 'Files nod combined'

      outfilestr = strsplit(files_combed, /extract)
      outfileroot = outfilestr[0]
      ;fileroot = str_replace(outfileroot,'.chopsub.fits', '')
      outfile = fileroot+'_'+file_dlam+'_'+file_dbet+file_chan+'.nodcomb.fits'
   
      
      fits_struct = create_struct('header', exthdr, 'data', comb_data, $
                                  'stddev', stddev_data)

      if ~file_test(outfile) then fxwrite, outfile, primehead , $
                                           /noupdate, /append
      ; write the inductosyn extension for this particular unique
      ; channel+dither position+indpos
      mwrfits, fits_struct, outfile, /silent
      ;not_comb = [0]
   endif else begin
      if verbose then begin
          print, 'only one of the following channel/dither/inductosyn found:'
          print, uniq_row
      endif

      if (n_elements(not_comb) lt 1) then begin
         not_comb = [file_for_this]
      endif else begin
         not_comb = [file_for_this, not_comb]
      endelse
            
   endelse
   ii += 1
endwhile

if (n_elements(not_comb) gt 0) then begin
    not_comb_u = not_comb[uniq(not_comb, sort(not_comb))]
    openw, lun, 'combine_nods_not_combined.txt', /get_lun
    printf, lun, not_comb_u
    free_lun, lun

endif

return, fits_struct

end
