;+
; NAME: fifi_ls_combine_nods
;
; PURPOSE: Combine nods of ramp-fitted, chop-subtracted data 
         
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: result=fifi_ls_combine_nods(FILEROOT, verbose=verbose) 
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
; EXAMPLE: foo=fifi_ls_combine_nods('00003_rawlw031')
;          fluxer, foo.data
;
; MODIFICATION HISTORY:
; 2015Mar22  jholt:  Created new routine to handle separate inductosyn
;                    extensions
; 
;-
function fifi_ls_combine_nods_asym, fileinfo, verbose, $
                               write_same_position_to_one


  ; need to go through and for each 'A', add nearby matching 'A's
  ; find the closest in time 'B' with matching nodpos/channel

  not_comb = []
  
if ~keyword_set(verbose) then verbose=0

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

; sort the filenames
filenames = filenames[sort(filenames)]

datalistA = list()
datalistB = list()
meantimesA = list()
meantimesB = list()
Btimes = list()
;------------------------------
i = 0
while (i lt n_elements(filenames)) do begin
   file = filenames[i]
  
   primehead = headfits(filenames[i], ext=0)

   nodstyle = strlowcase(strtrim(fxpar(primehead, 'NODSTYLE')))
   ;if (nodstyle eq 'asymmetric') then begin
   ;   result = fifi_ls_combine_nods_asym(filenames, verbose)
   ;   return, result
   ;endif
   
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

   ; read times
   utstart = strtrim(fxpar(primehead, 'UTCSTART'))
   utend = strtrim(fxpar(primehead, 'UTCEND'))
   uts = strsplit(utstart, ':', /extract)
   ute = strsplit(utend, ':', /extract)
   hours_start = uts[0] + uts[1]/60. + uts[2]/3600.
   hours_end = ute[0] + ute[1]/60. + ute[2]/3600.
   meantime = (hours_start + hours_end) / 2 

   ; read positions
   telra = strtrim(fxpar(primehead, 'TELRA'))
   teldec = strtrim(fxpar(primehead, 'TELDEC'))

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
      if (nod_beam eq 'A') then begin
          need2matchA = [channel, dlam_map, dbet_map, ind_pos]
      
         ; create master match array, if it isn't already there
         ; otherwise append to it
         if (i eq 0) and (k eq 1) then need2match_bigA = [need2matchA] else $
            need2match_bigA = [[need2match_bigA],[need2matchA]]      
         ; this is the list that has the need-to-match information and
         ; other relevant metadata and data needed to combine nods
         datalistA.Add, list(filenames[i], channel, dlam_map, dbet_map, $
                             ind_pos, data, sigma, primehead, $
                             exthdr, telra, teldec, meantime)
         meantimesA.Add, list(filenames[i], meantime)

      endif else begin
         
          if (i eq 0) and (k eq 1) then need2match_bigA = [need2matchA] else $
             need2match_big = [[need2match_bigA],[need2matchA]]
          
          datalistB.Add, list(filenames[i], channel, dlam_map, dbet_map, $
                              ind_pos, data, sigma, primehead, $
                                exthdr, telra, teldec, meantime)
          Btimes.Add, meantime
          ; this needs to be a pandas table or a sql database
                       
      endelse
      
      
      k += 1
   endwhile

   i += 1
   endwhile

   ; first find all the unique A/chan/nod/dlam/dbets

   size_n2m = size(need2match_bigA)

   if (size_n2m[0] gt 1) then begin
    ; first sort on inductosyn number (grating position)
    n2m = need2match_bigA[*,sort(need2match_bigA[3,*])]

    if verbose then print,'n2m = ', n2m
    rn = 10.7*n2m[0,*] + 14.2*n2m[1,*] + 16.3 * n2m[2,*] + n2m[3,*] / 23.

    uniq_indices = uniq(rn, sort(rn))

    uniq_n2m = n2m[*,uniq_indices]
    
    s = size(uniq_n2m)
    ; number of unique rows
    if s[0] gt 1 then n_uniq_rows = s[2] $
    else n_uniq_rows = 1
 
  endif else begin
    uniq_n2m = need2match_bigA
    n_uniq_rows = 1
    
  endelse

  ii = 0
  ; ii is looping over uniq channel+dither+grat combos for A nods
  while ii lt n_uniq_rows do begin
     uniq_row = uniq_n2m[*, ii]
     print,'uniq_row=',uniq_row
       u_channel = uniq_row[0]
       u_dlam = uniq_row[1]
       u_dbet = uniq_row[2]
       u_indpos = uniq_row[3]
   
       comb_data = dblarr(5, 5, 18)
       stddev_data_sqd = dblarr(5, 5, 18)

       combined = 0
       files_combed = ''
       ; go through A list and add, using just header matching
       kk = 0
       Bsubbed = list()
       while kk lt n_elements(datalistA) do begin 
           Arow = datalistA[kk]
           Afilename = Arow[0]
           Achannel = Arow[1]
           Adlam = Arow[2]
           Adbet = Arow[3]
           Aindpos = Arow[4]
           Adata = Arow[5]
           Asigma = Arow[6]
           Aprimehead = Arow[7]
           Aexthdr = Arow[8]
           Ara = Arow[9]
           Adec = Arow[10]
           Ameantime = Arow[11]
   
           if ((u_channel eq Achannel) and (u_indpos eq Aindpos)) and $
              (u_dlam eq Adlam) and  (u_dbet eq Adbet)  then begin

              exthdr = Aexthdr
              primehead = Aprimehead
              file_dlam = strtrim(fix(Adlam),2)
              file_dbet = strtrim(fix(Adbet),2)

              files_combed = files_combed+' '+Afilename

              comb_data = comb_data + Adata
              stddev_data_sqd = stddev_data_sqd + Asigma^2 
              combined = combined + 1
              file_for_this = Afilename
              print,'added file = ',file_for_this,' indpos=',Aindpos
              
           ;endif 
     
              nearestBtime_index = min(abs(Btimes.toArray() - Ameantime), idx)

             ; need to sort throught Btimes - Ameantime and find first match

              sorted_tdiff = sort(abs(Btimes.toArray() - Ameantime))
              print, abs(Btimes.toArray() - Ameantime)
              print, sorted_tdiff
              ;db_array = datalistB.toArray()
              s_dl_B = datalistB[sorted_tdiff]
              ll = 0
              while (ll lt n_elements(s_dl_B)) do begin
                 Brow = datalistB[ll]      
                 Bfilename = Brow[0]
                 Bchannel = Brow[1]
                 Bindpos = Brow[4]
                 Bdata = Brow[5]
                 Bsigma = Brow[6]
           
      
                 if ((u_channel eq Bchannel) and (u_indpos eq Bindpos)) then begin 
                 
                     files_combed = files_combed+' '+Bfilename
              
                     comb_data = comb_data - Bdata
                     stddev_data_sqd = stddev_data_sqd + Bsigma^2 
                     combined = combined + 1
                     file_for_this = Bfilename
                     print,'subtracted file = ',file_for_this,' indpos=',Bindpos
                     Bsubbed.Add, Bfilename
                 
                  endif else begin
                     print, 'need next closest B in time'
                     ;  bin = Value_Locate(Btimes.toArray(), Ameantime)
                  endelse
                  ll += 1
               endwhile
              
           endif
           kk += 1
        endwhile
       
        
       ; if more than one file were actually combined, write output
       if (combined gt 1) then begin
          if u_channel eq 1 then file_chan = '_sw' $
          else file_chan = '_lw'

          ; take sqrt of the sum of the squares
          stddev_data = sqrt(stddev_data_sqd)

          fxaddpar, exthdr,  'NODFILES', files_combed, 'Files nod combined'
      
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
              print, 'only one of the following chan/dither/inductosyn found:'
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
   
return, fits_struct
end
