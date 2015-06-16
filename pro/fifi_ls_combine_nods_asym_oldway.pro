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
function fifi_ls_combine_nods_asym_old_way, filenames, verbose


  ; need to go through and for each 'A', add nearby matching 'A's
  ; find the closest in time 'B' with matching nodpos/channel

  not_comb = []
  
if ~keyword_set(verbose) then verbose=0

; sort the filenames
filenames = filenames[sort(filenames)]

datalistA = list()
datalistB = list()

;------------------------------
i = 0
while (i lt n_elements(filenames)) do begin
   file = filenames[i]
  
   primehead = headfits(filenames[i], ext=0)

   nodstyle = strlowercase(strtrim(fxpar(primehead, 'NODSTYLE')))
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

   ; read times
   utstart = strtrim(fxpar(primehead, 'UTSTART'))
   utend = strtrim(fxpar(primead, 'UTEND'))
   uts = strsplit(utstart, ':', /extract)
   ute = strsplit(utend, ':', /extract)
   hours_start = uts[0] + uts[1]/60. + uts[2]/3600.
   hours_end = ute[0] + ute[1]/60. + ute[2]/3600.
   meantime = (hours_start + hours_end) / 2 

   ; read positions
   telra = strtrim(fxpar(primehead, 'TELRA'))
   teldec = strtrim(fxpar(primead, 'TELDEC'))

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
            need2match_big = [[need2match_bigA],[need2matchA]]      
         ; this is the list that has the need-to-match information and
         ; other relevant metadata and data needed to combine nods
         datalistA.Add, list(filenames[i], channel, dlam_map, dbet_map, $
                             ind_pos, nod_beam, data, sigma, primehead, $
                             exthdr, meantime, telra, teldec)
      endif else begin
          if (i eq 0) and (k eq 1) then need2match_bigA = [need2matchA] else $
            need2match_big = [[need2match_bigA],[need2matchA]]      
            datalistB.Add, list(filenames[i], channel, dlam_map, dbet_map, $
                              ind_pos, nod_beam, data, sigma, primehead, $
                                exthdr, meantime, telra, teldec)
      endelse
      
      
      k += 1
   endwhile
 
   i += 1
endwhile   


i = 0
while (i lt n_elements(filenames)-3) do begin

   file1 = filenames[i]
   file2 = filenames[i+1]
   file3 = filenames[i+2]
   file4 = filenames[i+3]

   ph1 = headfits(file1)
   ph2 = headfits(file2)
   ph3 = headfits(file3)
   ph4 = headfits(file4)

   ch1 = strtrim(fxpar(ph1, 'CHANNEL'))
   ch2 = strtrim(fxpar(ph2, 'CHANNEL'))
   ch3 = strtrim(fxpar(ph3, 'CHANNEL'))
   ch4 = strtrim(fxpar(ph4, 'CHANNEL'))
   
   nb1 = strtrim(fxpar(ph1, 'NODBEAM'))
   nb2 = strtrim(fxpar(ph2, 'NODBEAM'))
   nb3 = strtrim(fxpar(ph3, 'NODBEAM'))
   nb4 = strtrim(fxpar(ph4, 'NODBEAM'))

   telra1 = strtrim(fxpar(ph1, 'TELRA'))
   telra2 = strtrim(fxpar(ph2, 'TELRA'))
   telra3 = strtrim(fxpar(ph3, 'TELRA'))
   telra4 = strtrim(fxpar(ph4, 'TELRA'))
   teldec1 = strtrim(fxpar(ph1, 'TELDEC'))
   teldec2 = strtrim(fxpar(ph2, 'TELDEC'))
   teldec3 = strtrim(fxpar(ph3, 'TELDEC'))
   teldec4 = strtrim(fxpar(ph4, 'TELDEC'))

   utstart1 = strtrim(fxpar(ph1, 'UTSTART'))
   utstart2 = strtrim(fxpar(ph2, 'UTSTART'))
   utstart3 = strtrim(fxpar(ph3, 'UTSTART'))
   utstart4 = strtrim(fxpar(ph4, 'UTSTART'))
   utend1 = strtrim(fxpar(ph1, 'UTEND'))
   utend2 = strtrim(fxpar(ph2, 'UTEND'))
   utend3 = strtrim(fxpar(ph3, 'UTEND'))
   utend4 = strtrim(fxpar(ph4, 'UTEND'))
   
   primehead = ph1
   
   fits_info, file1, n_ext=n_ext1, /silent
   fits_info, file2, n_ext=n_ext2, /silent
   fits_info, file3, n_ext=n_ext3, /silent
   fits_info, file4, n_ext=n_ext4, /silent

   if (n_ext1 eq n_ext2) and (n_ext1 eq n_ext3) and $
      (n_ext3 eq n_ext4) then begin
      
      outfile = str_replace(file1, '.chopsub.fits', '.nodcomb.fits')
      outfile = str_replace(outfile,'_A_', '_')
      outfile = str_replace(outfile,'_B_', '_')
      
      if file_test(outfile) then file_delete, outfile
      outfile1 = str_replace(file1, '.chopsub.fits', '.nodcomb.fits')
      outfile2 = str_replace(file3, '.chopsub.fits', '.nodcomb.fits')
      outfile1 = str_replace(outfile1, '_A_', '_')
      outfile1 = str_replace(outfile1, '_B_', '_')
      outfile2 = str_replace(outfile2, '_A_', '_')
      outfile2 = str_replace(outfile2, '_B_', '_')

      if file_test(outfile1) then file_delete, outfile1
      if file_test(outfile2) then file_delete, outfile2
       
      k=1
      while (k lt n_ext1+1) do begin
          fs1 = mrdfits(file1, k, /silent)
          data1 = fs1.data
          exthdr1 = fs1.header
          stddev1 = fs1.stddev
          fs2 = mrdfits(file2, k, /silent)
          data2 = fs2.data
          exthdr2 = fs2.header
          stddev2 = fs2.stddev
          fs3 = mrdfits(file3, k, /silent)
          data3 = fs3.data
          exthdr3 = fs3.header
          stddev3 = fs3.stddev
          fs4 = mrdfits(file4, k, /silent)
          data4 = fs4.data
          exthdr4 = fs4.header
          stddev4 = fs4.stddev
 
          ; read inductosyn postion
          ind_pos1 = fxpar(exthdr1,'INDPOS')
          ind_pos2 = fxpar(exthdr2,'INDPOS')
          ind_pos3 = fxpar(exthdr3,'INDPOS')
          ind_pos4 = fxpar(exthdr4,'INDPOS')

          if (ind_pos1 eq ind_pos2) and (ind_pos1 eq ind_pos3) and $
          (ind_pos3 eq ind_pos4) then begin
              files_combed = file1+', '+file2+', '+file3+', '+file4
              comb_data = data1 + data2 + data3 + data4
              stddev_data = sqrt(stddev1^2 + stddev2^2 + stddev3^3 + stddev4^2)
              ;math is handled in chopsub, the Bs are already negative
         
              fxaddpar, exthdr1,  'NODFILES', files_combed, 'Files nod combined'
      
              outfile = str_replace(file1, '.chopsub.fits', '.nodcomb.fits')
              if ~file_test(outfile) then fxwrite, outfile, primehead, $
                                                   /noupdate, /append
              fits_struct = create_struct('header', exthdr1, 'data', comb_data,$
                                  'stddev', stddev_data)

              ; write the inductosyn extension for this particular unique
              ; channel+dither position+indpos
              mwrfits, fits_struct, outfile, /silent
           endif else if (ind_pos1 eq ind_pos2) and (ind_pos3 eq ind_pos4) then begin
              files_combed1 = file1+', '+file2
              files_combed2 = file3+', '+file4
              comb_data1 = data1 + data2
              comb_data2 = data3 + data4
              stddev_data1 = sqrt(stddev1^2 + stddev2^2)
              stddev_data2 = sqrt(stddev3^3 + stddev4^2)
              ;math is handled in chopsub, the Bs are already negative
         
              fxaddpar, exthdr1, 'NODFILES', files_combed1, 'Files nod combined'
              fxaddpar, exthdr3, 'NODFILES', files_combed2, 'Files nod combined'
              
              fits_struct1 = create_struct('header', exthdr1, 'data', $
                                           comb_data1, 'stddev', stddev_data1)
              fits_struct2 = create_struct('header', exthdr3, 'data', $
                                           comb_data2, 'stddev', stddev_data2)

              primehead3 = headfits(file3)

              if ~file_test(outfile1) then fxwrite, outfile1, primehead, $
                                                    /noupdate, /append
              if ~file_test(outfile2) then fxwrite, outfile2, primehead3, $
                                                    /noupdate, /append
              ; write the inductosyn extension for this particular unique
              ; channel+dither position+indpos
              mwrfits, fits_struct1, outfile1, /silent
              mwrfits, fits_struct2, outfile2, /silent
              fits_struct = fits_struct1
           endif
           
          k+=1
       endwhile

   endif
      
   i+=4
endwhile



         ;;files_combed=files_combed+' '+dl_filename
         ;math handled in subtract chops -jh
         ;if dl_nodbeam eq 'B' then begin
         ;   comb_data = comb_data - dl_data
         ;endif else if dl_nodbeam eq 'A' then begin
         ;;   comb_data = comb_data + dl_data
         ;endif else begin
         ;   print, 'nod is not A or B, not combining'
         ;endelse
         
         ;;stddev_data_sqd = stddev_data_sqd + dl_sigma^2 
         ;;combined = combined+1
         ;;file_for_this = dl_filename
    ;  endif

return, fits_struct
end
