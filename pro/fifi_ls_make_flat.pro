;+
; NAME: make_flat 
;
; PURPOSE: drizzle FIFI-LS raw data 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=make_flat(skydip1, skydip2)
;
; INPUTS: skydip1 - fits file of first sky dip
;         skydip2 - fits file of second sky dip (must have matching
;                   indpos and channel, and skydip2 must have OBSLAM
;                   or OBSBET greater than skydip1.
;
; OPTIONAL INPUTS:
;
; OUTPUTS: 1. result - 1 if successful, -1 if failure
;          2. FITS file containing normalised flatfield data
;
; OPTIONAL OUTPUTS: 
;
; KEYWORDS: 
;
; SIDE EFFECTS: 1. Adds astrometric keywords/values to FITS header.
;               2. SOFIA FITS keyword PROCSTAT gets updated to LEVEL_3.
;               3. Output FITS file now becomes prime array + binary table
;
; RESTRICTIONS: 1.
;
; PROCEDURE: 1. Subtract data in skydip2 from skydip1
;            2. Fit a 5th order poly to the difference
;            3. divide skydip difference by polyfit
;            4. Create FITS file and write results to disk.
;
; EXAMPLE:
; result=fifi_ls_make_flat('00524_114014_0010_Skydip200um118_55_lw.fits', $
;         '00523_113922_0009_Skydip200um118_25_lw.fits')
;
; MODIFICATION HISTORY:
; 2015may27  J.H. first version

;-
function fifi_ls_make_flat, skydip1, skydip2, sky1_ext, sky2_ext

  cd, current=flatdir
  
  
  all_ext = 0
  if ~keyword_set(sky1_ext) then all_ext=1
  if ~keyword_set(sky2_ext) then all_ext=1

  if (all_ext eq 0) then begin
     n_ext = 1
  endif else begin
    fits_info, skydip1, n_ext=n_ext, /silent
  endelse
  
  primehead = headfits(skydip1)
  channel = strtrim(fxpar(primehead, 'CHANNEL'))
  b_order = fxpar(primehead, 'G_ORD_B')
  primehead2 = headfits(skydip1)
  channel2 = strtrim(fxpar(primehead2, 'CHANNEL'))
  b_order2 = fxpar(primehead2, 'G_ORD_B')

  if (channel ne channel2) then begin
     print, ' channels do not match'
     return, -1
  endif

  br_num = 0
  if (channel eq 'BLUE') then begin
         if b_order eq 1 then br_num = 1 else $
            if b_order eq 2 then br_num = 2 else begin
            print, 'Invalid Blue grating order.'
            return, -1
         endelse
  endif
  
  tmp = str_replace(skydip1, 'all_chops.', '')
  tmp = str_replace(tmp, '.rampfits', '')
  outfile = str_replace(tmp, '.fits', '.flat.fits')

  outfile_sub = 'sub_'+skydip1+'_'+skydip2+'.fits'
  
  if ~file_test(outfile) then fxwrite, outfile, $
                                       primehead, /noupdate, /append
   if ~file_test(outfile_sub) then fxwrite, outfile_sub, $
                                           primehead, /noupdate, /append
  if (file_test('flat_ind.csv') gt 0) then begin
     openw, lun, 'flat_ind.csv', /get_lun, /append
  endif else begin
     openw, lun, 'flat_ind.csv', /get_lun
  endelse
    
  k=1
  while (k lt n_ext + 1) do begin
     if (all_ext eq 1) then begin
        sky1_ext = k
        sky2_ext = k
     endif
     
     fits1_struct = mrdfits(skydip1, sky1_ext, /silent, /unsigned)
     fits2_struct = mrdfits(skydip2, sky2_ext, /silent, /unsigned)
     
     data1 = fits1_struct.data
     stddev1 = fits1_struct.stddev
     exthdr1 = fits1_struct.header
         
     data2 = fits2_struct.data
     stddev2 = fits2_struct.stddev
     exthdr2 = fits2_struct.header

     indpos = fxpar(exthdr1,'INDPOS')
     indpos2 = fxpar(exthdr2, 'INDPOS')
     if (indpos ne indpos2) then begin
        print, 'ERROR - cannot subtract different grating positions'
        return, -1
     endif
     
     flat = dblarr(5,5,16)
     i = 0
     while (i lt 5) do begin
         j = 0
         while (j lt 5) do begin
            flat_sub = data1[i, j, 1:16] - data2[i, j, 1:16]

            flat_sub_no_nans = flat_sub[where(finite(flat_sub))]

            coeffs = robust_poly_fit(indgen(n_elements(flat_sub_no_nans)), $
                                     flat_sub_no_nans, $
                                     5, yfit, sig, /double)
            newpoly = poly(indgen(16), coeffs)

            ;newpoly = yfit
            plot, flat_sub, psym=4
            oplot, newpoly
            wait, 0.1

            flat[i, j, *] = flat_sub / newpoly 
            
            if ~(finite(coeffs[0])) then begin
               print, 'bad fit'
               print,flat_sub
               wait, 1.0
            endif
            
            j+=1
         endwhile
         i+=1
     endwhile

     spat_flat = flat[*,*,8]
     ;spat_flat = median(flat, dimension=3)
     
        ;----------------------------------------------------------
     ; update fits header, write grating position to output file
     ;----------------------------------------------------------
     fits_struct = create_struct('header', exthdr1, 'spec_flat', flat, $
                                 'spat_flat', spat_flat, 'stddev', stddev1)

    ; if ~file_test(outfile) then fxwrite, outfile, $
     ;                                      primehead, /noupdate, /append
     mwrfits, fits_struct, outfile, /silent    

     print,'outfile=',outfile
     od = strtrim(outfile, 2)
     fd = strtrim(flatdir, 2)
     full_flat = fd+'/'+od
     ff = strcompress(full_flat) 
     br_num = string(br_num)
     
     printf, lun, strtrim(ff, 2)+', '+strtrim(indpos, 2)+ $
             ', '+strtrim(sky1_ext, 2)+', '+strtrim(br_num, 2)

     mwrfits, flat_sub, outfile_sub, /silent
     
     k+=1

  endwhile

  free_lun, lun
  
  return, 1
end
