;+
; NAME: fifi_ls_show_indwaves 
;
; PURPOSE: fifi_ls_show_indwaves FIFI-LS raw data 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=fifi_ls_show_indwaves(infile, flat)
;
; INPUTS: infile - root of the file to be fifi_ls_show_indwaves,
;                  including the
;                  path to it
;
; OPTIONAL INPUTS: flat - flatfield name, if not set the code will
;                         look for a flat_ind.csv file and use that to
;                         match the inductosyn positions in the extensions
;
; OUTPUTS: 1. Hypercube with flatfielded flux in one cube and lambdas
;             in the other cube
;          2. FITS file containing flat corrected flux in ext1, lambda
;          in ext2, xoffsets in ext3, yoffsets in ext4
;
; OPTIONAL OUTPUTS: 
;
; SIDE EFFECTS: 
;
; RESTRICTIONS: 
;
; PROCEDURE: 1. Go through each input file inductosyn extension
;            2. Either read the flatfield provided or find the
;            matching flat field extension using config file.
;            3. Divide by the flat in both spatial and spectral dimensions.
;            4. Create FITS file and write results to disk.
;
; EXAMPLE: fifi_ls_show_indwaves
;
; MODIFICATION HISTORY:
; 2015may30  J.H. first version

;-
pro fifi_ls_show_indwaves
;

     flat_info = read_csv('/Users/jholt/fifi_testdata/sky_dips-FIFI#10/outdir/flat_ind.csv')
     flatnames = flat_info.field1
     indposes = flat_info.field2
     ks = flat_info.field3
     br_nums = flat_info.field4

     if (file_test('flat_ind_waves.csv') gt 0) then begin
         openw, lun, 'flat_ind_waves.csv', /get_lun, /append
     endif else begin
         openw, lun, 'flat_ind_waves.csv', /get_lun
     endelse
    
     for i=0, n_elements(indposes)-1 do begin
        indpos = indposes[i]
        if br_nums[i] eq 0 then begin
           waves = fifi_ls_new_wave(indpos, [2015,01,01])
        endif else if br_nums[i] eq 1 then begin
           waves = fifi_ls_new_wave(indpos, [2015,01,01], blue='B1')
        endif else begin
           waves = fifi_ls_new_wave(indpos, [2015,01,01], blue='B2')
        endelse
        
         minw=min(waves)
         maxw=max(waves)
         printf, lun, strtrim(indpos, 2), ', ', strtrim(minw,2), ', ', strtrim(maxw,2)
     endfor
     free_lun, lun
  end

  
