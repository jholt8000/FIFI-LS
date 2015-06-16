;+
; NAME: fifi_ls_apply_flat 
;
; PURPOSE: fifi_ls_apply_flat FIFI-LS raw data 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=fifi_ls_apply_flat(infile, flat)
;
; INPUTS: infile - root of the file to be fifi_ls_apply_flat,
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
; EXAMPLE: result=fifi_ls_apply_flat('00121_221031_FullScan-G_STRT_R-752500_lw')
;
; MODIFICATION HISTORY:
; 2015may30  J.H. first version

;-
function fifi_ls_apply_flat, infile, flat
;

  use_csv = 0
  if ~keyword_set(flat) then use_csv = 1
  
  if use_csv then begin
     flat_info = read_csv('/Users/jholt/fifi_testdata/sky_dips-FIFI#10/outdir/flat_ind.csv')
     flatnames = flat_info.field1
     indposes = flat_info.field2
     ks = flat_info.field3
     br_nums = flat_info.field4
  endif

  outfile = str_replace(infile, '.nodcomb', '')
  outfile = str_replace(outfile, '.lambda','')
  outfile = str_replace(outfile, '_xy','')
  
  outfile = str_replace(outfile, '.fits', '.ff.fits') 
  fits_info, infile, n_ext=n_ext, /silent
  primehead = headfits(infile)
     
  ; if the output file already exists, delete it
  if file_test(outfile) gt 0 then file_delete, outfile

  ; write the file with the primary header and make it appendable 
  fxwrite, outfile, primehead, /noupdate, /append

  k = 1
  while (k lt n_ext+1) do begin
     
     fits_struct = mrdfits(infile, k, /silent, /unsigned)
     data = reform(fits_struct.data, 5, 5, 16)
     lambda = reform(fits_struct.lambda, 5, 5, 16)
     stddev = reform(fits_struct.stddev, 5, 5, 16)
     exthdr = fits_struct.header
     indpos_sci = fxpar(exthdr, 'INDPOS')
     
     xs = reform(fits_struct.xs, 5, 5, 16)
     ys = reform(fits_struct.ys, 5, 5, 16)

     j=0
     print, n_elements(flatnames)
     while (j lt n_elements(flatnames)-1) do begin
        flatname = flatnames[j]
        indpos_flat = indposes[j]
        print,'flatname=',flatname,' indpos_sci=',indpos_sci, k
        print,'indpos_flat=',indpos_flat, k
        print,'diff=',abs(indpos_sci - indpos_flat), k
        ext = ks[j]
        if (abs(indpos_sci - indpos_flat) lt 150000) then begin
           print,'using inpos flat=', indpos_flat, ' for indpos_sci=',indpos_sci, k
           flat_struct = mrdfits(flatname, ext, /silent)
           flat_spec = flat_struct.spec_flat
           flat_spat = flat_struct.spat_flat
        endif
        j += 1
     endwhile

     flat_corr = dblarr(5, 5, 16)
     
     ll = 0
     while (ll lt 5) do begin
        mm = 0
        while (mm lt 5) do begin

           flat_corr[ll, mm, *] = data[ll,mm,*] / flat_spec[ll,mm,*]           
            
           mm += 1
        endwhile
        
        ll += 1
     endwhile

     nn = 0
     while (nn lt 16) do begin
        ;plot, data[2,2,*],psym=4
        ;wait,0.5
        ;oplot, flat_corr[2,2,*], psym=6
        ;print,'old flat corr[*,*,nn]=',flat_corr[*, *, nn]
        ;wait,0.5       
        flat_corr[*, *, nn] = flat_corr[*, *, nn] / flat_spat

        ;print,'new flat corr=',flat_corr[*, *, nn]
        ;print,'flat_spat=',flat_spat
        ;oplot, flat_corr[2, 2, *], psym=5
        
        ;wait,0.1
        
                
        nn += 1
     endwhile
    ; 
        
     fits_struct = create_struct('header', exthdr, 'lambda', lambda, 'data', $
                              flat_corr, 'stddev', stddev, 'xs', xs, 'ys', ys)

     mwrfits, fits_struct, outfile, /silent

     k += 1
  endwhile
  
  return, fits_struct
end
