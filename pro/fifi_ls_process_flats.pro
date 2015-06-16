;+
; NAME: fifi_ls_process_flats 
;
; PURPOSE: process a directory of sky-dips into flats 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=fifi_ls_process_flats(flatdir)
;
; INPUTS: flatdir - where raw sky dip data exists
;
; OPTIONAL INPUTS: outdir - where to write output flats
;
; OUTPUTS: 1. FITS files containing flatfield data [flux, weight, sigma] as a
;             a primary array, and the data header of the first data readout
;             an extension binary table.
;             File name: *.drizzle.L3.FITS
;             File content: 
;              prime header and prime array 
;                [# of x pixels,# of y pixels,# of lambda pixels,3]
;              binary table = data header of the first data readout
;
; OPTIONAL OUTPUTS: 
;
; KEYWORDS: CONFIGFILENAME - name of the file defining drizzle parameters,
;                            including the path to it.  If unused, a fallback
;                            file will be used.
;
; SIDE EFFECTS: 1. Adds astrometric keywords/values to FITS header.
;               2. SOFIA FITS keyword PROCSTAT gets updated to LEVEL_3.
;               3. Output FITS file now becomes prime array + binary table
;
; RESTRICTIONS: 1. x,y values are absolute mm when FITS keyword OBJ_NAME is
;                  TELSIM, but x,y offset in mm for other cases. 
;
; PROCEDURE: 
;
; EXAMPLE: result=fifi_ls_process_flats('00121_221031_FullScan-G_STRT_R-752500_lw')
;
; MODIFICATION HISTORY:
; 2015may28  J. Holt first version

;-
function fifi_ls_process_flats, flatdir

  ; move into flatdir
  cd, flatdir

  ; setup outdir
  outdir = flatdir+'/outdir/'
  if ~file_test(outdir, /DIRECTORY) then file_mkdir, outdir

  filenames=file_search('*Skydip*.fits')

  datalist = list()

  cd, current=rawdir
  
  i = 0
  while (i lt n_elements(filenames)) do begin

      cd, rawdir
      foo = fifi_ls_split_grating_and_chop(filenames[i], outdir)

      cd, outdir
      
      split_file0 = str_replace(filenames[i], '.fits', '.chop0.fits')
      split_file1 = str_replace(filenames[i], '.fits', '.chop1.fits')

      fits_info, split_file0, N_ext=n_ext0, /silent
      fits_info, split_file1, N_ext=n_ext1, /silent

      primehead = headfits(split_file0)

      ; read offset positions
      obslam = fxpar(primehead, 'OBSLAM')
      obsbet = fxpar(primehead, 'OBSBET')

      channel = strtrim(fxpar(primehead, 'CHANNEL'))
      if channel eq 'BLUE' then channel=1 else channel=0

      outfilename1 = str_replace(split_file0, '.chop0', '.all_chops')
      fxwrite, outfilename1, primehead, /noupdate, /append

      k = 1
      ; concatanate chops instead of subtracting
      while (k lt n_ext0+1) do begin

          f0 = mrdfits(split_file0, k, f0_exthdr, /unsigned, /silent)
          f1 = mrdfits(split_file1, k, f1_exthdr, /unsigned, /silent)

          hdr0=f0_exthdr
          hdr1=f1_exthdr

          data0=f0
          data1=f1

          ; read inductosyn postion
          ind_pos = fxpar(hdr0,'INDPOS')
          
          data_alles = [[[data0]], [[data1]]]
          
          hdr_alles = hdr0
 
          ;fits_struct = create_struct('header', hdr_alles, 'data', $
           ;                           data_alles)

          mwrfits, data_alles, outfilename1, hdr_alles, /silent
          ;mwrfits, fits_struct, outfilename1, /silent

          datalist.Add, list(filenames[i], ind_pos, channel, obslam, obsbet, k)
          k+=1
          
      endwhile

      ;now fit ramps
      infile1 = str_replace(filenames[i],'.fits', '.all_chops.fits')

      print,'fitting ramps for ', infile1
      foo = fifi_ls_fit_ramps(infile1, outdir, s2n=10, verbose=0)
      i+=1
      
  endwhile
  
  ; now need to go through and find matching channel & indpos
  ; and different OBSBETs
  jj = 0
  while (jj lt n_elements(filenames)) do begin
     dl_row = datalist[jj]
     dl_filename = dl_row[0]
     print,'dl_filename=', dl_filename,' jj=',jj
     dl_ind_pos = dl_row[1]
     dl_channel = dl_row[2]
     dl_obslam = dl_row[3]
     dl_obsbet = dl_row[4]
     dl_ext = dl_row[5]
     kk = 0
     while kk lt n_elements(datalist) do begin
         kdl_row = datalist[kk]
         kdl_filename = kdl_row[0]
         
         kdl_ind_pos = kdl_row[1]
         kdl_channel = kdl_row[2]
         kdl_obslam = kdl_row[3]
         kdl_obsbet = kdl_row[4]
         kdl_ext = kdl_row[5]
         
         if (dl_ind_pos eq kdl_ind_pos) and $
            (dl_channel eq kdl_channel) then begin
            if (dl_obslam gt kdl_obslam) or (dl_obsbet gt kdl_obsbet) then begin
               skydip1 = str_replace(dl_filename, '.fits', $
                                     '.all_chops.rampfits.fits')
               skydip2 = str_replace(kdl_filename, '.fits', $
                                     '.all_chops.rampfits.fits')
               print,'making flat for ',skydip1,' and ',skydip2, 'using extensions ', dl_ext,' and ',kdl_ext

               foo = fifi_ls_make_flat(skydip1, skydip2, dl_ext, kdl_ext)

               ; this makes flats n_ext times for each file
               ; 
            endif
         endif
         kk+=1
      endwhile
     
     jj+=1
  endwhile
 

  return, jj

end
