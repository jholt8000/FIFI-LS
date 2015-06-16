;+
; NAME: fifi_ls_comb_grating_scans 
;
; PURPOSE: Combine separate grating positions in MEF FIFI-LS raw data
; into a single extension. Will produce an unevenly spaced wavelength
; grid if multiple grating (inductosyn) extensions in input file
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=fifi_ls_comb_grating_scans(filename,
; sort_this, plot_this, write_hypercube)
;
; INPUTS: filename: name of file to be combined
;
; OPTIONAL INPUTS:
;                  sort_this = 1: sort the wavlength and fluxes in
;                                 order of wavelength, recommended
;                  plot_this = 1: create plots of the center spaxel
;                  write_hypercube = 1: write just the hypercube of
;                  wavelength/flux to be read with DS9.
;
; OUTPUTS: 1. result - array with wavelength and flux, if successful, -1 if failure
;          2. FITS file containing combined data  as a structure with
;          data.header as the extension header, data.data as the flux cube,
;          data.lambda as the wavelength cube, data.stddev as the
;          error cube, data.xs and data.ys as the detector spatial
;          offsets
;
;         
;
; OPTIONAL OUTPUTS: 
;
; KEYWORDS: 
;
; SIDE EFFECTS: 1. 
;               2. 
;               3. 
;
; RESTRICTIONS: 1. Needs a lambda and spatial calibrated input file
;
; PROCEDURE: 1. Read input file extensions and build master lambda and
;               flux arrays
;            2. Sort flux and wavelength, if sort_this=1
;            3. Create FITS file and write results to disk.
;
; EXAMPLE: result=fifi_ls_comb_grating_scans('file.lambda_xy.fits')
;
; MODIFICATION HISTORY:
; 2015may01  jholt: first version

;-
function fifi_ls_comb_grating_scans, infile, sort_this, plot_this, $
                                     write_hypercube

  if ~keyword_set(plot_this) then plot_this = 1
  if ~keyword_set(sort_this) then sort_this = 1
  if ~keyword_set(write_hypercube) then write_hypercube = 1
  ; TODO -- write n_elements check for when eq 0
  
  fits_info, infile, n_ext=n_ext, /silent

  outfile = str_replace(infile, '.fits', '.scancomb.fits')
  
  if plot_this then begin
      device, decomposed=0
      ;loadct, 13
      
      ncol=255/(n_ext)
      colors = ncol*indgen(n_ext+1)+ncol

      ;psym_arr = indgen(7) + 1
      ;psym_arr[2] = 6
      ;psyms=[psym_arr, psym_arr, psym_arr, psym_arr, psym_arr, psym_arr]
      psyms=['*','+','s','*','p','d','o','h','x','+','*','s','d','+','*','+','s','*','p',$
             'd','o','h','x','+','*','s','d','+']
      colors=['red','orange','black','green','blue','purple','cyan','red','orange','black',$
              'green','blue','purple','cyan','red','orange','black','green','blue','purple',$
              'cyan','red','orange','black','green','blue','purple','cyan']
  endif
  
  primehead = headfits(infile)
  
  k=1

  while (k lt n_ext+1) do begin
     
     fits_struct = mrdfits(infile, k, /silent, /unsigned)
     data = fits_struct.data
     ;data = reform(data, 5, 5, 16)     
     lambda = fits_struct.lambda
     ;lambda = reform(lambda, 5, 5, 16)
 
     stddev = fits_struct.stddev
     exthdr = fits_struct.header
     xs = reform(fits_struct.xs, 5, 5, 16)
     ys = reform(fits_struct.ys, 5, 5, 16)

     ; first go through and make master arrays for all flux, lamba, x, y
     ;in file (this combines all extensions/grat positions)
     i=0
     while (i lt 5) do begin
        j=0
        while (j lt 5) do begin
           if (k eq 1) and (i eq 0) and (j eq 0) then  begin
              print,'size l1=', size(lambda[i, j, *])
        
              lambda_pix = lambda[i, j, *]
              flux_pix = data[i, j, *]
              xs_pix = xs[i, j, *]
              ys_pix = ys[i, j, *]
              stddev_pix = stddev[i,j,*]
              
           endif else begin
              
              lambda_pix = [lambda_pix, lambda[i, j, *]]
              flux_pix = [flux_pix, data[i, j, *]]
              xs_pix = [xs_pix, xs[i, j, *]]
              ys_pix = [ys_pix, ys[i, j, *]]
              stddev_pix = [stddev_pix, stddev[i, j, *]]

              if plot_this then begin
                 if (i eq 2) and (j eq 2) then begin
                    if k eq 1 then begin
                       loadct,0
                   
                       p=plot(lambda[2, 3, *], data[2, 3, *], symbol='*', $
                              xrange=[min(lambda[2, 3, *])-0.03, $
                                      max(lambda[2, 3,*])+0.1],$
                             ; yrange = [-30, 160], $
                             yrange=[min(data[2, 3,*])-3, max(data[2, 3,*])+2],$
                             title='mult flatfield', xtitle='wavelength (micron)',$
                             color='black',linestyle=6)
                    endif else begin
                       loadct,13
                       
                       p=plot(lambda[2, 3,*], data[2, 3,*], symbol=psyms[k], $
                              color=colors[k],linestyle=6,/overplot)
                     
                    endelse
                 endif
                 
                    
             endif ; end plotting
              
           endelse

           j+=1
           
        endwhile
        i+=1
        
     endwhile

     k+=1
  endwhile

  p.Save,"center_spaxel.png",border=10,resolution=300
  
  rlp = reform(lambda_pix, 5, 5, 16*n_ext)
  rfp = reform(flux_pix, 5, 5, 16*n_ext)
  rstd = reform(stddev_pix, 5, 5, 16*n_ext)
  rxs = reform(xs_pix, 5, 5, 16*n_ext)
  rys = reform(ys_pix, 5, 5, 16*n_ext)
  ;forplot = median(rfp, dimension=3)
  ;tv,forplot
 
  hypercube = dblarr(5, 5, 16*n_ext, 5)
  i=0
  while (i lt 5) do begin
     j=0
     while (j lt 5) do begin
        lp = rlp[i, j, *]
        fp = rfp[i, j, *]
        st = rstd[i, j, *]
        x = rxs[i, j, *]
        y = rys[i, j, *]
        
        if sort_this then begin
            hypercube[i, j, *, 0] = lp[sort(lp)]
            hypercube[i, j, *, 1] = fp[sort(lp)]
            hypercube[i, j, *, 2] = st[sort(lp)]
            hypercube[i, j, *, 3] = x[sort(lp)]
            hypercube[i, j, *, 4] = y[sort(lp)]
        endif else begin
            hypercube[i, j, *, 0] = lp
            hypercube[i, j, *, 1] = fp
            hypercube[i, j, *, 2] = st
            hypercube[i, j, *, 3] = x
            hypercube[i, j, *, 4] = y
         endelse
        
        j+=1
     endwhile
     i+=1
  endwhile  

  ys = hypercube[*, *, *, 4]
  xs = hypercube[*, *, *, 3]
  stddev = hypercube[*, *, *, 2]
  fluxes = hypercube[*, *, *, 1]
  lambdas = hypercube[*, *, *, 0]
  if ~file_test(outfile) then fxwrite, outfile, $
                                       primehead, /noupdate, /append

  ;----------------------------------------------------------
  fits_struct = create_struct('header', exthdr, 'data', fluxes, $
                              'lambda', lambdas, 'stddev', stddev, $
                              'xs', xs, 'ys', ys)

  mwrfits, fits_struct, outfile, /silent    

  if write_hypercube then begin
      of1=str_replace(infile, '.fits', '.lambda_flux_stddev_hypercube.fits')
      fxwrite, of1, primehead, hypercube, /noupdate
  endif
  
  return, fits_struct
end
