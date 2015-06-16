;+
; NAME: drizzle 
;
; PURPOSE: drizzle FIFI-LS raw data 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=drizzle(FILEROOT,CONFIGFILENAME=CONFIGFILENAME,verbose=verbose)
;
; INPUTS: FILEROOT - root of the file to be drizzled, including the path to it
;
; OPTIONAL INPUTS:
;
; OUTPUTS: 1. result - 1 if successful, -1 if failure
;          2. FITS file containing drizzled data [flux, weight, sigma] as a
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
; RESTRICTIONS: 1. x, y values are absolute mm when FITS keyword OBJ_NAME is
;                  TELSIM, but x, y offset in mm for other cases. 
;
; PROCEDURE: 1. Read input file.
;            2. Define output grid (4-D cube) based on values in either
;               CONFIGFILENAME or the fallback file.
;            3. Drizzle data.
;            4. Create FITS file and write results to disk.
;
; EXAMPLE: result=fifi_ls_wave_interpolate('00121_221031_FullScan-G_STRT_R-752500_lw')
;
; MODIFICATION HISTORY:
; 2015may01  J.H. first version

;-
function fifi_ls_wave_interpolate_old, infile

  device, decomposed = 0
  loadct, 13
  
  fits_info, infile, n_ext = n_ext, /silent

  ncol = 255 / (n_ext)
  colors = ncol*indgen(n_ext+1)+ncol

  psym_arr = indgen(7) + 1
  psym_arr[2] = 6
  psyms = [psym_arr, psym_arr, psym_arr, psym_arr, psym_arr, psym_arr]
  
 ; f2 =  [0.948657,    1.07551 ,   1.05191,    1.00710 ,  0.900825 ,  0.965355,   0.964521 ,   1.10460 ,   1.01417,    1.09404 ,  0.953356,   0.982923,   0.905014  ,  1.03718,    1.09085,   0.918496]

  ;f2 = 1/f2
  ;  f3 = [f2, f2,f2,f2,f2,f2,f2,f2,f2,f2,f2,f2,f2,f2]
  fits_info, infile, n_ext = n_ext, /silent

  primehead = headfits(infile)
  channel = strtrim(fxpar(primehead, 'CHANNEL'), 2)

  k = 1

  while (k lt n_ext + 1) do begin
     
     fits_struct = mrdfits(infile, k, /silent, /unsigned)
     data = fits_struct.data
     data = reform(data, 5, 5, n_elements(data)/25)
     
     lambda = fits_struct.lambda
     lambda = reform(lambda, 5, 5, n_elements(lambda)/25)
     stddev = fits_struct.stddev
     exthdr = fits_struct.header
     xs = reform(fits_struct.xs, 5, 5, 16)
     ys = reform(fits_struct.ys, 5, 5, 16)
     
     i = 0
     while (i lt 5) do begin
        j = 0
        while (j lt 5) do begin
           if (k eq 1) and (i eq 0) and (j eq 0) then  begin
              lambda_pix = lambda[i, j, *]            
              flux_pix = data[i, j, *]
              xs_pix = xs[i, j, *]
              ys_pix = ys[i, j, *]
               
           endif else begin
              lambda_pix = [lambda_pix, lambda[i, j, *]]
              flux_pix = [flux_pix, data[i, j, *]]
              xs_pix = [xs_pix, xs[i, j, *]]
              ys_pix = [ys_pix, ys[i, j, *]]             
              if (i eq 2) and (j eq 2) then begin
                 if (k eq 1) then begin
                  plot, lambda[2, 2, *], data[2, 2, *], psym = 1, $
                       xrange = [min(lambda[2, 2, *]), max(lambda[2, 2, *])], $
                        yrange = [min(data[2, 2, *])-1, $
                                  max(data[2, 2, *])-0.1], color = colors[1]
                  oplot, lambda[2, 2, *], data[2, 2, *], color = colors[1]
                 endif else begin
                    oplot, lambda[2, 2, *], data[2, 2, *], psym = psyms[k], $
                           color = colors[k]
                   oplot, lambda[2, 2, *], data[2, 2, *], color = colors[k]
                        
                 endelse
               
              endif else if i eq 2 and j eq 2 then begin
                 oplot, lambda[2, 2, *], psym = psyms[k], color = colors[k]
              endif                             
           endelse
           j+ = 1
           
        endwhile
        i+ = 1
        
     endwhile
     k+ = 1
     
  endwhile

  print, size(lambda_pix)
  rlp = reform(lambda_pix, 5, 5, n_elements(lambda_pix)/25)
  rfp = reform(flux_pix, 5, 5, n_elements(flux_pix)/25)
  ;rlp = reform(lambda_pix, 5, 5, 16*n_ext)
  ;rfp = reform(flux_pix, 5, 5, 16*n_ext)

  ;t1_rfp = rfp[2, 2, *]
  ;t1_rlp = rlp[2, 2, *]
  ;window, 1
  ;plot, rlp[2, 2, *], rfp[2, 2, *], psym = 2  , title = 'no flat'
  ;print, 's1_22 = ', size(rfp[2, 2, *, *])

  ;print, 's1 = ', size(rfp)
  ;rlp = reform(lambda_pix, 5, 5, 16, n_ext)
  ;rfp = reform(flux_pix, 5, 5, 16, n_ext)
  ;print, 's2_22 = ', size(rfp[2, 2, *, *])

 ; print, 's2 = ', size(rfp)
;  test1 = reform(rfp[2, 2, *, *]) / f3

  ;print, 's3 = ', size(rfp)
  ;print, 's3_22 = ', size(rfp[2, 2, *, *])
  ;rlp = reform(rlp, 5, 5, 16*n_ext)
  ;rfp = reform(rfp, 5, 5, 16*n_ext)

  
  ;for k = 0, n_ext-1 do begin
    ; rfp[2, 2, *, *] = reform(rfp[2, 2, *, *]) * f3
     ;print, k = 
  ;endfor
  
  ;window, 2  
  ;plot, rlp[2, 2, *], test1, psym = 2  , title = 'div by flat '

;  hypercube = dblarr(5, 5, 2, 16*n_ext)
  hypercube = dblarr(5, 5, 2, n_elements(flux_pix)/25)
  ; have to divide flat for each ind pos before we sort hypercube
  
  ; don't really need to make hypercube, can use rfp and rlp
  i = 0
  while (i lt 5) do begin
     j = 0
     while (j lt 5) do begin
        lp = rlp[i, j, *]
        fp = rfp[i, j, *]

        print, size(hypercube[i, j, 1, *])
        
        hypercube[i, j, 0, *] = lp[sort(lp)]
        hypercube[i, j, 1, *] = fp[sort(lp)]
       
        j+ = 1
     endwhile
     i+ = 1
  endwhile

  newarr = dblarr(5, 5, n_elements(hypercube[0, 0, 1, *])/2.5, n_elements(hypercube[0, 0, 0, *])/2.5)

  ; make a master grid
  wmin = min(hypercube[*, *, 0, *])
  wmax = max(hypercube[*, *, 0, *])
  fwhm = 0.0478
  sample_size = 5
  
 ; dw = fwhm / sample_size
  dw = 0.02 / sample_size
  nw = (wmax - wmin) / dw + 1
  
  master_grid_jh = findgen(nw) * dw + wmin
  master_grid_bv = findgen(109.188) * 0.004 + 51.5879
  
  master_grid =  master_grid_jh
    
  print, 'master_grid_jh = ', master_grid_jh
  print, 'master_grid_bv = ', master_grid_bv
  
  new_cube = dblarr(5, 5, 2, n_elements(master_grid))

  window, 4
  !P.MULTI = [0, 5, 5]
  ;!X.MARGIN = !X.MARGIN/2.
  ;!Y.MARGIN = !Y.MARGIN/2.
  loadct, 0
  i = 0
  while (i lt 5) do begin
     j = 0
     while (j lt 5) do begin
                                ;newarr[i, j, *]
        
        farr = fifi_ls_compflux(wobs = hypercube[i, j, 0, *], $
                                fobs = hypercube[i, j, 1, *], $
                        master_grid = master_grid, FITORDR = 5.0, FACTR = 2.5, $
                        SAMPL = sample_size, SFACT = 2.0);, /NOFILT)

        new_cube[i, j, 0, *] = master_grid
        new_cube[i, j, 1, *] = farr
        
        plot, hypercube[i, j, 0, *], hypercube[i, j, 1, *], psym = 3, yrange = [-1, max(hypercube[*, *, 1, *])]
       ; oplot, master_grid, farr
;oplot, hypercube[i, j, 3, *], fout, psym = 5, color = 3
       ; if i eq 1 and j eq 2 then begin
       ;    farr_12 = farr
       ;    master_grid_12 = master_grid
       ;    openw, lun, 'new12_noisy.dat', /get_lun
       ;   for ll = 0, n_elements(hypercube[i, j, 1, *])-1 do begin
       ;       printf, lun, hypercube[1, 2, 0, ll], hypercube[1, 2, 1, ll]
       ;    endfor
       ;    close, lun
       ;    free_lun, lun
       ; endif
        
        j+ = 1
    endwhile
  i+ = 1
     
endwhile
  !P.MULTI = 0
  ;window, 3
 
  ;plot, hypercube[2, 2, 0, *], hypercube[2, 2, 1, *], psym = 4, title = 'with flat';, yrange = [-1, 2]
  ;oplot, new_cube[2, 2, 0, *], new_cube[2, 2, 1, *], psym = 4
  ;oplot, new_cube[2, 2, 0, *], new_cube[2, 2, 1, *], psym = 10

  ;for oo = 0, n_elements(hypercube[2, 2, 1, *])/16.-1 do begin
  ;    oplot, hypercube[2, 2, 0, *], hypercube[2, 2, 1, *]/ f2, psym = 4
   ;endfor
  
  plots, [!X.CRANGE[0], !X.CRANGE[1]], [0.0, 0.0], linestyle = 2

  ;!P.MULTI = 0
  ;window, 1
  ;plot, hypercube[2, 2, 1, *], hypercube[2, 2, 0, *], psym = 3
  ;oplot, master_grid_22, farr_22, psym = 5
  ;oplot, master_grid_22, farr_22, psym = 10
 
  x0 = fxpar(primehead, 'OBSLAM')
  y0 = fxpar(primehead, 'OBSBET')  
  rot = fxpar(primehead, 'DET_ANGL')
;want ra in deg = hour * 15 + min / 4 + sec /240
  DA      = fxpar(primehead, 'DET_ANGL') * !DTOR ; ccw array rotation
  
  outfile = str_replace(infile, '.fits', '.lgrid.fits')

  ;fits_struct = create_struct('header', exthdr, 'data', data, $
  ;                                'lambda', lambda, 'stddev', stddev)
  if (file_test(outfile) gt 0) then file_delete, outfile

  if ~file_test(outfile) then fxwrite, outfile, primehead, /noupdate, /append
  mwrfits, new_cube, outfile, /silent

  outfile = str_replace(infile, '.fits', '.lgrid_fluxonly.fits')

  ;fits_struct  = create_struct('header', exthdr, 'data', data, $
  ;                                'lambda', lambda, 'stddev', stddev)

  if (file_test(outfile) gt 0) then file_delete, outfile
  
  mwrfits, reform(new_cube[*, *, 1, *], 5, 5, n_elements(master_grid)), outfile, primhead, /silent

  return, new_cube
  end
