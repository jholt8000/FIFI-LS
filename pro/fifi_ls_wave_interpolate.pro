;+
; NAME: wave_interpolate
;
; PURPOSE: Interpolate unevenly spaced FIFI-LS wavelength slices onto
; regular grid
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=fifi_ls_wave_interpolate(filename)
;
; INPUTS: filename - name of the file to be interpolated, including
;                    the path to it
;
; OPTIONAL INPUTS:
;
; OUTPUTS: 1. result - structure, -1 if failure
;          2. FITS file containing evenly spaced data [flux, weight, sigma] as a
;             a primary array, and the data header of the first data readout
;             an extension binary table.
;             File name: *.wave_interp.FITS
;             File content: 
;              prime header and prime array 
;                [# of x pixels,# of y pixels,# of lambda pixels,3]
;              binary table = data header of the first data readout
;
; OPTIONAL OUTPUTS: 
;
; KEYWORDS: 
;
; SIDE EFFECTS: 1. 
;               2. 
;               3. 
;
; RESTRICTIONS: 1. 
;
; PROCEDURE: 1. Read input file.
;            2. Define output grid based on values in either data.lambda
;            3. Interpolate data at each spaxel
;            4. Create FITS file and write results to disk.
;
; EXAMPLE: result=fifi_ls_wave_interpolate('00121_221031_FullScan-G_STRT_R-752500_lw')
;
; MODIFICATION HISTORY:
; 2015may01  jholt:  first version

;-
function fifi_ls_wave_interpolate, infile

  primehead = headfits(infile)
  in_struct = mrdfits(infile, 1)
  flux = in_struct.data
  lambda = in_struct.lambda
  xs = in_struct.xs
  ys = in_struct.ys
  stddev = in_struct.stddev
  print,'size_std=',size(stddev)
  hypercube = dblarr(5, 5, n_elements(flux)/25, 3)
  print, n_elements(flux)
  
  print, size(flux)
  hypercube[*, *, *, 0] = lambda ;reform(lambda, 5, 5, n_elements(lambda[*,*]) / 25)
  hypercube[*, *, *, 1] = flux ;reform(flux, 5, 5, n_elements(flux) / 25)

  hypercube[*, *, *, 2] = stddev
  
  ; make a master grid
  wmin = min(hypercube[*,*,*,0])
  wmax = max(hypercube[*,*,*,0])
  sample_size = 2.5
  pixscale = 0.02
 
  ;dw = pixscale / sample_size
  fwhm = 0.0478
  sampl = 2.5
  
  novsamp = fwhm*sampl/pixscale
  dw = fwhm / novsamp
  nw = (wmax - wmin) / dw + 1
  
  master_grid = findgen(nw) * dw + wmin
  ;master_grid_bv = findgen(109.188) * 0.004 + 51.5879
  
  new_cube = dblarr(5, 5, 2, n_elements(master_grid))

  window,0
  !P.MULTI = [0,5,5]
  !X.MARGIN = !X.MARGIN/2.
  !Y.MARGIN = !Y.MARGIN/2.
  loadct,1
  i=0
  while (i lt 5) do begin
     j=0
     while (j lt 5) do begin

        bill = 1

        print, ' starting on i, j =', i, j
        if bill eq 1 then begin
           farr = compflux_arrs(hypercube[i, j, *, 0], $
                                hypercube[i, j, *, 1], $
                                ERR=hypercube[i, j, *, 2], $
                                FITORDR=5, WMIN=wmin, WMAX=wmax, $
                             FACTR=2.0, SAMPL=sample_size, $
                             master_grid=master_grid, $
                               /NEGCUT)
        endif else begin
           farr = fifi_ls_compflux(hypercube[i, j, *, 0], $
                                hypercube[i, j, *, 1], $
                                ERR=hypercube[i, j, *, 2], $
                                master_grid=master_grid, FITORDR=5, $
                                FACTR=2.0, SAMPL=sample_size, $
                               /NEGCUT)
        endelse
        new_cube[i,j,0,*] = master_grid
        new_cube[i,j,1,*] = farr

        plot, hypercube[i,j,*,0], hypercube[i,j,*,1], psym=3, yrange=[-1,median(hypercube[*,*,*,1])+30]
        oplot, master_grid, farr,psym=10

        print, 'ended i, j=', i, j
        
        j+=1
    endwhile
  i+=1
     
endwhile
  !P.MULTI=0

  
  ;plots, [!X.CRANGE[0],!X.CRANGE[1]], [0.0,0.0], linestyle=2

  x0=fxpar(primehead,'OBSLAM')
  y0=fxpar(primehead,'OBSBET')  
  rot = fxpar(primehead,'DET_ANGL')
;want ra in deg = hour * 15 + min / 4 + sec /240
  DA      = fxpar(primehead,'DET_ANGL') * !DTOR ; ccw array rotation
  
  outfile=str_replace(infile, '.fits', '.lgrid.fits')

  ;fits_struct = create_struct('header', exthdr, 'data', data, $
  ;                                'lambda', lambda, 'stddev', stddev)
  if (file_test(outfile) gt 0) then file_delete, outfile

  if ~file_test(outfile) then fxwrite, outfile, primehead, /noupdate, /append
  mwrfits, new_cube, outfile, /silent

  outfile=str_replace(infile, '.fits', '.lgrid_fluxonly.fits')

  ;fits_struct = create_struct('header', exthdr, 'data', data, $
  ;                                'lambda', lambda, 'stddev', stddev)

  if (file_test(outfile) gt 0) then file_delete, outfile
  
  mwrfits, reform(new_cube[*,*,1,*],5,5,n_elements(master_grid)), outfile, primhead,  /silent

  return, new_cube
  end
