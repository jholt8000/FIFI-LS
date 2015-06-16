;+
; NAME: call_drizzle 
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
; RESTRICTIONS: 1. x,y values are absolute mm when FITS keyword OBJ_NAME is
;                  TELSIM, but x,y offset in mm for other cases. 
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
function fifi_ls_call_drizzle, infile

  ; make slices for each lambda
  

  primehead = headfits(infile)

  cube = mrdfits(infile, 1, /silent, /unsigned)

  l=0
  while (l lt n_elements(cube[1,1,0,*])) do begin
     plane = reform(cube[*,*,1,l],5,5)
     lambda = cube[0,0,0,l] ;lambda same along plane
     print,plane
     print,lambda
     outfile = str_replace(infile, '.fits', strtrim(string(lambda),2)+'.fits')
     print,'outfile=',outfile
     ;fxwrite, outfile, primehead, /noupdate, /append
     mwrfits, plane, outfile, primehead, /silent
     l+=1
     
  endwhile

  return, cube
 end
