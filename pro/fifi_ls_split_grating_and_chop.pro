;+
; NAME: split_grating_and_chop 
;
; PURPOSE: Split a FITS raw data file into two FITS files based on
; chop cycle with each of the grating postions in multiple FITS extensions
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: RESULT=split_grating_and_chop(FILENAME, outdir, verbose)
;
; INPUTS: FILENAME - Name of the FITS raw data file, including the directory
;                    path if not in the current working directory.
;         outdir - name of the output path
;
; OUTPUTS: 1. RESULT - data structure
;          2. FITS files containing raw data of a certain chopper phase.
;             File name: : *.chop#.fits, where # is 0 or 1
;             File content:
;              prime header
;              for each inductosyn (grating) postion an extension with
;              extension header, with INDPOS as a keyword
;              binary table - (header + [26,18] data array) x (# of frames)
;
; SIDE EFFECTS: 1. Removes any previously splitted files (this is to avoid 
;                  appending data to existing files). 
;               2. Removes any partial ramps at beginning or end of
;               file. Does it still do this? 
;               3. Removes any "unpaired" chop plateaus.  2POINT chop tested; 
;                  4POINT chop NOT tested. Does it still do this?
;
; RESTRICTIONS: FIFILS keyword CH_SCHEME must be 2POINT 
;               4POINT not currently working
;
; PROCEDURE: 1. Read input file, a raw data FITS file.
;            2. Remove any partial ramps, or unpaired chop plateaus
;            3. Split data according to chopper phase, which is determined
;               from the FIFILS keyword CH_SCHEME and split according
;               to multiple inductosyn positions. 
;            4. Write output files to disk.
;
; EXAMPLE: result = split_grating_and_chop('00001_lw.fits', 'outdir')
;
; MODIFICATION HISTORY:
; 25feb08 k.n.  started coding...
;
; 12mar08 k.n.  code runs for cases where (one chop position) = (one ramp)
;
; 12mar08 k.n.  started modifying code to accomodate cases where
;               (one chop position) ne (one ramp)
;
; 24sep08 k.n.  bug fixes; added CVS version in prime header
;
; 25jun09 k.n.  square brackets ([]) now allowed in file name
;
; 27oct14 k.n.  1. added code to remove any partial ramps at beginning or
;                  end of file
;               2. added code to remove any "unpaired" chop plateaus.  Note
;                  that while 2POINT chop is tested, 4POINT chop is
;                  not.
;
; 17feb15 J.H.  Update mrdfits call to use filename instead of file
;               unit to work with new astrolib
;
; 2015Mar22 JH  Overhaul and rewrite


function fifi_ls_split_grating_and_chop, filename, outdir, verbose = verbose

if ~keyword_set(verbose) then verbose = 0
;--------------
;-- help file
;--------------
if n_params() eq 0 then begin
  doc_library, 'fifi_ls_split_grating_and_chop'
  return, -1
endif

;------------------
;-- sanity checks 
;------------------
;---- does the raw FITS file exist?  if not, return -1 and exit
if file_test(filename) ne 1 then begin
  message, 'split_grating_and_chop.pro: input FITS file does not exist.'
  return, -1
endif

;---- has the raw FITS file already been split?  if so, remove them to avoid
;---- appending data to existing file
basename = strmid(filename, 0, strlen(filename)-4)
x = strsplit(filename, '/', /extract)
base_filename = strmid(x[n_elements(x)-1], 0, strlen(x[n_elements(x)-1])-4)

fitsfile0 = outdir+'/'+base_filename+'chop0.fits'
fitsfile1 = outdir+'/'+base_filename+'chop1.fits'

foo = file_search(outdir+'/'+base_filename+'chop*fits')
if foo[0] ne '' then file_delete, foo[0]

;---------------------------
;-- define some parameters - these are locations in the data table
; header to find chop info
;---------------------------
P_FLAG = 3L
P_SMPL_CNT = 4L
P_RMP_CNT = 5L

;----------------------------
;-- read FITS raw data file
;----------------------------
primehead = headfits(filename)

data = mrdfits(filename, 1, /unsigned,/silent)

C_SCHEME = strtrim(fxpar(primehead, 'C_SCHEME'), 2)	; chopper scheme
CHOP_LN = fxpar(primehead, 'C_CHOPLN')	; choplength
CHOP_LN = fix(strtrim(CHOP_LN, 2))
channel = fxpar(primehead, 'CHANNEL')

CHPFREQ = fxpar(primehead, 'CHPFREQ') ; chopping frequency in Hz

if (C_SCHEME eq '2POINT') then begin
    C = CHOP_LN * 1.
endif else begin
    print, 'Invalid chopper scheme. Exiting split_cycle.pro.'
    return, -1
endelse   

if (size(channel, /tn) eq 'INT') then begin
      redBlue = ((data[0].header)[P_FLAG]  and 2) / 2 ; yields 0 for red and
                                ;        1 for blue
      if (redBlue eq 0) then channel = 'RED'
      if (redBlue eq 1) then channel = 'BLUE'
endif
   
;case channel of
if (strtrim(channel, 2) eq 'RED') then begin
     G_STRT = fxpar(primehead, 'G_STRT_R'); grating/inductosyn starting position
     G_PSUP = fxpar(primehead, 'G_PSUP_R'); # of grating pos. going up
     G_SZUP = fxpar(primehead, 'G_SZUP_R'); step size on the way up
     G_PSDN = fxpar(primehead, 'G_PSDN_R'); # of grating pos. coming down
     G_SZDN = fxpar(primehead, 'G_SZDN_R'); step size on the way down
     G_CYC = fxpar(primehead, 'G_CYC_R')	; # of grating sweeps (up-down)
     C_CYC = fxpar(primehead, 'C_CYC_R')	; # of chop cycles / grat. pos.
     ramplength = fxpar(primehead, 'RAMPLN_R') ; ramplength
     fxaddpar, primehead, 'CHANNEL', 'RED'
endif else if (strtrim(channel, 2) eq 'BLUE') then begin
     G_STRT = fxpar(primehead, 'G_STRT_B'); grating starting position
     G_PSUP = fxpar(primehead, 'G_PSUP_B'); # of grating positions on the up 
     G_SZUP = fxpar(primehead, 'G_SZUP_B'); step size on the way up
     G_PSDN = fxpar(primehead, 'G_PSDN_B'); # of grating positions coming down
     G_SZDN = fxpar(primehead, 'G_SZDN_B'); step size on the way down
     G_CYC = fxpar(primehead, 'G_CYC_B')	; # of grating sweeps (up-down)
     C_CYC = fxpar(primehead, 'C_CYC_B') ; # of chop cycles / grating position
     ramplength = fxpar(primehead, 'RAMPLN_B') ; ramplength
     fxaddpar, primehead, 'CHANNEL', 'BLUE'
endif else begin
        print, 'Invalid spectral channel.  Exiting split_grating_and_chop'
        return, -1
endelse

;---------------------------------------------------------------------
;  Sanity check: 
;  if C and ramplength are not an integer multiple of each other, exit 
;---------------------------------------------------------------------
if (C ge ramplength) and ~(fix(C/ramplength)*ramplength eq C) then begin
  print, 'Choplength not an integer multiple of ramplength.  Exiting split_cycle.pro.'
  return, -1
endif

if (ramplength ge C) and ~(fix(ramplength/C)*C eq ramplength) then begin
  print, 'Ramplength not an integer multiple of choplength.  Exiting split_cycle.pro.'
  return, -1
endif

;------------------------------------------------------
; remove any partial ramps @ beginning and end of file
; 24oct2014
;------------------------------------------------------
i = 0L
first_ramp = (data[0].header)[P_RMP_CNT]
while (data[i]).header[P_SMPL_CNT] ne 0 do i = i + 1
head = i
anz_records = n_elements(data)
last_ramp = (data[anz_records-1]).header[P_RMP_CNT]
last_frame = (data[anz_records-1]).header[P_SMPL_CNT]
if last_frame eq ramplength-1 then data = data[head:*] $ ; remove 1st ramp only
else begin
while (data[i]).header[P_RMP_CNT] lt last_ramp do i = i + 1
  tail = i
  data = data[head:tail]     ; remove 1st and last partial ramps
endelse

;-----------------------------------------------------------------------
; trim off any "unpaired" chop plateaus - 4POINT chop is UNTESTED!
; 27oct2014
; assumption: 1. file starts with chop0
;             2. chop0 and a following chop1 (2POINT) / chop 2 (4POINT)
;                makes a pair (chop1 and chop3 makes another pair for 
;                4POINT chop scheme)
;             3. no missing ramps except at beginning or
;                end of file
;-----------------------------------------------------------------------
; # of ramps per chop phase - or # of ramps to co-add in fit_ramps
n =  fix(C/ramplength) 
if n gt 0 then begin
  index = (data.header)[P_RMP_CNT, *]/n
  chop_scheme = fix(strmid(C_SCHEME, 0, 1))   ; 2 or 4
  chop_phase = index mod chop_scheme   ; chop phase of each readout 
  chop_phase_max = max(chop_phase)
  head = 0
  while chop_phase[head] ne 0 do head +=  1 
  tail = 1L
  while chop_phase[-tail] ne chop_phase_max do tail +=  1
  data = data[head:-tail] 
endif

anz_records = n_elements(data)        ; # of readouts
nramps = anz_records/ramplength       ; # of ramps

nramps_per_ind_pos = nramps / (G_PSUP+G_PSDN)

nchops = nramps / (CHOP_LN/ramplength)
nchops_per_ind_pos = nchops / (G_PSUP+G_PSDN)

if verbose then begin
   print, '# of readouts = ', anz_records
   print, 'ramplength = ', ramplength
   print, 'total number of ramps = ', nramps
   print, 'number of ramps per grating position = ', nramps_per_ind_pos
   print, 'number of chops per grating position = ', nchops_per_ind_pos 
endif

; from demodulate
new_pos = G_STRT+indgen(G_PSUP)*G_SZUP

n = C/ramplength

;------------------------------------------
;-- separate the chop cycles: 2POINT chop
;------------------------------------------
if strtrim(C_SCHEME, 2) eq '2POINT' then begin 
  ;---- case 1: ramplength = C, one ramp per chop; does this happen?
   if (ramplength eq C) then begin
      ; why is this different than case 3?
    q0 = where((data.header)[P_RMP_CNT, *] mod 2 eq 0)
    q1 = where((data.header)[P_RMP_CNT, *] mod 2 eq 1)
    if verbose then print, '2POINT chop, case 1: ramplength = CHOP_LN'
  endif

  ;---- case 2: ramplength = 2*n * C, n = 1, 2, 3, ...
  if (ramplength gt C) and (ramplength/C mod 2 eq 0) then begin
    chop = fix((data.header)[P_SMPL_CNT]/C) mod 2	
    q0 = where(chop eq 0) 
    q1 = where(chop eq 1) 
    if verbose then print, '2POINT chop, case 2: ramplength = 2*n*CHOP_LN where n = 1, 2, 3, I do not handle this case...'
  endif

  ;---- case 3: C = n * ramplength, n =  2, 3, 4,  ...
  if (ramplength lt C) then begin
     n = C/ramplength                ; # of ramps per chop phase     
     
    ;chop = fix(indgen(nramps)/n) mod 2	; chop phase @ each ramp	
    chop = fix((data.header)[P_RMP_CNT, *]/n) mod 2 
    q0 = where(chop eq 0)		; ramp number in chop phase 0 
    q1 = where(chop eq 1)		; ramp number in chop phase 1
                                
    size_chop_grat_bin = nchops_per_ind_pos / n
    if verbose then print, 'size of grating chop bin = ', size_chop_grat_bin
    ; for each size_chop_grat bin in the final split up chop0 and chop1 split 
    ;
    G_POS = G_STRT + indgen(G_PSUP)*G_SZUP

    ind_array0 = intarr(G_PSUP, size_chop_grat_bin*CHOP_LN)
    ind_array1 = intarr(G_PSUP, size_chop_grat_bin*CHOP_LN)

    if file_test(fitsfile0) gt 0 then file_delete, fitsfile0
    if file_test(fitsfile1) gt 0 then file_delete, fitsfile1
    
    fxaddpar, primehead, 'CHOPNUM', 0, " Chop number"
    fxwrite, fitsfile0, primehead, /noupdate, /append
    fxaddpar, primehead, 'CHOPNUM', 1, " Chop number"    
    fxwrite, fitsfile1, primehead, /noupdate, /append

    nodpos = fxpar(primehead, 'NODPOS')
    
    i = 0

    while i lt G_PSUP do begin
       ind_array0[i, *] = q0[i*size_chop_grat_bin*CHOP_LN:(i+1)*(size_chop_grat_bin*CHOP_LN)-1]
       ind_array1[i, *] = q1[i*size_chop_grat_bin*CHOP_LN:(i+1)*(size_chop_grat_bin*CHOP_LN)-1]
    
       fitsfile0 = outdir+'/'+base_filename+'chop0.fits'
       fitsfile1 = outdir+'/'+base_filename+'chop1.fits'
       
       mkhdr, exthdr, data[ind_array0[i, *]].data
       
       fxaddpar, exthdr, 'INDPOS', new_pos[i], " Inductosyn position "
       fxaddpar, exthdr, 'CHOPNUM', 0, " Chop number"
       fxaddpar, exthdr, 'NODPOS', nodpos, " Nod position"

       mwrfits, data[ind_array0[i, *]].data, fitsfile0, exthdr, /silent
       
       fxaddpar, exthdr, 'CHOPNUM', 1, "Chop number"
       mwrfits, data[ind_array1[i, *]].data, fitsfile1, exthdr, /silent
       i += 1
    endwhile
    
   
 endif

  data0 = data[q0]

  result  = 1

endif

return, result

end

