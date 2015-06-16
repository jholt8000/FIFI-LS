;+
; NAME: fifi_ls_offset_xy 
;
; PURPOSE: return each module's (x,y) positions from date of observation
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: result=fifi_ls_offset_xy(DATE,BLUE=BLUE)
;
; INPUTS: DATE - date of observation input as a vector, ex. [2009,3,23] 
;
; OPTIONAL INPUTS:
;
; KEYWORDS: BLUE - specify blue channel (RED is assumed if keyword not set)
;
; OUTPUTS: result - a 2-D array (2x25) containing the (x,y) offsets of each 
; module (currently in mm). 
;
; OPTIONAL OUTPUTS:
;
; SIDE EFFECTS: 
;
; RESTRICTIONS: 
;
; PROCEDURE:
;
; EXAMPLE: result=offset_xy([2009,06,17],/red)
;
; MODIFICATION HISTORY:
; 14apr09  k.n.  started coding...
;
; $Id: offset_xy.pro, v 1.4 2014-03-24 23:55:49 klein Exp $
;-
function fifi_ls_offset_xy, date, blue = blue 

;-------------
;  help file 
;-------------
if n_params() lt 1 then begin
  doc_library, 'fifi_ls_offset_xy'
  return, -1
endif

;----------------
; sanity checks
;----------------
; RED or BLUE?

if ~keyword_set(BLUE) then begin 
;  print, 'offset_xy.pro: BLUE keyword not set; RED channel assumed.'
  channel = 'R'
endif else begin
  channel = 'B'
endelse

; is DATE a 3-integer vector and reasonable?
if size(date, /n_elements) ne 3 then begin
  print, 'offset_xy.pro: DATE must be a 3-element vector.  Exiting.'
  return, -1
endif

;------------------------------------------
; Select the appropriate calibration file
;------------------------------------------
; convert date of calibration and observation to Reduced Julian Dates
orgdate = date		; save the 3-element array
juldate, date, obsdate	; date is replaced by a 6-element array
date = orgdate		; restore the original array

;------------------------------------------------
; define location of calibration parameter files
;------------------------------------------------
 fifi_ls_path = file_dirname(file_dirname(file_which( $
                 'fifi_ls_split_grating_and_chop.pro')), /MARK)
 fifi_ls_data = filepath('data', ROOT_DIR = fifi_ls_path) + path_sep()

calfilename_r = fifi_ls_data+['poscal_20081212_r.txt', 'poscal_20140201_r.txt']	
;calfilename_r = fifi_ls_data+['poscal_20081212_r.txt', 'poscal_20091212_r.txt'] 
calfilename_b = fifi_ls_data+['poscal_20120401_b.txt', 'poscal_20140201_b.txt']
; BLUE the first blue file is just a dummy
;caldate_r = [[2008, 12, 12], [2009, 12, 12]]
caldate_r = [[2008, 12, 12], [2014, 02, 01]]
caldate_b = [[2012, 04, 01], [2014, 02, 01]]

case channel of
  'B': begin
       calfilename = calfilename_b & caldate = caldate_b
       end
  'R': begin
       calfilename = calfilename_r & caldate = caldate_r
       end
endcase

jd = dblarr(n_elements(caldate)/3)
for i = 0, n_elements(caldate)/3-1 do begin
  juldate, caldate[*, i], temp
  jd[i] = temp
endfor

; select the appropriate (closest by date) calibration file
;
;diff = abs(obsdate-jd)
;q = where(abs(obsdate-jd) eq min(abs(obsdate-jd)))

; select the appropriate (latest before the obs date) calibration file
q = max(where((jd ge obsdate) eq 0)) > 0
; if obsdate is before the first cal.  file, then use the first one. 

; read calibration file
;print, 'offset_xy.pro: calibration file '+calfilename[q]+' read.'
readcol, calfilename[q], x, y, /silent

result = fltarr(2, 25)
result[0, *] = x
result[1, *] = y


;STOP
return, result 

end
