;+
; NAME: fifi_ls_bad_pixels 
;
; PURPOSE: return dead pixels [module #, channel #] from date of observation
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: result=bad_pixels(DATE,BLUE=BLUE)
;
; INPUTS: DATE - date of observation input as a vector, ex. [2009,3,23] 
;
; OPTIONAL INPUTS:
;
; KEYWORDS: BLUE - specify blue channel (RED is assumed if keyword not set)
;
; OUTPUTS: result - a 2-D array (2 x # of lines in bad pixel file) 
; containing the (module,channel) pair of each dead pixel.
;
; OPTIONAL OUTPUTS:
;
; SIDE EFFECTS: 
;
; RESTRICTIONS: 
;
; PROCEDURE:
;
; EXAMPLE: result=fifi_ls_bad_pixels([2009,06,17],/blue)
;
; MODIFICATION HISTORY:
; 15apr09  k.n.  started coding...
;
; 05feb14  k.n.  added reference to bad pixel files from 20131202 
;                lab calibration
; 
; $Id: bad_pixels.pro,v 1.5 2014-02-05 22:46:41 kaori Exp $
;-
function fifi_ls_bad_pixels, date, blue=blue 

;-------------
;  help file 
;-------------
if n_params() lt 1 then begin
  doc_library,'bad_pixels'
  return, -1
endif

;----------------
; sanity checks
;----------------
; RED or BLUE?

if ~keyword_set(BLUE) then begin 
;  print,'bad_pixels.pro: BLUE keyword not set; RED channel assumed.'
  channel='R'
endif else begin
  channel='B'
endelse

; is DATE a 3-integer vector and reasonable?
if size(date,/n_elements) ne 3 then begin
  print,'bad_pixels.pro: DATE must be a 3-element vector.  Exiting.'
  return, -1
endif

;----------------------
;
;----------------------
; convert date of calibration and observation to Reduced Julian Dates
orgdate=date		; save original 3-element array
juldate,date,obsdate	; date gets replaced by a 6-element array
date=orgdate		; restore the original array

;------------------------------------------------
; define location of calibration parameter files
;------------------------------------------------
 fifi_ls_path = file_dirname(file_dirname(file_which( $
                 'fifi_ls_split_grating_and_chop.pro')), /MARK)
 fifi_ls_data = filepath('data', ROOT_DIR=fifi_ls_path) + path_sep()

calfilename_r=fifi_ls_data+['badpixels_20081212_r.txt','badpixels_20120613_r.txt', $
             'badpixels_20131202_r.txt']	; RED
;calfilename_r=fifi_ls_data+['poscal_20081212_r.txt','poscal_20091212_r.txt']	; RED  
calfilename_b=fifi_ls_data+['badpixels_20120201_b.txt','badpixels_20131202_b.txt']	; BLUE
;caldate_r=[[2008,12,12],[2009,12,12]]
caldate_r=[[2008,12,12],[2012,06,13],[2013,12,02]]
caldate_b=[[2012,02,01],[2013,12,02]]; BLUE 

case channel of
  'B': begin
       calfilename=calfilename_b & caldate=caldate_b
       end
  'R': begin
       calfilename=calfilename_r & caldate=caldate_r
       end
endcase

jd=dblarr(n_elements(caldate)/3)
for i=0,n_elements(caldate)/3-1 do begin
  juldate,caldate[*,i],temp
  jd[i]=temp
endfor

; select the appropriate (closest by date) calibration file
;
;diff=abs(obsdate-jd)
;q=where(abs(obsdate-jd) eq min(abs(obsdate-jd)))

; select the appropriate (latest before the obs date) calibration file
q=max(where((jd ge obsdate) eq 0)) > 0	; if obsdate is before the first cal. 
				 	; file, then use the first one. 

; read calibration file
;print, 'bad_pixels.pro: calibration file '+calfilename[q]+' read.'
readcol,calfilename[q],x,y,format='I,I',/silent

result=intarr(2,n_elements(x))
result[0,*]=x
result[1,*]=y


;STOP
return,result 

end
