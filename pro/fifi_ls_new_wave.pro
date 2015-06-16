
; NAME: fifi_ls_wave 
;
; PURPOSE: calculate wavelength from inductosyn position and date of
; observation 
;
; CATEGORY: FIFI LS data reduction
;
; CALLING SEQUENCE: result=fifi_ls_wave(IND,DATE,BLUE=BLUE)
;
; INPUTS: IND - inductosyn position
;         DATE - date of observation input as a vector, ex. [2009,3,23]
;
; OPTIONAL INPUTS:
;
; KEYWORDS: BLUE - specify blue channel (RED is assumed if keyword not set)
;
; OUTPUTS: result - a 2-D array (25x16) containing the wavelength of each pixel 
; for each module in microns 
;
; OPTIONAL OUTPUTS:
;
; SIDE EFFECTS: the dummy/open spectral channels will be removed 
;
; RESTRICTIONS: 
;
; PROCEDURE:
;
; EXAMPLE: result=fifi_ls_wave(100000.,[2009,06,17],/red)
;
; MODIFICATION HISTORY:
; 24mar09  k.n.  started coding...
;
; 03mar14  k.n.  introducing new calibration file format - one common 
;                (spectral channel independent) and three channel-dependent
; 
; $Id: wave.pro,v 1.7 2014-03-03 16:26:30 kaori Exp $
;-
function fifi_ls_new_wave, ind, date, blue=blue 

;-------------
;  help file 
;-------------
if n_params() lt 2 then begin
  doc_library,'wave'
  return, -1
endif

;----------------
; sanity checks
;----------------
; RED or BLUE?

if ~keyword_set(BLUE) then begin
;  print,'fifi_ls_wave: BLUE keyword not set; RED channel assumed.'
  channel='R' 
endif else begin
  if (strupcase(blue) ne 'B1') AND (strupcase(blue) ne 'B2') then begin
    print, 'fifi_ls_wave: BLUE keyword must be B1 or B2. Exiting.'
    return, -1
  endif else channel=strupcase(blue)
endelse

; is IND a reasonable number?
if ind lt 0 then begin		; what's the min/max IND can be?
  print, 'fifi_ls_wave: IND must be positive.  Exiting.'	
  return, -1
endif

; is DATE a 3-integer vector and reasonable?
if size(date, /n_elements) ne 3 then begin
  print,'fifi_ls_wave: DATE must be a 3-element vector.  Exiting.'
  return, -1
endif

;------------------------------------------------
; define location of calibration parameter files
;------------------------------------------------
 fifi_ls_path = file_dirname(file_dirname(file_which( $
                 'fifi_ls_split_grating_and_chop.pro')), /MARK)
 fifi_ls_data = filepath('data', ROOT_DIR=fifi_ls_path) + path_sep()

;--------------------------------------------------------
; after 2014-02-01: select appropriate calibration files 
; -------------------------------------------------------
;orgdate=date        ; save the 3-element array
;juldate,date,obsdate    ; date is replaced by a 6-element array
;date=orgdate        ; restore the original array

;juldate, [2015,01,01], transdate
obsdate = date
; three dates are 20150101, 20140201, and before

if obsdate[0] ge 2014 then begin

   calfilename = fifi_ls_data + ['wavecal_20150420.txt']

   case channel of
      'R': begin
         isofffilename = fifi_ls_data + ['wavecalISUOFF_20150420_r.tsv']
         cn = 0
         m = 1
         end
      'B1': begin
         isofffilename = fifi_ls_data + ['wavecalISUOFF_20150420_b1.tsv']
         cn = 1
         m = 1         
         end
      'B2': begin
         isofffilename = fifi_ls_data + ['wavecalISUOFF_20150420_b2.tsv']
         cn = 2
         m = 2
         end
    endcase

   readcol, calfilename, g0, GOFF, a, ISF, gama, PS, POFF, scalecorr, $
         format='d, d, d, d, d, d, d, d', /silent
   readcol, isofffilename, module, ISOFF, format='I, D', /silent
                                ; module vs ISOFF

   pix = indgen(16) + 1     ; spectral pixel number, 1 to 16
   result = dblarr(25,16)

   ; here's the equation
   for module = 0, 24 do begin
     
     phi = 2d*!pi*ISF[cn]*((double(ind)+ISOFF[module])/2d^24)

     sign = (pix - POFF[cn])/abs((pix - POFF[cn])) ; +1 or -1
     delta = (pix - 8.5) * PS[cn] + sign * (pix - POFF[cn])^2 * scalecorr[cn]
     
     SlitPos = 25 - 6 * fix(module/5) + (module mod 5) ; module = 0 to 24
    
     g = g0[cn] * cos(atan((slitPos - GOFF[cn])/a[cn]))
     
     lambda = 1000. * (g/m) *(sin(phi+gama[cn]+delta) + sin(phi-gama[cn]))     
     result[module,*] = lambda
  endfor
  
endif else if obsdate[0] gt 2013 then begin
   calfilename = fifi_ls_data + ['wavecal_20140201.txt']  
   isoff_r = fifi_ls_data + ['wavecalISUOFF_20140201_r.tsv']
   isoff_b1 = fifi_ls_data + ['wavecalISUOFF_20140201_b1.tsv']
   isoff_b2 = fifi_ls_data + ['wavecalISUOFF_20140201_b2.tsv']

   readcol, calfilename, ga, gb, gc, ISF, gama, PS, POFF, scalecorr, $
         format='d, d, d, d, d, d, d, d', /silent

   case channel of
  'R': begin
     isofffilename = fifi_ls_data + ['wavecalISUOFF_20140201_r.tsv']
     cn = 0  &  m=1
     end
  'B1': begin
     isofffilename  = fifi_ls_data + ['wavecalISUOFF_20140201_b1.tsv']
     cn = 1  &  m=1
     end
  'B2': begin
     isofffilename = fifi_ls_data + ['wavecalISUOFF_20140201_b2.tsv']
     cn = 2  &  m=2
     end
   endcase

   ; each variable contains values for R, B1, and B2
   readcol, isofffilename, module, ISOFF, format='I, D',/silent  

   pix = indgen(16) + 1     ; spectral pixel number, 1 to 16
   result = dblarr(25, 16)

   ; here's the equation
   for module = 0, 24 do begin

     ;print, ISF[cn], ISOFF[module], pix, POFF[cn], PS[cn], scalecorr[cn]
     
     phi = 2*!pi*ISF[cn]*((ind+ISOFF[module])/2d^24)
     sign = (pix - POFF[cn])/abs((pix - POFF[cn])) ; +1 or -1
     delta = (pix - 8.5) * PS[cn] + sign * (pix - POFF[cn])^2 * scalecorr[cn]

     ;print, phi, sign, delta
     
     SlitPos = 25 - 6 * fix(module/5) + (module mod 5) ; module = 0 to 24
     ; SlitPos = 25 26 27 28 29 19 20 21 22
     ; 23 13 14 15 16 17 7 8 9 10 11 1 2 3
     ; 4 5 from slit geometery, pixels are really spaxel number -jh
  
     g = ga[cn] * SlitPos^2 + gb[cn] * SlitPos + gc[cn]

     lambda = 1000. * (g/m) *(sin(phi+gama[cn]+delta) + sin(phi-gama[cn]))
     result[module,*] = lambda
  endfor
   
endif else begin
 
   ;-------------------------------------------------------
   ; Pre-20140201: select the appropriate calibration file
   ;-------------------------------------------------------
   ; convert date of calibration and observation to Reduced Julian Dates
   orgdate=date                 ; save the 3-element array
   juldate,date,obsdate         ; date is replaced by a 6-element array
   date=orgdate                 ; restore the original array

   calfilename_r=fifi_ls_data+['wavecal_20081212_r.txt','wavecal_20120201_r.txt'] ; RED  
   calfilename_b=fifi_ls_data+['wavecal_20120401_b.txt'] ; BLUE dummy calibration file; RK 2012-07-02
   caldate_r=[[2008,12,12],[2012,02,01]]
   caldate_b=[[2012,04,01]]     ;BLUE dummy calibration file; RK 2012-07-02

   case channel of
      'B1': begin         ;--- before 20140201 there is no difference between B1 and B2
         calfilename=calfilename_b & caldate=caldate_b
      end
      'B2': begin
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

   ; select the appropriate (latest before the obs date) calibration file
   q=max(where((jd ge obsdate) eq 0)) > 0 ; if obsdate is before the first cal. 
                                ; file, then use the first one. 

   ; read calibration file
   readcol,calfilename,g,off,pixscal,scalcorr2,scalcorr3,format='d,l,d,d,d',/silent

   ; here's the equation
   result=dblarr(25,16)
   m=1
   gama=1.2*!pi/180
   pix=indgen(16)+1

   for module = 0, 24 do begin
      phi=(ind+off[module])/2d^16*360./256.*!pi/180.
      delta=(pix-8.5)*pixscal[module]+(pix-8.5)^2*scalcorr2[module]+ $
            (pix-8.5)^3*scalcorr3[module]
      lambda=1000.*g[module]/m*(sin(phi+gama+delta)+sin(phi-gama))

      result[module,*]=lambda
   endfor

endelse

return,result 

end
