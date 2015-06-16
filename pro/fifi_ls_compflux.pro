
; MODIFICATION HISTORY:
; 2015apr20  bvacca: first version

function fifi_ls_compflux, wobs, fobs, ERR=eobs, $
                           master_grid=master_grid, $
                   FITORDR=fitordr, FACTR=factr, SAMPL=sampl, $
                   THRESH=thresh, FWHM=fwhm, NEGCUT=negcut, NEGTHRESH=negthresh
                           
if ~keyword_set(FITORDR) then fitordr = 5	; Fitting order
if ~keyword_set(FACTR) then factr = 2.0		; Size of region considered for fitting: w-FWHM/factr to w+FWHM/factr
if ~keyword_set(SAMPL) then sampl = 2.5		; sampling interval: FWHM/sampl
if ~keyword_set(FWHM)  then fwhm  = 0.0478	; FWHM in microns
if ~keyword_set(THRESH) then thresh = 5.0	; robust threshold
if ~keyword_set(NEGTHRESH) then negthresh=thresh	; threshold for identifying negative values

eps = 0.01
pixscale = 0.02
novsamp = fwhm * sampl / pixscale

findx = where(finite(fobs))
wtmp = wobs[findx]
etmp = eobs[findx]
ftmp = fobs[findx]
wobs = wtmp
eobs = etmp
fobs = ftmp

windx = sort(wobs)
fobs = fobs[sort(wobs)]
eobs = eobs[sort(wobs)]
wobs = wobs[sort(wobs)]

;wobs = wtmp
;fobs = ftmp
;eobs = etmp
nobs = n_elements(wobs)

if keyword_set(NEGCUT) then begin
   print,'in negcut'
   robuststats, fobs, thresh, fmean, fvar, fstddev, /SILENT
   ibad = where(fobs lt -1.0*negthresh*fstddev, COMPLEMENT=ingood)
   wtmp = wobs[ingood]
   ftmp = fobs[ingood]
   etmp = eobs[ingood]
   wobs = wtmp
   fobs = ftmp
   eobs = etmp
   nobs = n_elements(wobs)
endif

;farr = master_grid


;;;

;wmin = min(wave)
;wmax = max(wave)

;dw   = fwhm/novsamp
;nw   = (wmax - wmin)/dw + 1
;print, 'DW = ',dw, fwhm/sampl

;warr = findgen(nw)*dw + wmin
;warr = findgen(nw)*dw + wmin + (dw/2.0)
;farr = fltarr(nw)
;ferr = fltarr(nw)
;master_grid = warr
;;;


delta = fwhm/factr
;print, 'Delta = ',delta,' Novsamp = ',novsamp
rms = 0.0

warr = master_grid
nw = n_elements(warr)
farr = fltarr(nw)
ferr = fltarr(nw)

 for i = 0,nw-1 do begin

    indx = where((wobs ge warr[i]-delta) and (wobs le warr[i]+delta), ncnt)
    if (i gt 0) and isarray(savindx) then begin
       idum = where(indx eq savindx, icnt)
       if (icnt eq ncnt) then print, 'Warning: FACTR too large'
    endif
    win  = wobs[indx]
    fin  = fobs[indx]
    ein  = eobs[indx]


    if (ncnt ge 1) then begin
       if (ncnt gt 1) then begin
           robuststats, fin, thresh, fmean, fvar, fstd, fskew, fkurt, OGOODBAD=ogood, /SILENT
           nigood = n_elements(ogood)
           igood  = where(ogood eq 1)
       endif else begin
           igood = indx
           nigood = 1
       endelse

;    print, 'Numbers ',i, ncnt, nigood

       if (nigood ge 1) then begin 
         if (nigood ge fitordr+1) then order = fitordr else order = nigood - 1
         if keyword_set(ERR) then begin
;              coefs = robustpoly1d(win[igood], fin[igood], order, thresh, eps, YERR=ein[igood], RMS=rms, /SILENT)
             coefs = poly_fit1d(win[igood], fin[igood], order, YERR=ein[igood], RMS=rms, /SILENT)
         endif else begin
;              coefs = robustpoly1d(win[igood], fin[igood], order, thresh, eps, RMS=rms, /SILENT)
             coefs = poly_fit1d(win[igood], fin[igood], order, RMS=rms, /SILENT)
         endelse

 	 if ((warr[i] lt min(win[igood])) or (warr[i] gt max(win[igood]))) then farr[i] = mean(fin[igood]) else farr[i] = poly1d(warr[i],coefs)
         ferr[i] = rms

;    print, min(win),max(win)
       endif else begin
          print, 'Number of good points eq 0'
       endelse 
       savindx = indx
    endif else begin
       print, 'ncnt = ', ncnt
       print, ' i = ', i
       print, 'Number of good points eq 0'
  endelse 

endfor
  
   
;window, 3

;print, 'err=', err
;plot, wave, flux, psym=4, xtitle='Wavelength (microns)', ytitle='Flux'
;if keyword_set(ERR) then oploterr, wave, flux, err
;oplot, wobs, fobs, psym=4, color=100
;oplot, wobs, ffilt, psym=4, color=5
;oplot, master_grid, farr, psym=4, color=200
;oplot, master_grid, farr, psym=10, color=200
;errplot, master_grid, farr-ferr,farr+ferr, color=200

;oplot, wout, fout, psym=10, color=5

;plots, [!X.CRANGE[0],!X.CRANGE[1]], [0.0,0.0], linestyle=2
;wait, 0.5
return, farr
;return, fout

end

