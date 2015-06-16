

PRO fifi_smoothing_tests, hypercube

  lambdas=reform(hypercube[2,2,1,*])
  fluxes = reform(hypercube[2,2,0,*])

  slambdas = lambdas[sort(lambdas)]
  sfluxes = fluxes[sort(lambdas)]
  
  x2=slambdas
  
  y2=sfluxes
   
  ; Sampling interval:
  dl = 0.01

     
  ; Display the first plots and define the layout
  p1 = PLOT(x2, y2, NAME='Signal Noise', 'k4*', LAYOUT=[1,2,1], DIMENSIONS=[1000,800])
   
  p2 = PLOT(x2, SMOOTH(y2, 10, /EDGE_TRUNCATE), /OVERPLOT, $
     COLOR=[255, 0, 0], NAME='Smooth (width 10)')
  savgolFilter = SAVGOL(3, 3, 0, 4)
   
  p3 = PLOT(x2, CONVOL(y2, savgolFilter, /EDGE_TRUNCATE), /OVERPLOT, $
     COLOR=[0, 0, 255], $
            NAME='Savitzky-Golay (width 6, 4th degree)')

  savgolFilter2 = SAVGOL(5,5,0,4)
  p33 = PLOT(x2, CONVOL(y2, savgolFilter2, /EDGE_TRUNCATE), /OVERPLOT, $
     COLOR=[0, 255, 0], $
            NAME='Savitzky-Golay (width 10, 4th degree)')

  nlambdas = (max(x2) - min(x2)) / dl
  lgrid = findgen(nlambdas+1)*dl + min(x2)

  p4 = PLOT(x2, y2, 'k4*', layout=[1,2,2], /current)
  p5 = PLOT(lgrid, interpol(y2,x2,lgrid, /nan), '*', $
            COLOR=[250,100,50], NAME='interpol linear',/overplot)
  ;p6 = PLOT(lgrid, interpol(y2,x2,lgrid,/spline, /nan), 'S', $
  ;          COLOR=[0,250,50], NAME='interpol spline',/overplot)
  p7= PLOT(lgrid, interpol(y2,x2,lgrid,/lsquadratic, /nan), 'D', $
            COLOR=[0,0,250], NAME='interpol lsquadratic',/overplot)

  
   
  l1 = LEGEND(TARGET=[p1, p2, p3,p33,p5,p7], $
     /DEVICE, /AUTO_TEXT_COLOR,position=[400,450])
   
END

pro poly_interp, hypercube

  lambdas=reform(hypercube[2,2,1,*])
  fluxes = reform(hypercube[2,2,0,*])

  x = lambdas[sort(lambdas)]
  y = fluxes[sort(lambdas)]
  delta_lambda=0.02

  nlambdas=(max(x) - min(x)) / delta_lambda

  lambda_grid = findgen(nlambdas+1)*delta_lambda + min(x)

  new_fluxes=dblarr(nlambdas)
  
  for i=0, nlambdas-1 do begin
     polint, x, y, lambda_grid[i], new_fluxes[i], delta_lambda
     print, new_fluxes[i]
     print, lambda_grid[i]
  endfor

  p1 = plot(x,y, 'k4*')
  p2 = plot(lambda_grid, new_fluxes,/overplot)
end

pro fst, hypercube

  lambdas=reform(hypercube[2,2,1,*])
  fluxes = reform(hypercube[2,2,0,*])

  slambdas = lambdas[sort(lambdas)]
  sfluxes = fluxes[sort(lambdas)]
  dl=0.01
  x2=slambdas
  
  y2=sfluxes
  
  ; Display the first plots and define the layout
  ;p1 = PLOT(x2, y2, NAME='Signal Noise', 'k4*', LAYOUT=[1,2,1], DIMENSIONS=[1000,800])

  sm33 = SMOOTH(y2, 33, /EDGE_TRUNCATE)
  p2 = PLOT(y2-sm33,  COLOR=[255, 0, 0], NAME='Smooth (width 33)', $
           layout=[1,3,1], dimensions=[1000,1000])

  savgolFilter = SAVGOL(3, 3, 0, 4)
  sg44=convol(y2, savgolFilter, /edge_truncate)
  p3 = PLOT(y2-sg44, /OVERPLOT, $
     COLOR=[0, 0, 255], $
            NAME='Savitzky-Golay (width 4, 4th degree)')

  savgolFilter2 = SAVGOL(5,5,0,4)
  sg104 = CONVOL(y2, savgolFilter2, /EDGE_TRUNCATE)
  p33 = PLOT(y2-sg104, /OVERPLOT, COLOR=[0, 255, 0], $
            NAME='Savitzky-Golay (width 10, 4th degree)')

  nlambdas = (max(x2) - min(x2)) / dl
  lgrid = findgen(nlambdas+1)*dl + min(x2)

 ; p4 = PLOT(y2-x2, 'k4*', layout=[1,2,2], /current)
  ilin=interpol(y2,x2,lgrid, /nan)
  print, 'ilin=',ilin
  print, 'y2=',y2
  print,'lgrid=',lgrid
print,'x2=',x2
  p5 = PLOT(y2-ilin[x2], '*', $
            COLOR=[250,100,50], NAME='interpol linear', $
            layout=[1,3,2],/current)

  ilsquad= interpol(y2,x2,lgrid,/lsquadratic, /nan)
  p7= PLOT(y2- ilsquad[x2], 'D', $
            COLOR=[0,0,250], NAME='interpol lsquadratic',/overplot)


  
  l1 = LEGEND(TARGET=[p2, p3,p33,p5,p7], $
      /DEVICE, /AUTO_TEXT_COLOR,position=[400,450])



  p6=plot(x2,y2,layout=[1,3,3],'r4D',/current)
END

pro gaussplot, hypercube

  lambdas=reform(hypercube[2,2,1,*])
  fluxes = reform(hypercube[2,2,0,*])

  slambdas = lambdas[sort(lambdas)]
  sfluxes = fluxes[sort(lambdas)]
  
  x=slambdas

  lambda1 = min(x)
  print,lambda1
  y=sfluxes

  yfit = gaussfit(x, y, coeff, NTERMS=3)

  print,'coeff=',coeff
  sigma = coeff[2]
  fwhm = 2 * sqrt(2*alog(2)) * sigma
  delta_lambda = fwhm/2.
  print,'sigma=',sigma
  I0 = 1./(sqrt(2.*!pi*sigma^2))

  print,'I0=',I0

  p1=plot(x,y,'rd*', title='sigma='+strtrim(sigma,2)+' delta lambda = '+strtrim(delta_lambda,2))
  p2=plot(x,yfit,thick=2,/overplot)

  topj = dblarr(n_elements(y))
  Pd=0
  yy2 = dblarr(n_elements(y))

  nlambdas=(max(x) - min(x)) / delta_lambda

  lambda_grid = findgen(nlambdas+1)*delta_lambda + min(x)

  print, qromb('gauss', lambda1, lambda1+delta_lambda)
  
  ;for i=0,n_elements(x)-1 do begin
     
     ;P = I0 * exp(-(x-x[i])^2 / (2*sigma^2))

     ;print,'l1=', (i*delta_lambda)+lambda1, '  l2=',lambda1+(i+1)*delta_lambda
     ;p_discrete = qromb('gauss', (i*delta_lambda)+lambda1, lambda1+(i+1)*delta_lambda)

     ;print,'p_discrete=',p_discrete,' i=',i
     
     ;yy2[i] = y[i] / p_discrete
     
     ;Pd = Pd + I0 * exp(-(x-x[i])^2 / (2*sigma^2))
     ;topj[i] = y[i] / P
     
     ;plot,P
     ;oplot,y,x-coeff[1],psym=3
     ;wait,0.01
   ;endfor

  ;p3=plot(x, yy2, 'rd*')
end

function gauss, x
  return, 19.62 * exp(-( (x-63.170383)^2) / (2*0.020332663^2))
end


pro nyq, hypercube
  lambdas=reform(hypercube[2,2,1,*])
  fluxes = reform(hypercube[2,2,0,*])

  slambdas = lambdas[sort(lambdas)]
  sfluxes = fluxes[sort(lambdas)]
  
  x=slambdas
  
  y=sfluxes
  
  plot,x,y,psym=10
end

pro deconv_psf, hypercube
  lambdas=reform(hypercube[2,2,1,*])
  fluxes = reform(hypercube[2,2,0,*])

  slambdas = lambdas[sort(lambdas)]
  sfluxes = fluxes[sort(lambdas)]
    
end
