           


           ; case B
           ; case B2
           ; case B3

           ; case C
           ; case D
           

           
           ; get the left hand side
           if ((((wtd_x[j] - deltax) / 2.) lt x1) and (wtd_x[j] gt x0)) then begin
              ; first case if point is on the left of the flux
              if (wtd_x[j] - x0) lt deltax / 2. then begin
                 Lx1[j] = wtd_x[j] - x0
              ; point is on the right of the flux
              endif else if (wtd_x[j] - x1) lt deltax / 2. then begin
                 Lx1[j] = x1 - ( wtd_x[j] - deltax /2.)
              endif else begin
                 Lx1[j] = deltax / 2.
              endelse             
           endif else Lx1[j] = 0
           
           ; get the right hand side
           if (((wtd_x[j] + deltax / 2.) gt x0) and $
               (wtd_x[j] lt x1)) then begin
              if (x1 - wtd_x[j]) lt (deltax / 2.) then begin                 
                 Lx2[j] = x1 - wtd_x[j]
              endif else begin
                 Lx2[j] = deltax /2.
                 print,'deltax/2.=',deltax / 2.
              endelse              
           endif else Lx2[j] = 0

        endfor
     endif

     y0 = yes[spax] - deltay / 2.
     y1 = yes[spax] + deltay / 2.
     ycenter = where(y_grid ge y0 and y_grid le y1, count)
  
     if count ne 0 then begin
        if ycenter[0] lt 1 then begin
           ystart = 0
        endif else begin
           ystart = ycenter[0]-1
        endelse
        if ycenter[-1]+1 gt n_elements(y_grid) then begin
           yend = ycenter[-1]
        endif else begin
           yend = ycenter[-1] + 1
        endelse  
          
        wtd_y = y_grid[ystart:yend]
        Ly1 = dblarr(n_elements(wtd_y))
        Ly2 = dblarr(n_elements(wtd_y))
        print,'wtd_y = ', wtd_y
        print, 'y_grid=',y_grid
        for j=0, n_elements(wtd_y)-1 do begin
           
           ; get the left hand side
           print,'wtyd=',wtd_y[j],'dy=',deltay,'y1=',y1,'y0=',y0
           if (((wtd_y[j] - (deltay / 2.)) lt y1) and $
               (wtd_y[j] gt y0)) then begin
              ; 
              if (wtd_y[j] - y0) lt deltay / 2. then begin
                 Ly1[j] = wtd_y[j] - y0
              ; 
              endif else if ((wtd_y[j] - y1) lt (deltay / 2.)) then begin
                 Ly1[j] = y1 - ( wtd_y[j] - deltay /2.)
                 print,'here'
                 print, 'ly1[j]=',Ly1[j]
              endif else begin
                 Ly1[j] = deltay / 2.
              endelse             
           endif else Ly1[j] = 0
           print,'Ly1[j]=',Ly1[j],' j =',j           
           ; get the right hand side
           if (((wtd_y[j] + deltay / 2.) gt y0) and $
               (wtd_y[j] lt y1)) then begin
              if (y1 - wtd_y[j]) lt (deltay / 2.) then begin                 
                 Ly2[j] = y1 - wtd_y[j]
              endif else if (y0 - wtd_y[j]) lt (deltay / 2.)  then begin
                 Ly2[j] = wtd_y[j] + (deltay / 2.) - y0
              endif else begin
                 Ly2[j] = deltay /2.
              endelse              
           endif else Ly2[j] = 0

        endfor
     endif

     lambda0 = lambdas[spax] - deltalambda / 2.
     lambda1 = lambdas[spax] + deltalambda / 2.
     lambdacenter = where(lambda_grid ge lambda0 and lambda_grid le lambda1, count)
  
     if count ne 0 then begin
        if lambdacenter[0] lt 1 then begin
           lambdastart = 0
        endif else begin
           lambdastart = lambdacenter[0]-1
        endelse
        if lambdacenter[-1]+1 gt n_elements(lambda_grid) then begin
           lambdaend = lambdacenter[-1]
        endif else begin
           lambdaend = lambdacenter[-1] + 1
        endelse  
          
        wtd_lambda = lambda_grid[lambdastart:lambdaend]
        Llambda1 = dblarr(n_elements(wtd_lambda))
        Llambda2 = dblarr(n_elements(wtd_lambda))
        print,'wtd_lambda = ', wtd_lambda
        print, 'lambda_grid=',lambda_grid
        for j=0, n_elements(wtd_lambda)-1 do begin
           
           ; get the left hand side
           if ((((wtd_lambda[j] - deltalambda) / 2.) lt lambda1) and (wtd_lambda[j] gt lambda0)) then begin
              ; first case if point is on the left of the flulambda
              if (wtd_lambda[j] - lambda0) lt deltalambda / 2. then begin
                 Llambda1[j] = wtd_lambda[j] - lambda0
              ; point is on the right of the flulambda
              endif else if (wtd_lambda[j] - lambda1) lt deltalambda / 2. then begin
                 Llambda1[j] = lambda1 - ( wtd_lambda[j] - deltalambda /2.)
              endif else begin
                 Llambda1[j] = deltalambda / 2.
              endelse             
           endif else Llambda1[j] = 0
           
           ; get the right hand side
           if (((wtd_lambda[j] + deltalambda / 2.) gt lambda0) and $
               (wtd_lambda[j] lt lambda1)) then begin
              if (lambda1 - wtd_lambda[j]) lt (deltalambda / 2.) then begin                 
                 Llambda2[j] = lambda1 - wtd_lambda[j]
              endif else begin
                 Llambda2[j] = deltalambda /2.
                 print,'deltalambda/2.=',deltalambda / 2.
              endelse              
           endif else Llambda2[j] = 0

        endfor
     endif

     
