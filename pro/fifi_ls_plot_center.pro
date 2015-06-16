
pro fifi_ls_plot_center
  fn = file_search('*lambda_xy*')
  for i=0, n_elements(fn)-1 do begin
     d1=mrdfits(fn[i])
     plot,d1[2,2,*,0],d1[2,2,*,1],psym=4
     print,'filename=',fn[i]
     wait,1
     
  endfor
end


pro all_frames,q
  for i=0, 191 do begin
     ximgtool,q[*,*,i]
     wait,0.07
     print,i
  endfor
end

  
