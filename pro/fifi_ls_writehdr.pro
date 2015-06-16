function fifi_ls_writehdr, data, header, HISTORY=history, CANCEL=cancel

  cancel = 0

  ; make blank header from data
  mkhdr, hdr, data

  ; check type of input header hash
  t = size(header, /TYPE)
  if t ne 11 || obj_class(header) ne 'HASH' then begin
     print, 'Input header variable must be a hash.'
     cancel = 1
     return, hdr
  endif

  ; get comments for keys
  if header.haskey('COMMENTS') then begin
     comment = header['COMMENTS']
  endif else comment = hash()

  ; get the keys from the hash and sort them
  hdr_keys = header.keys()
  idx = sort(hdr_keys.toarray())
  hdr_keys = hdr_keys[idx]

  ; mark a few keys that should not be added to
  ; the header
  skip_keys = ['COMMENTS', 'DATAPATH', 'PKGPATH']
  skip = bytarr(n_elements(hdr_keys))
  foreach k,skip_keys do begin
     skip_idx = where(hdr_keys eq k, count)
     if count gt 0 then skip[skip_idx] = 1B
  endforeach

  ; and a few filename keywords to take
  ; the base name from
  filename_keys = ['BPM','FILENAME','LINFILE']
  file = bytarr(n_elements(hdr_keys))
  foreach k,filename_keys do begin
     file_idx = where(hdr_keys eq k, count)
     if count gt 0 then file[file_idx] = 1B
  endforeach

  ; add all keys from header structure
  foreach key,hdr_keys,iter do begin
     if skip[iter] then continue
     if file[iter] then begin
        val = file_basename(header[key]) 
     endif else val = header[key]
     if comment.haskey(strupcase(key)) then begin
        com = comment[key] 
     endif else com = ''
     sxaddpar, hdr, key, val, com
  endforeach

  if n_elements(history) gt 0 then begin
     sxaddhist, history, hdr
  endif

  return,hdr
end
