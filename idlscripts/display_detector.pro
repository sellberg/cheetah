pro display_detector, xarray, yarray, raw, wid, fout=fout

if n_elements(wid) eq 0 then wid = 0

pixsize = 109.92e-6  ;1.
dsize = 1024

sx = size(xarray)
sy = size(yarray)

xmax = max(xarray/pixsize)
xmin = min(xarray/pixsize)
ymax = max(yarray/pixsize)
ymin = min(yarray/pixsize)

xarraypix = xarray/pixsize - xmin
yarraypix = yarray/pixsize - ymin

nx = long(xmax - xmin) + 2
ny = long(ymax - ymin) + 2

detector = fltarr(nx,ny)
for i =0,sx[1]-1 do begin
for j =0,sx[2]-1 do begin
 xcoord = round(xarraypix[i,j]) 
 ycoord = round(yarraypix[i,j]) 
 detector[xcoord,ycoord] = raw[i,j] ;1.
endfor
endfor

display, congrid(detector,dsize,dsize)^0.3, wid=wid
;display, detector^0.3, wid=wid
if n_elements(fout) ne 0 then write_png, fout, bytscl(detector^0.3)

end