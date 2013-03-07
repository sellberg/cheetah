
;; the metrology data measured the corners of the silicon sensor. This is
;; slightly larger than the asic. So the measurements have to be corrected
;; to get the correct positions for the corners of the asic

function asic_position_from_sensor_position, input

input = reform(input)
output = input*0.
pix = 109.92

temp = shift(input,0,-1)
len = sqrt((temp[0,*] - input[0,*])^2 + (temp[1,*] - input[1,*])^2)
xlen = (temp[0,*] - input[0,*])
ylen = (temp[1,*] - input[1,*])

il = where(len lt 30000,complement=ilc)

lshift = len*0.

lshift[il] = (len[il] - 185*pix)/2
lshift[ilc] = (len[ilc] - 391*pix)/2

yshift = ylen*lshift/len
xshift = xlen*lshift/len

tempx = shift(xshift,1)
tempy = shift(yshift,1)

output[0,0:3] = input[0,0:3] + xshift[0:3] - tempx[0:3] 
output[1,0:3] = input[1,0:3] + yshift[0:3] - tempy[0:3]

;output[0,2:3] = input[0,2:3] - xshift[2:3] + tempx[2:3] 
;output[1,2:3] = input[1,2:3] - yshift[2:3] + tempy[2:3]

return, output
end

;; Main program

pro make_pixel_map_from_LCLS_metrology_data


datapath = "D:\work\water LCLS 2011\new metrology\"
flist = file_search(datapath+"CSPAD2-Alignment-PostRun3_Q*.csv")
slist = size(flist)

print, flist
print, slist

pixsize = 109.92
qdim = fltarr(4,4)

;; read in the metrology data from .csv files and store them in an array
for i =0,slist[1]-1 do begin
  qlist = read_ascii(flist[i],delimiter=",")
  tag = tag_names(qlist)
  sqlist = size(qlist.(tag[0]))
  if i eq 0 then qlistsave = fltarr(4,3,sqlist[2])

  test = fltarr(32)
  ;; Give each quadrant an appropriate rotation
  if i eq 0 then begin
    test = qlist.(tag[0])[1,0:31]
    qlist.(tag[0])[1,0:31] = -qlist.(tag[0])[2,0:31]
    qlist.(tag[0])[2,0:31] = test[*]
    qlist.(tag[0])[1,0:31] += -min(qlist.(tag[0])[1,0:31])
    qlist.(tag[0])[2,0:31] += -min(qlist.(tag[0])[2,0:31])    
  endif
    
  if i eq 2 then begin
    test = qlist.(tag[0])[1,0:31]
    qlist.(tag[0])[1,0:31] = qlist.(tag[0])[2,0:31]
    qlist.(tag[0])[2,0:31] = -test[*]
    qlist.(tag[0])[1,0:31] += -min(qlist.(tag[0])[1,0:31])
    qlist.(tag[0])[2,0:31] += -min(qlist.(tag[0])[2,0:31])    
  endif
  
  if i eq 3 then begin
    qlist.(tag[0])[1,0:31] = -qlist.(tag[0])[1,0:31] 
    qlist.(tag[0])[2,0:31] = -qlist.(tag[0])[2,0:31]
    qlist.(tag[0])[1,0:31] += -min(qlist.(tag[0])[1,0:31])
    qlist.(tag[0])[2,0:31] += -min(qlist.(tag[0])[2,0:31])    
  endif
  
  qlistsave[i,0,*] = qlist.(tag[0])[1,0:31]
  qlistsave[i,1,*] = qlist.(tag[0])[2,0:31] 
  qlistsave[i,2,*] = qlist.(tag[0])[3,0:31] 
  
  qdim[i,0] = min(qlist.(tag[0])[1,0:31],/NAN)
  qdim[i,1] = max(qlist.(tag[0])[1,0:31],/NAN)
  qdim[i,2] = min(qlist.(tag[0])[2,0:31],/NAN)
  qdim[i,3] = max(qlist.(tag[0])[2,0:31],/NAN)
endfor

xmax = max(qdim[*,1],/NAN)
xmin = min(qdim[*,0],/NAN)
ymax = max(qdim[*,3],/NAN)
ymin = min(qdim[*,2],/NAN)
print, xmin, xmax, ymin, ymax

nx = (xmax-xmin)/pixsize
ny = (ymax-ymin)/pixsize
  
offset = 850  
image = fltarr(nx+offset+300,ny+offset+300)

xquad = [0-16   ,offset-20 ,offset+12  ,+25] +100
yquad = [offset-15 ,offset+22 ,0+20       ,-10] +100

for i =0,3 do begin ;;slist[1]-1 do begin
  xshift = xquad[i]
  yshift = yquad[i]
  for j=0,7 do begin
    m = 4*j
    n = 4*(j+1) -1
     
    xmax2 = (max(qlistsave[i,0,m:n],/NAN) - xmin - 2*pixsize)/pixsize  +xshift
    xmin2 = ((min(qlistsave[i,0,m:n],/NAN) - xmin)/pixsize)            +xshift
    ymax2 = ((max(qlistsave[i,1,m:n],/NAN) - ymin - 2*pixsize)/pixsize)+yshift
    ymin2 = ((min(qlistsave[i,1,m:n],/NAN) - ymin)/pixsize)            +yshift
    if finite(xmax2) and finite(xmin2) and finite(ymax2) and finite(ymin2) then begin 
        image[xmin2:xmax2,ymin2:ymax2] = 1.
        output = asic_position_from_sensor_position(qlistsave[i,0:1,m:n])
        for iii=0,3 do print, i, j, " 0", qlistsave[i,0,iii], output[0,iii]
        for iii=0,3 do print, i, j, " 1", qlistsave[i,1,iii], output[1,iii]
        qlistsave[i,0:1,m:n] = output[0:1,0:3]
    endif
  endfor  
endfor  
display, congrid(image,512,512), wid=1

;; fix up the dimensions
qlistsave[*,0,*] +=  - (max(qlistsave[*,0,*]) + min(qlistsave[*,0,*]))/2
qlistsave[*,1,*] +=  - (max(qlistsave[*,1,*]) + min(qlistsave[*,1,*]))/2

;; add in quadrant positions from Rick's geometry
quadshift = fltarr(2,4)
quadshift[0,0]  = -8.80731
quadshift[0,1]  = -4.71673
quadshift[0,2]  = -5.07504
quadshift[0,3]  = -13.0048
quadshift[1,0]  = -6.89487
quadshift[1,1]  = -5.22273
quadshift[1,2]  = 2.76495
quadshift[1,3]  = 5.84621
for i=0,3 do begin
for j=0,1 do begin
  qlistsave[i,j,*] += quadshift[j,i]
endfor
endfor


;stop

;; ______________________________________________
;; put these coordinates into a pixel mask

datapath2 = "D:\Projects CFEL\cpad_geometry\"
geometry_file = "cspad_pixelmap_with_split.h5"
filename = datapath2+geometry_file
xarray2 = read_h5(filename,field='x') /pixsize   ;; 1e-6 /pixsize
yarray2 = read_h5(filename,field='y') /pixsize   ;; 1e-6 /pixsize
xarray = xarray2 ;-yarray2
yarray = yarray2 ;xarray2


rows = 194
cols = 185

zarray = fltarr(8*rows,8*cols)

qorder = [0,1,2,3]
;;order = [0,1,2,3,5,4,6,7] ;;  [1,0,3,2,4,5,7,6] ;; 
order = lonarr(4,8)
; for the rotation case
order[0,*] = [1,0,3,2,4,5,7,6]
order[1,*] = [1,0,3,2,4,5,7,6]
order[2,*] = [1,0,3,2,4,5,7,6]
order[3,*] = [1,0,3,2,4,5,7,6]


qtemp = fltarr(2,sqlist[2]-1)
for i = 0,3 do begin
  qtemp[*,*] = shift(qlistsave[i,0:1,*],0,0,(qorder[i]-1)*8)
  xshift = xquad[i]
  yshift = yquad[i]
for j = 0,7 do begin
  m = 4*j
  n = 4*(j+1) -1
   xmax2 = (max(qlistsave[i,0,m:n],/NAN) - xmin - 2*pixsize)/pixsize  +xshift
   xmin2 = ((min(qlistsave[i,0,m:n],/NAN) - xmin)/pixsize)            +xshift
   ymax2 = ((max(qlistsave[i,1,m:n],/NAN) - ymin - 2*pixsize)/pixsize)+yshift
   ymin2 = ((min(qlistsave[i,1,m:n],/NAN) - ymin)/pixsize)            +yshift
  
   xmin_a1 = qorder[i]*(2*rows)
   xmax_a1 = xmin_a1 + rows -1
   ymin_a1 = order[qorder[i],j]*cols
   ymax_a1 = ymin_a1 + cols - 1
   
   xmin_a2 = qorder[i]*(2*rows) + rows
   xmax_a2 = xmin_a2 + rows -1
   ymin_a2 = ymin_a1
   ymax_a2 = ymax_a1
   
   
   xmin4 = min(xarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1])
   ymin4 = min(yarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1])
   
   xmax4 = max(xarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1])
   ymax4 = max(yarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1])
   
   which_asic = 0
   if max(xarray[xmin_a2:xmax_a2,ymin_a2:ymax_a2]) gt max(xarray[xmin_a1:xmax_a1,ymin_a1:ymax_a1]) then which_asic=1
   
  if finite(xmax2) and finite(xmin2) and finite(ymax2) and finite(ymin2) then begin
   ; if which_asic eq 1 then begin
      xarray[xmin_a1:xmax_a1,ymin_a1:ymax_a1] += -xmin4 + xmin2 -(xmax-xmin)/(2*pixsize) ;;-ymin2 + ymin2 +(ymax-ymin)/(2*pixsize) ;;
      yarray[xmin_a1:xmax_a1,ymin_a1:ymax_a1] += -ymin4 + ymin2 -(ymax-ymin)/(2*pixsize) ;;+xmin2 -(xmax-xmin)/(2*pixsize) ;;
    
    ;  xarray[xmin_a2:xmax_a2,ymin_a2:ymax_a2] += -xmax4 + xmax2 -(xmax-xmin)/(2*pixsize) ;;-ymin2 + ymin2 +(ymax-ymin)/(2*pixsize) ;;
    ;  yarray[xmin_a2:xmax_a2,ymin_a2:ymax_a2] += -ymax4 + ymax2 -(ymax-ymin)/(2*pixsize) ;;+xmin2 -(xmax-xmin)/(2*pixsize) ;;
    ;endif else begin
    ;  xarray[xmin_a1:xmax_a1,ymin_a1:ymax_a1] += -xmax4 + xmax2 -(xmax-xmin)/(2*pixsize)
    ;  yarray[xmin_a1:xmax_a1,ymin_a1:ymax_a1] += -ymax4 + ymax2 -(ymax-ymin)/(2*pixsize) 
      xarray[xmin_a2:xmax_a2,ymin_a2:ymax_a2] += -xmin4 + xmin2 -(xmax-xmin)/(2*pixsize) ;;-ymin2 + ymin2 +(ymax-ymin)/(2*pixsize) ;;
      yarray[xmin_a2:xmax_a2,ymin_a2:ymax_a2] += -ymin4 + ymin2 -(ymax-ymin)/(2*pixsize) ;;+xmin2 -(xmax-xmin)/(2*pixsize) ;; 
    ;endelse
     
     
    ;; add in some rotation
    angle_av = 0.
    for ith=0,3 do begin
      ith2 = (ith+1) mod 4 
      sign = (ith mod 2)*2 -1
      if abs(qlistsave[i,0,m+ith2]-qlistsave[i,0,m+ith]) lt abs(qlistsave[i,1,m+ith2]-qlistsave[i,1,m+ith]) then begin
        flag = 1
        flag2 = 0
      endif else begin
        flag = 0
        flag2 = 1
      endelse
    
      angle = -sign * atan(  (qlistsave[i,flag2,m+ith2]-qlistsave[i,flag2,m+ith]) / (qlistsave[i,flag,m+ith2]-qlistsave[i,flag,m+ith])   )
      angle_av += angle
      print, angle, angle_av/4.
    endfor
   angle_av *= 1./4.
   temp = xarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1]
   xarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1] = cos(angle_av)*temp - sin(angle_av)*yarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1]
   yarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1] = sin(angle_av)*temp + cos(angle_av)*yarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1]   
 
   ;; fill out z values
   zarray[xmin_a1:xmax_a2,ymin_a1:ymax_a1] = total(qlistsave[i,2,m:n])/4.
  
  endif ;check for finite values  
endfor
endfor

;; fix some missing panels
i = 3
j = 4
m = 4*j
n = 4*(j+1) -1
j2 = 5
m2 = 4*j2
n2 = 4*(j2+1) -1
xmin3 = qorder[i]*(2*rows)
xmax3 = xmin3 + 2*rows -1
ymin3 = order[qorder[i],j]*cols
ymax3 = ymin3 + cols - 1
xmin4 = qorder[i]*(2*rows)
xmax4 = xmin4 + 2*rows -1
ymin4 = order[qorder[i],j2]*cols
ymax4 = ymin4 + cols - 1
yarray[xmin3:xmax3,ymin3:ymax3] = yarray[xmin4:xmax4,ymin4:ymax4] + cols + 5
xarray[xmin3:xmax3,ymin3:ymax3] = xarray[xmin4:xmax4,ymin4:ymax4] 

i = 0
j = 5
m = 4*j
n = 4*(j+1) -1
j2 = 4
m2 = 4*j2
n2 = 4*(j2+1) -1
xmin3 = qorder[i]*(2*rows)
xmax3 = xmin3 + 2*rows -1
ymin3 = order[qorder[i],j]*cols
ymax3 = ymin3 + cols - 1
xmin4 = qorder[i]*(2*rows)
xmax4 = xmin4 + 2*rows -1
ymin4 = order[qorder[i],j2]*cols
ymax4 = ymin4 + cols - 1
xarray[xmin3:xmax3,ymin3:ymax3] = xarray[xmin4:xmax4,ymin4:ymax4] - (cols + 5)
yarray[xmin3:xmax3,ymin3:ymax3] = yarray[xmin4:xmax4,ymin4:ymax4] 

i = 3
j = 7
m = 4*j
n = 4*(j+1) -1
j2 = 6
m2 = 4*j2
n2 = 4*(j2+1) -1
xmin3 = qorder[i]*(2*rows)
xmax3 = xmin3 + 2*rows -1
ymin3 = order[qorder[i],j]*cols
ymax3 = ymin3 + cols - 1
xmin4 = qorder[i]*(2*rows)
xmax4 = xmin4 + 2*rows -1
ymin4 = order[qorder[i],j2]*cols
ymax4 = ymin4 + cols - 1
xarray[xmin3:xmax3,ymin3:ymax3] = xarray[xmin4:xmax4,ymin4:ymax4] - (cols + 5)
yarray[xmin3:xmax3,ymin3:ymax3] = yarray[xmin4:xmax4,ymin4:ymax4] 

;; fix up the dimensions
xarray +=  - (max(xarray) + min(xarray))/2
yarray +=  - (max(yarray) + min(yarray))/2

;; fix up the dimensions
xarray *= pixsize/1.e6
yarray *= pixsize/1.e6

;filename3 = "D:\work\water LCLS 2011\cpad_geometry\r0589-RawSum.h5"
filename3 = "F:\Work back up\LCLSWater2011\cspad_water\r0005-RawSum.h5"
raw = read_h5(filename3)
dsize = 512
display, congrid(raw,dsize,dsize)^0.3, wid=3

display_detector, xarray, yarray, raw, 0 

outname = datapath+"CSPAD2-Alignment-PostRun3_pixelmap.h5"

;;
  ;; Save it
  ;;
  print,'Writing data to file: ',outname
  fid = H5F_CREATE(outname) 

  datatype_id = H5T_IDL_CREATE(xarray) 
  dataspace_id = H5S_CREATE_SIMPLE(size(xarray,/DIMENSIONS)) 
  dataset_id = H5D_CREATE(fid,'x',datatype_id,dataspace_id) 
  H5D_WRITE,dataset_id,xarray
  H5D_CLOSE,dataset_id   
  H5S_CLOSE,dataspace_id 
  H5T_CLOSE,datatype_id 

  datatype_id = H5T_IDL_CREATE(yarray) 
  dataspace_id = H5S_CREATE_SIMPLE(size(yarray,/DIMENSIONS)) 
  dataset_id = H5D_CREATE(fid,'y',datatype_id,dataspace_id) 
  H5D_WRITE,dataset_id,yarray
  H5D_CLOSE,dataset_id   
  H5S_CLOSE,dataspace_id 
  H5T_CLOSE,datatype_id 

  datatype_id = H5T_IDL_CREATE(zarray) 
  dataspace_id = H5S_CREATE_SIMPLE(size(zarray,/DIMENSIONS)) 
  dataset_id = H5D_CREATE(fid,'z',datatype_id,dataspace_id) 
  H5D_WRITE,dataset_id,zarray 
  H5D_CLOSE,dataset_id   
  H5S_CLOSE,dataspace_id 
  H5T_CLOSE,datatype_id 
  
  H5F_CLOSE,fid 



end