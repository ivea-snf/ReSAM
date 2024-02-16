PRO setColorTable,sw=sw
  COMMON colortable,red,green,blue,blackindex,redindex,greenindex,blueindex,ncolors
  IF n_elements(red) EQ 0 THEN init=1 ELSE init=0
  IF not keyword_set(sw) THEN sw=0 
 
    red=  [255,153, 77,  0,  0,  0, 64,194,255,255,255,255,230,191,153,  0]
    green=[255,204,102, 51,128,153,179,232,255,224,168,105, 38,  0,  0,  0]
    blue= [255,255,255,179,115,  0, 64, 71,  0,  0, 77, 26,  0,  0,  0,  0]

  ncolors=n_elements(red)
  red=rebin(red,2*ncolors)
  red=red(0:29) & red(29:*)=0
  green=rebin(green,2*ncolors)
  green=green(0:29) & green(29:*)=0
  blue=rebin(blue,2*ncolors)
  blue=blue(0:29) & blue(29:*)=0
  red(1)=red(2)
  green(1)=green(2)
  blue(1)=blue(2)

  tvlct,red,green,blue
  blackindex=n_elements(red)-1
  redindex=12*2
  greenindex=5*2
  blueindex=2*2
  ncolors=n_elements(red)
  !P.COLOR=blackindex
END  

COMMON colortable,red,green,blue,blackindex,redindex,greenindex,blueindex,ncolors

SET_PLOT,'ps'
;SET_PLOT,'x'

DEVICE,  xsize=20 , ysize= 18, /color,bits=8, font_size=15.0,yoffset=2,filename='virus.ps'


setcolortable
!p.font=0

;=== declarations ===============================================

c=''








 NS=2001
NZ=Ns
sars = FLTARR(nz)
iav = FLTARR(nz)

time = FLTARR(nz)
times = FLTARR(nz)
sars2 = FLTARR(nz)
iav2 = FLTARR(nz)



close,1
openr, 1, 'fort.22'

close,2
openr, 2, '../25_PNAS_RH90_v1.1/fort.22'

i=0
repeat begin
 READF,1, xx,yy
time(i)=xx
iav(i)=yy




i=i+1
endrep until eof(1)
       close,1

nt1=i-1
i=0
repeat begin


READF,2, xx,yy
sars(i)=yy
times(i)=xx



i=i+1
endrep until eof(2)
       close,1

nt=i-1


print,nt,nt1



plot, [1,1], [1,1], ytitle=' pH ', xtitle='time [s] ',xrange=[0, 3600],yrange=[7,11.5],ystyle=1,xstyle=1,/nodata,charsize=1.5,title='Initial diameter 50 micron'


oplot,time(0:nt1), iav(0:nt1), thick=8, color=blueindex
oplot,times(0:nt), sars(0:nt), thick=8, color=redindex

xx1=200
xx2= 200+180
yy= 10.4
dy=.5
xx3= xx2+50
oplot, [xx1,xx2],[yy,yy],line=0, thick=5, color=blueindex
xyouts,xx3,yy,'RH 50%'
yy=yy+.35
oplot, [xx1,xx2],[yy,yy],line=0, thick=5, color=redindex
xyouts,xx3,yy,'RH 90%'
yy=yy+.35
oplot, [xx1,xx2],[yy,yy],line=0, thick=5, color=greenindex
xyouts,xx3,yy,'RH 99.4%'


close,1
openr, 1, 'fort.22'

close,1
openr, 1, '../25_PNAS_RH99_v1.1/fort.22'

i=0
repeat begin
 READF,1, xx,yy
time(i)=xx
iav(i)=yy




i=i+1
endrep until eof(1)
       close,1

nt1=i-1
i=0
oplot,time(0:nt1), iav(0:nt1), thick=8, color=greenindex


IF ( !D.NAME EQ 'PS') THEN DEVICE,/CLOSE
end
















