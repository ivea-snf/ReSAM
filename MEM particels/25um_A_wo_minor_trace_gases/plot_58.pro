
PRO setColorTable,sw=sw
  COMMON colortable,red,green,blue,blackindex,redindex,greenindex,blueindex,ncolors
ncolors=57
n256=256

red=intarr(n256)
green=intarr(n256)
blue=intarr(n256)

close,1
openr,1,'/home/bluo/IVEA/Reid/color.dat'
for i  =  ncolors-1,0,-1 do begin
readf,1,rr,gg,bb
red(i)= rr
green(i) = gg
blue(i)= bb
endfor


blackindex= 251
red(blackindex) =0
blue(blackindex) =0
green(blackindex) =0
redindex= 252
red(redindex) =255
green(redindex) =0
blue(redindex) =0

blueindex= 253
red(blueindex) =0
green(blueindex) =0
blue(blueindex) =255

greenindex= 254
red(greenindex) =130
green(greenindex) =130
blue(greenindex) =130

tvlct,red,green,blue

!P.COLOR=blackindex
END  

COMMON colortable,red,green,blue,blackindex,redindex,greenindex,blueindex,ncolors

SET_PLOT,'ps'
;SET_PLOT,'x'

DEVICE,  xsize=20 , ysize= 18, /color,bits=8, font_size=15.0,yoffset=2,filename='resam.ps'

setcolortable
!p.font=0

;=== declarations ===============================================

c=''







;  ==== plotting the temperature
close,1

file='output_shells_M.dat'

 NT=20001
NS=200
close,1
file='output_shells_M.dat'
openr,1,file
npsolid=5
xsat=fltarr(npsolid)


readf,1,xx1,xx2,xx3,xx4,np
close,1
file='output_shells_M.dat'


color=fltarr(2*Ns)

x = FLTARR(NT,ns)
xsolid = FLTARR(NT,ns)

y= FLTARR(2*ns)

time = FLTARR(NT)
rad = FLTARR(NT)

aw =fltarr(NT,nS)
m =fltarr(NT,nS,NP+1+npsolid)

gammaH =fltarr(NT,nS)
titerA =fltarr(NT,nS)

x18 = FLTARR(NP)
ymin18 = FLTARR(NP+1+npsolid)
ymax18 = FLTARR(NP+1+npsolid)

ns2 = intARR(Nt)


NAme=strarr(NP+6)
      name=['water activity','molality H2O', 'molality organics small', 'molality total acetic acid', 'molality total ammonium', 'molality total CO!D2!N', 'pH', 'molality OH!U-!N', 'molality CH!D3!NCOO!U-!N', 'molality CH!D3!NCOOH', 'molality NH!D4!NCH!D3!NCOO', 'molality organics big', 'molality NH!D4!N!U+!N', 'molality CO!D2!N', 'molality CO!D3!N!U-2!N', 'molality HCO!D3!N!U-!N', 'molality Na!U+!N', 'molality Cl!U-!N', 'molality NO!D3!N!U-!N','X!U+!N','Y!U-!N','log!D10!N(titer Flu) ','log!D10!N(titer SARS) ','molality NH!D3!N', 'H3PO4','H2PO4-','HPO4-2', 'HSO4-','SO4-2', 'NaCl sauration ratio']

if NP eq 29 then       name=['water activity','molality H2O', 'molality organics small', 'molality total acetic acid', 'molality total ammonium', 'molality total CO!D2!N', 'pH', 'molality OH!U-!N', 'molality CH!D3!NCOO!U-!N', 'molality CH!D3!NCOOH', 'molality NH!D4!NCH!D3!NCOO', 'molality organics big', 'molality NH!D4!N!U+!N', 'molality CO!D2!N', 'molality CO!D3!N!U-2!N', 'molality HCO!D3!N!U-!N', 'molality Na!U+!N', 'molality Cl!U-!N', 'molality NO!D3!N!U-!N','X!U+!N','Y!U-!N','log!D10!N(titer Flu) ','log!D10!N(titer SARS) ','molality NH!D3!N', 'H3PO4','H2PO4-','HPO4-2', 'HSO4-','SO4-2', 'm-NaCl crystal', 'NaCl sauration ratio']


if np eq 28 then      name=['water activity','molality H2O', 'molality organics small', 'molality total acetic acid', 'molality total ammonium', 'molality total CO!D2!N', 'pH', 'molality OH!U-!N', 'molality CH!D3!NCOO!U-!N', 'molality CH!D3!NCOOH', 'molality NH!D4!NCH!D3!NCOO', 'molality organics big', 'molality NH!D4!N!U+!N', 'molality CO!D2!N', 'molality CO!D3!N!U-2!N', 'molality HCO!D3!N!U-!N', 'molality Na!U+!N', 'molality Cl!U-!N', 'molality NO!D3!N!U-!N','X!U+!N','Y!U-!N','log!D10!N(titer Flu) ','log!D10!N(titer SARS) ','molality NH!D3!N', 'H3PO4','H2PO4-','HPO4-2', 'HSO4-','SO4-2', 'NaCl sauration ratio']


if NP eq 32 then       name=['water activity','molality H2O', 'molality organics small', 'molality total acetic acid', 'molality total ammonium', 'molality total CO!D2!N', 'pH', 'molality OH!U-!N', 'molality CH!D3!NCOO!U-!N', 'molality CH!D3!NCOOH', 'molality NH!D4!NCH!D3!NCOO', 'molality organics big', 'molality NH!D4!N!U+!N', 'molality CO!D2!N', 'molality CO!D3!N!U-2!N', 'molality HCO!D3!N!U-!N', 'molality Na!U+!N', 'molality Cl!U-!N', 'molality NO!D3!N!U-!N','X!U+!N','Y!U-!N','log!D10!N(titer Flu) ','log!D10!N(titer SARS) ','molality NH!D3!N', 'H3PO4','H2PO4-','HPO4-2', 'HSO4-','SO4-2', 'Oxalic acid total', 'OA-', 'OA-2','m-NaCl crystal', 'NaCl sauration ratio']


if NP eq 33 then       name=['water activity','molality H2O', 'molality organics small', 'molality total acetic acid', 'molality total ammonium', 'molality total CO!D2!N', 'pH', 'molality OH!U-!N', 'molality CH!D3!NCOO!U-!N', 'molality CH!D3!NCOOH', 'molality NH!D4!NCH!D3!NCOO', 'molality organics big', 'molality NH!D4!N!U+!N', 'molality CO!D2!N', 'molality CO!D3!N!U-2!N', 'molality HCO!D3!N!U-!N', 'molality Na!U+!N', 'molality Cl!U-!N', 'molality NO!D3!N!U-!N','X!U+!N','Y!U-!N','log!D10!N(titer Flu) ','log!D10!N(titer SARS) ','molality NH!D3!N', 'H3PO4','H2PO4-','HPO4-2', 'HSO4-','SO4-2', 'Oxalic acid total ', 'OA-', 'OA-2','NH42-OA fraction','NaCl fraction', 'NaCl sauration ratio','NH42-OA saturation ratio']

if NP eq 36 then       name=['water activity','molality H2O', 'molality organics small', 'molality total acetic acid', 'molality total ammonium', 'molality total CO!D2!N', 'pH', 'molality OH!U-!N', 'molality CH!D3!NCOO!U-!N', 'molality CH!D3!NCOOH', 'molality NH!D4!NCH!D3!NCOO', 'molality organics big', 'molality NH!D4!N!U+!N', 'molality CO!D2!N', 'molality CO!D3!N!U-2!N', 'molality HCO!D3!N!U-!N', 'molality Na!U+!N', 'molality Cl!U-!N', 'molality NO!D3!N!U-!N','X!U+!N','Y!U-!N','log!D10!N(titer Flu) ','log!D10!N(titer SARS) ','molality NH!D3!N', 'H3PO4','H2PO4-','HPO4-2', 'HSO4-','SO4-2', 'Oxalic acid total ', 'OA-', 'OA-2','NH!D4!NHC!D2!NO!D4!N fraction','(NH!D4!N)!D2!NC!D2!NO!D4!N fraction','NaHC!D2!NO!D4!N fraction','Na!D2!NC!D2!NO!D4!N fraction','NaCl fraction', 'NH!D4!NHC!D2!NO!D4!N Sauration Ratio','(NH4)!D2!NC!D2!NO!D4!N Saturation Ratio', 'NaHC!D2!NO!D4!N Saturation Ratio' ,'Na!D2!NC!D2!NO!D4!N Saturation Ratio','NaCl Saturation Ratio']

print, 'np  = ', np,n_elements(name)


klog=intarr(NP+1+npsolid)

for k = 0, np do begin
klog(k)=1
ymin18(k)=1D-8
ymax18(k)=1D1
endfor

klog(0)=0
klog(1)=0
klog(6)=0
klog(3)=1
klog(4)=1

klog(Np+1)=1   ; Saturation

;ymax18(Np+1)=100

Ymin18(0)=0
Ymax18(0)=1

Ymin18(1)=0
Ymax18(1)=60

Ymin18(3)=1D-7

Ymin18(4)=1D-7


Ymin18(6)=2.5
Ymax18(6)=7.5

Ymin18(7)=1D-14
Ymax18(7)=1d-1


xmin=1D-5
xmax= 1D6
if k ge np then begin
klog(k)=1
ymin18(k)=.1
;ymin18(Np+1)=.1
endif

zmin= 0
zmax=.03

; READF,1, timee, radd, nss,x00

close,1
openr,1,'output_shells_M.dat'


i=0

repeat begin

 READF,1, timee, radd, nss,x00


  time(i)=timee
  rad(i)=radd

if (radd*1D4*1.1 gt zmax and timee le xmax) then zmax= radd*1D4*1.1 
  ns2(i)=nss


x(i,0)=x00*1D4


for j = 0, nss-1 do begin
;if np le 32 then readf,1,ii, xx,aww, x18,smisch, gammahh,phh,xxs
if np eq 29  then  readf,1,ii, xx,aww, x18,smisch, gammahh,phh,xxs, sAO
if np eq 36  then  readf,1,ii, xx,aww, x18, gammahh,phh,xxs, xsat




;              write(3,'(I6,501E14.5)') I,
;     &        x(I+1),aw,(M(j),J=1,NP),smisch, gamahh,pHH



x(i,j+1)= xx*1D4 ; um
xsolid(i,j)= xxs*1D4 ; um



aw(i,j)=aww
m(i,j,0) = aww

for k =0,np-1 do begin
m(i,j,k+1) = x18(k)
endfor

gammaH(i,j)=gammahh

m(i,j,6)=10d0^(-phh)

;m(i,j,6)=gammahh* x18(5)



m(i,j,29)=m(i,j,29)+m(i,j,30)+m(i,j,31)

xsolidd=0
for  kk = np-npsolid, np-1 do begin
xsolidd=xsolidd +x18(kk)
m(i,j,kk+1) = 0d0
endfor

if xsolidd gt 1D-30 then begin
for  kk = np-npsolid, np-1 do begin
m(i,j,kk+1)=x18(kk)/xsolidd
endfor
endif


for kk=0,npsolid-1 do m(i,j,np+kk+1)=xsat(kk)




xiav=m(i,j,21)
if xiav ge 1D-7 then m(i,j,21) = alog(xiav)/alog(10d0)
if xiav le 1D-7 then m(i,j,21) = -7
sarsa=m(i,j,22)
if sarsa ge 1D-7 then m(i,j,22) = alog(sarsa)/alog(10d0)
if sarsa le 1D-7 then  m(i,j,22) = -7
endfor



i=i+1

endrep until eof(1)
       close,1

nt=i-1

if (timee gt xmax) then xmax=timee

print, xmin,xmax


;plot aw

k=0
irepeat =0
ii=1
if np eq 33 then ii =2





for k=0,np+npsolid  do begin

klog(np+2)=1


dmin=1000
dmax=0

for i=2,nt-1 do begin
nz=ns2(i)

for j=0,nz-1 do begin
if m(i,j,k) lt dmin then dmin=m(i,j,k)
if m(i,j,k) gt dmax then dmax=m(i,j,k)
endfor
endfor

if k gt np then begin
ymin=.01
klog(k)=1

endif


if dmin eq 0 then begin
; dmin=dmax*1D-10
klog(k)=1
endif

klog(3)=1
klog(4)=1
klog(2)=1


ymin=dmin


if klog(k) eq 1  then ymax=dmax*2
if klog(k) eq 0  then ymax=dmax*1.1

if klog(k) eq 1 then ymin=1D-7
if k eq 7 then ymin=1D-12

if k eq 3 then ymin=1D-6

if k eq 2 then ymin=1D-6
if k eq 4 then ymin=1D-6
if k eq 16 or k eq 17   then begin
klog(k)=0
endif


;!if k eq np+2 then ymax=100


if k eq 1 then begin
ymin=ymin18(k)
ymax=ymax18(k)
endif

if k eq 6 then begin
ymin=ymin18(k)
ymax=ymax18(k)
endif

if k eq 3 then print, 'ymin ,ymax,', ymin, ymax
if ymax le 1D-20 then ymin=1D-20
if ymax le 1D-20 then ymax=1D-10


if klog(k) eq 1 then ymin =alog(ymin)
if klog(k) eq 1 then ymax =alog(ymax)
 if k eq 0 then ymin=ymin18(k)
 if k eq 0 then ymax=ymax18(k)
 if k eq 6 then ymin=ymin18(k)
 if k eq 6 then ymax=ymax18(k)
 if k eq 6 then klog(k)=0


;if k eq np+1 then ymax=alog(20)
if k eq np+1 then klog(k)=1
if k ge np+1 then ymin=alog(.1)

if k ge np+1 and ymax le alog(10.) then ymax=alog(10.)

if k ge np-npsolid+1 and k le  np then begin
klog(k)= 1
ymin=alog(1D-6)
ymax=alog(2.)
endif

 
if k ge 21 and k le 22 then begin
ymin=-5
ymax=1
klog(k)=0

print, 'titer',  name(k)
endif

repeatpH:

!P.POSITION=[0.15,0.15,0.76,0.95]

!y.style    = 1
!x.style    = 1
!y.title    = 'radius [!9m!Xm]'
!x.title    = 'time [s]'

close,1


zmin=0

 plot_oi, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[zmin,zmax],/nodata,charsize=1.5

  

for i=0,nt-1 do begin

nz=ns2(i)

NZ1=0

ytime= fltarr(4*ns)

l=0
for j= 0 , nz-1 do begin

ytime(l)=time(i)
y(l) = x(i,j)
xxx=xsolid(i,j)



if xxx gt x(i,j)+1D-10   then begin
color(l) = 100


if k ge np-npsolid+1 and k le  np and  m(i,j, k) ge 1D-30 THEN color(l) = ncolors  *(alog(m(i,j, k)) -ymin) /(ymax-ymin)


l=l+1
ytime(l)=time(i)
color(l) = 100

if k ge np-npsolid+1 and k le  np and  m(i,j, k) ge 1D-30 THEN color(l) = ncolors  *(alog(m(i,j, k)) -ymin) /(ymax-ymin)

if j eq 0 and k ge np-npsolid+1 and k le  np THEN print,i,j, color(l), m(i,j, k)

y(l) = xxx


l=l+1
ytime(l)=time(i)
y(l) = xxx
endif


c1=m(i,j,k)
if klog(k) eq 1 then c1 =alog(c1)
if k eq 6 then begin
c1=m(i,j,k)
c1= - alog(c1)/alog(10d0)




endif
xjj=(ncolors)*(c1-ymin)/(ymax-ymin)
if k ge np-npsolid+1 and k le  np THEN XJJ=0

if xjj le 0 then xjj =0
color(l) = xjj
l=l+1
ytime(l)=time(i)
color(l) = xjj
y(l) = x(i,j+1)



l=l+1


endfor



nz1=0
l=l-1

if ytime(nz1) ge xmin and ytime(nz1) le xmax  then plots, ytime(0:l), y(0:l), thick=5, color=color(0:l)

if ytime(nz1) ge xmin and ytime(nz1) le xmax  then plots, [ytime(nz1), ytime(nz1)], [0d0, x(i,0)], thick=nthick, color=greenindex


endfor

if  k eq  np then begin

endif



i=0

	 plot_oi, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[zmin,zmax],/nodata,charsize=1.5,/noerase


xjs = FLTARR(2,ncolors)
xjcolors=FLTARR(2,ncolors)
xs=findgen(2)


ys=findgen(ncolors)/(ncolors*1d0)*(ymax-ymin)+ymin
xjcolors(0,*)= ((ys-ymin)/(ymax-ymin))*(ncolors)
xjcolors(1,*)= ((ys-ymin)/(ymax-ymin))*(ncolors)
ys=ys


 
!x.style    = 4
!y.style    = 1
!y.ticks    = 0
!x.ticks    = 1
!x.title    = ' '
!y.title    = name(k)






;=== ploting scale =====================================================


!P.POSITION=[.9,.2,.95,.95]

ymin1=(ymin)
ymax1=(ymax)

if k eq 3 then print, 'ff', ymin,ymax

if klog(k) eq 1 then ymin1 =exp(ymin)
if klog(k) eq 1 then ymax1 =exp(ymax)

print, name(k), ymin1,ymax1
if klog(k) eq 0 then contour,/noerase,xjs,xs,ys,/nodata,yrange=[ymin1,ymax1],charsize=1.25

if klog(k) eq 1 then contour,/noerase,xjs,xs,ys,/nodata,charsize=1.25,yrange=[ymin1,ymax1],/ylog

TV, xjcolors,$                                
    !x.window(0), !y.window(0),$
    xsize = !x.window(1) - !x.window(0),$
    ysize = !y.window(1) - !y.window(0),/norm

;=============
if klog(k) eq 0 then cONtoUR,/noerase,xjs,xs,ys,/nodata,charsize=1.25,yrange=[ymin1,ymax1]
if klog(k) eq 1 then cONtoUR,/noerase,xjs,xs,ys,/nodata,charsize=1.25,yrange=[ymin1,ymax1],/ylog


!y.style    = 1
!x.style    = 1
!y.title    = 'radius [um]'
!x.title    = 'time [s]'
!y.ticks    = 0
!x.ticks    = 0
if irepeat le 2 and k eq 6 then begin
irepeat=irepeat+1
ymin=1
ymax=6
if irepeat eq 2 then begin
name(6) = 'pH'
ymax=11
ymin=1

for i=0,nt-1 do begin

nz=ns2(i)
for j= 0 , nz-1 do begin
m(i,j,k) = m(i,j,k) ;*gammaH(i,j)
endfor
endfor

endif


goto, repeatph
endif

endfor
close,1
openr,1,'output_virus.dat'

i=0
timetiter=fltarr(NT+11)
titer=fltarr(NT+11)
titers=fltarr(NT+11)

T=fltarr(NT+11)
Tdrop=fltarr(NT+11)
ttiter=fltarr(NT+11)


repeat begin
readf,1, xxx,sflu,ssars,xnsolid ,td,tt

t(I)=tt
tdrop(I)=td
timetiter(I)=xxx
titer(I)=sflu
titers(I)=ssars
i=i+1
endrep until eof(1)
nt0=i
close,1
i=0
ymin=100
itime0=-1




!P.POSITION=[0.2,0.15,0.95,0.95]
 plot_oi, timetiter,t, xrange=[xmin,xmax], $
  yrange=[285,310d0],charsize=1.5,xtitle='Time [s]', ytitle = 'T [K] ', thick=3
  oplot, timetiter,tdrop,  line=2,thick=3
  
ymin=1D-4
ymax=1D3
if ymin le 1D-10 then ymin = 1D-10
plot_io, timetiter,titer, xrange=[0,1E4], $
  yrange=[ymin,2d0],charsize=1.5,xtitle='Time [s]', ytitle = ' log!D10!N(titer numbber)', thick=3

ntc=NT0
oplot, timetiter(0:Ntc),titer(0:Ntc),line=0,thick=6,color=blueindex
oplot, timetiter(0:Ntc),titers(0:Ntc),line=0,thick=6,color=redindex

ymina=min(titer(0:NTC))
ymina1=min(titers(0:NTC))
if ymina1 le ymina then ymina=ymina1
	
ymaxa=max(titer(0:NTC))
timemaxa=max(ttiter(0:NTC))

plot_io, timetiter,titer, xrange=[0,timemaxa], $
  yrange=[ymina,ymaxa],charsize=1.5,xtitle='Time [s]', ytitle = ' log!D10!N(titer numbber)', thick=3
oplot, timetiter(0:Ntc),titer(0:Ntc),line=0,thick=6,color=blueindex
oplot, timetiter(0:Ntc),titers(0:Ntc),line=0,thick=6,color=redindex



close,1
openr,1,'output_vapour.dat'
pa=fltarr(NT)
pNH3=fltarr(NT)
pco2=fltarr(NT)
pHNO3=fltarr(NT)
pHcl=fltarr(NT)
pOA=fltarr(NT)
ppOA=fltarr(NT)

ppa=fltarr(NT)
ppNH3=fltarr(NT)
ppco2=fltarr(NT)
ppHNO3=fltarr(NT)
ppHcl=fltarr(NT)
time0=fltarr(NT)

close,2
openr,2,'output_partial.dat'

for i =0,nt-1 do begin
readf,1,aa,aww, paa,pnh3a,phno3a,pHCLa,pco2a,ttt,poaa
pa(i)=paa
pnh3(i)=pnh3a
phno3(i)=phno3a
phcl(i)=phcla
pco2(i)=pco2a
poa(i)=poaa


readf,2,aa,aww, paa,pnh3a,phno3a,pHCLa,pco2a, mfss,mfsll,awiin, poaa

ppa(i)=paa
ppnh3(i)=pnh3a
pphno3(i)=phno3a
pphcl(i)=phcla
ppco2(i)=pco2a
ppoa(i)=poaa

endfor
!P.POSITION=[0.2,0.15,0.95,0.95]
	ymin=1D-4
	ymax=10D3
 plot_oo, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[ymin,ymax],/nodata,charsize=1.5,xtitle='Time [s]', ytitle = ' CH!D3!NCOOH [ppbv]'
  oplot, time, pa, thick=5
  oplot, time, ppa, thick=5, color=redindex
  xx1=10
  xx2=4*xx1
  xx3=8*xx1
  y1= ymin* exp(.8*alog(ymax/ymin) )
  y2= ymin* exp(.9*alog(ymax/ymin) )
oplot, [xx1,xx2],[y1,y1], thick=5
oplot, [xx1,xx2],[y2,y2], thick=5, color=redindex
xyouts, .99*xx3,y1, 'Vapour pressure'
xyouts, .99*xx3,y2, 'Partial pressure', color=redindex



	ymin=1D-4
	ymax=1D3
 plot_oo, [1,2],[2,2], xrange=[xmin,xmax], $
yrange=[ymin,ymax],/nodata,charsize=1.5,xtitle='Time [s]', ytitle = '  NH!D3!N [ppbv]'
  
  oplot, time, pNH3, thick=5
  oplot, time, ppNH3, thick=5, color=redindex


  xx1=xmin*2
  xx2=4*xx1
  xx3=8*xx1
  y1= ymin* exp(.1*alog(ymax/ymin) )
  y2= ymin* exp(.2*alog(ymax/ymin) )

oplot, [xx1,xx2],[y1,y1], thick=5
oplot, [xx1,xx2],[y2,y2], thick=5, color=redindex
xyouts, .99*xx3,y1, 'Vapour pressure'
xyouts, .99*xx3,y2, 'Partial pressure', color=redindex
      ymin=1D-3
	ymax=60D0

 plot_oo, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[ymin,ymax],/nodata,charsize=1.5,xtitle='Time [s]', ytitle = ' HNO!D3!N [ppbv]'
  oplot, time, pHNO3, thick=5
  oplot, time, ppHNO3, thick=5, color=redindex

  xx1=xmin*2
  xx2=4*xx1
  xx3=8*xx1
  y1= ymin* exp(.7*alog(ymax/ymin) )
  y2= ymin* exp(.8*alog(ymax/ymin) )

oplot, [xx1,xx2],[y1,y1], thick=5
oplot, [xx1,xx2],[y2,y2], thick=5, color=redindex
xyouts, .99*xx3,y1, 'Vapour pressure'
xyouts, .99*xx3,y2, 'Partial pressure'


  
 plot_oo, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[1D-3,ymax],/nodata,charsize=1.5,xtitle='Time [s]', ytitle = '  HCl [ppbv]'
  oplot, time, pHCL, thick=5
  oplot, time, ppHCl, thick=5, color=redindex
xx1=xmin*2
  xx2=4*xx1
  xx3=8*xx1

y1= ymin* exp(.8*alog(ymax/ymin) )
  y2= ymin* exp(.9*alog(ymax/ymin) )

oplot, [xx1,xx2],[y1,y1], thick=5
oplot, [xx1,xx2],[y2,y2], thick=5, color=redindex
xyouts, .99*xx3,y1, 'Vapour pressure'
xyouts, .99*xx3,y2, 'Partial pressure', color=redindex


ymax=500
 plot_oo, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[1D-11,ymax],/nodata,charsize=1.5,xtitle='Time [s]', ytitle = '  Oxalic acid [ppbv]'
  oplot, time, poa, thick=5
  oplot, time, ppoa, thick=5, color=redindex
xx1=xmin*2
  xx2=4*xx1
  xx3=8*xx1

y1= ymin* exp(.8*alog(ymax/ymin) )
  y2= ymin* exp(.9*alog(ymax/ymin) )

oplot, [xx1,xx2],[y1,y1], thick=5
oplot, [xx1,xx2],[y2,y2], thick=5, color=redindex
xyouts, .99*xx3,y1, 'Vapour pressure'
xyouts, .99*xx3,y2, 'Partial pressure', color=redindex




	ymin=100
	ymax=1D4
 plot_oo, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[ymin,ymax],/nodata,charsize=1.5,xtitle='Time [s]', ytitle = ' CO!D2!N [ppmv]'
  oplot, time, pco2, thick=5

    oplot, time, ppco2, thick=5, color=redindex,line=2

  xx1=xmin*2
  xx2=4*xx1
  xx3=8*xx1
  y1= ymin* exp(.1*alog(ymax/ymin) )
  y2= ymin* exp(.2*alog(ymax/ymin) )

oplot, [xx1,xx2],[y1,y1], thick=5
oplot, [xx1,xx2],[y2,y2], thick=5, color=redindex
xyouts, .99*xx3,y1, 'Vapour pressure'
xyouts, .99*xx3,y2, 'Partial pressure'


; plot Na+ and Cl- concnetration in particles relative to initial Concentration

yna=fltarr(nt+1)
ycl=fltarr(nt+1)
ynat=fltarr(nt+1)
yclt=fltarr(nt+1)

xnp=fltarr(NP+1)
close,1
openr,1,'output_shells_n.dat'

for i=0,nt do begin

xna=0d0
xcl=0d0
xnat=0d0
xclt=0d0



for  j=0, ns2(i)-1 do begin
readf,1,timee,xnp

xna=xna+xnp(16)
xcl=xcl+xnp(17)
xclt=xclt+xnp(17)+xnp(np)
xnat=xnat+xnp(16)+xnp(np)+xnp(np-1)*2+xnp(np-2)

endfor

if i eq 0 then xna0= xna
if i eq 0 then xcl0= xcl
yna(i) = xna/xna0
ycl(i) = xcl/xcl0
ynat(i) = xnat/xna0
yclt(i) = xclt/xcl0


endfor
close,1

ymin= .9* min(ycl)
ymax=1.01
xmin=1D-3
 plot_oI, [1,2],[2,2], xrange=[xmin,xmax], $
  yrange=[ymin,ymax],/nodata,charsize=1.5,xtitle='Time [s]', ytitle = 'relative # of Cl!U-!N Na ions'
  oplot, time, yNa, thick=5,line=0,color=redindex
  oplot, time, yNat, thick=7,line=2,color=redindex
  oplot, time, ycl, thick=5,line=0
  oplot, time, yclt, thick=7,line=2

  xx1=2*xmin
  xx2=4*xx1
  xx3=8*xx1
  y1= ymin* exp(.1*(ymax-ymin) )
dd=.07
oplot, [xx1,xx2],[y1,y1], thick=5, line=0,color=redindex
xyouts, .99*xx3,y1, 'Na!U+!N in liquid',color=redindex

y2= y1* exp(dd*(ymax-ymin) )
print,y2

xyouts, .99*xx3,y2, 'Na!U+!N in liquid and solids'
oplot, [xx1,xx2],[y2,y2], thick=7, line=2, color=redindex

y2= y2* exp(dd*(ymax-ymin) )
xyouts, .99*xx3,y2, 'Cl!U-!N in liquid'
oplot, [xx1,xx2],[y2,y2], thick=5, line=0

y2= y2* exp(dd*(ymax-ymin) )
xyouts, .99*xx3,y2, 'Cl!U-!N in liquid and solid'
oplot, [xx1,xx2],[y2,y2], thick=7, line=2



IF ( !D.NAME EQ 'PS') THEN DEVICE,/CLOSE
end 







