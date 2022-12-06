c     The main program ResAM

      
c     Input files: input.dat, ML.dat and species.dat
c     for the electrodynamic balance (EDB) simulation (imode_output=2) + input_eff.dat 
c     An example is provided in folder RH50_1um_v1.0
      
      
c     input.dat
      

c line 1: comment
C line 2: radius in cm
c line 3: Pressure in hpa           
c line 4: 5 , 0 
c          number of shells, 
c          eddy diffusion coefficient in cm2/s

c line 5: 0 
c         output mode
c         0: (log time), pH , initialize with ML.dat 
c         1: EDB linear time, no pH calculation initialize with ML.dat 
c         2: linear output time
c         3: not used
c         4: logtime, for initial composition, can adjust 19 X+ and 20 Y- for the desired pH /so
c         5: as 1, but first output time 1E-5 s
c         10+ 0-5, as 0-5, but solid only in the center (ES&T version)
      
      
c     line 6 mode for pH calculation
c     0: default , precise calculation of pH
c     1: many buffers are present
c     2: for the EDB modelling the pH is assumed to be constant (pH=7) 
      


c line 7 to 7+N
c          time in seconds
c          relative humidity (RH, 0.5 = 50%)
c          temperature in K
c          acetic acid in the gas phase in ppb
c          NH3 in the gas phase in ppb
c          CO2 in ppm
c          HNO3 in ppb
c          HCl in ppb : > 0 or =0, evaporation considered ; < 0 : no HCl evaporation

c     species.dat
      
      

      implicit real*8 (a-h,m, o-z)
      
      parameter (NSMM=1000)
      real*8 xn1(NSMM),xn2(NSMM), M1,M2,x(NSMM+1),xns(NSMM)
      parameter (np=29)         ! number of species
      character *60 name(np), filename
      character *200 text

      real*8 dl_factor2(2,np),dl_factor(NP)
      common /DL/ DL_factor2 ,deltaxgas,jmin,nmin
      
      common  /kout/xfcn2

      real*8 awshell(NSMM),sNaClshell(NSMM)
      common /awin/ awin
      common/solid/xnsolid
      common /output/imode_output,idiff,imode_pH
      


c     species, given by species.dat
      
c     1:  H2O
c     2:  organic 460 molar mass
c     3:  acetic acid total
c     4:  sum of NH3, NH4, and NH4CH3COO
c     5:  sum of CO2 + HCO3-, and CO3--
c     6:  H+
c     7:  OH-
c     8:  CH3COO-
c     9:  undisocciated acetic acid
c     10: NH4CH3COO
c     11: organics which can form protein crystals 
c     12: NH4+
c     13: CO2(aq)
c     14: CO3--
c     15: HCO3-
c     16: Na+
c     17: Cl-
c     18: NO3-
c     19: other cations than H+, Na+, and NH4+ 
c     20: other anions than Cl-, HCO3-, CO3--, CH3COO-, NO3-
c     21: Virions inactivation with aH+   Influenza virus
c     22: Virions inactivation with aH+   SARS-CoV-2
c     23: NH3(aq)
c     24: H3PO4
c     25: H2PO4-
c     26: HPO4--
c     27: HSO4-
c     28: SO4--
c     29: solid NaCl (the solid shell is put below the liquid, it its assumed that the liquid diffusion is only affected by the presence of the solid due to the size effect ,i.e. the gradient).

      logical ex

      real*8 ntiter(NSMM)
      real*8 flux(NSMM,np+1)
      real*8 f(NSMM)
      real*8 f3(NSMM)
      real*8 f4(NSMM)
      real*8 ml(NP) ! molality of the species
      real*8 m(NP)  ! molality of the species
      real*8 xn(NSMM,NP),  xnnew(NSMM,NP)
      real*8 ml6shell(NSMM)
      real*8 xhshell(NSMM),gammaclshell(NSMM)
      common /awshell/ awshell,xhshell,ml6shell,sNaClshell,gammaclshell
      
      real*8 mm(NP) ! molar mass
      real*8 mv(NP) ! molar volume
      integer IZC(NP),ishell(NSMM) !charges of all ions taken from species.dat
      common /m/ mm, mv, izc

c     xn1: moles of H2O
c     xn2: moles of sucrose
c     xn3: moles of acetic acid
c     xn4: moles of NH4 
      real*8 xn3(NSMM),xN4(NSMM)
      parameter (Ntimes=80000)
            
      real*8 timetr(ntimes),rhtr(ntimes),time1(1),temp1(1)
      real*8 ttr(ntimes) , timeliq(ntimes)
      real*8 timetr0(ntimes),rhtr0(ntimes),gash2otr0(ntimes)
      real*8 ttr0(ntimes),gasaatr(ntimes),gasamtr(ntimes)
      real*8 gasaatr0(ntimes),gasamtr0(ntimes)
      real*8 gasHNO3tr0(ntimes),gasco2tr0(ntimes)
      real*8 gasHCLtr0(ntimes)
      real*8 y1(1),outputtime(ntimes)

      common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &     partvapHCl
      
      real*8 MLL(25000), awL(25000), gcl(25000),gna(25000)
      common /lookup/ MLL,awL,gcl,gna
      real*8 timeEff(100),solidm(100)
      common /enhance/radius,r1,enh_factor
      common /Ienhance/Ienh,iscenter
      
      common /timeeq/ timeeq
      common /pi/pi

      real*8 xhhn(5001)
      real*8 fluxs(NP)
      real*8 times

 
      N1=1
      iscenter=0
      solidsm=0
      Nxhh=1
      loop=0
      Ieffcycle=1

      sigma = 72     ! surface tension erg/cm2 of water
      aweff = 0.569  ! efflorescence water activity (not valid for output mode 2 and 12), this value can be modified if the NaCl effloresces at a different aw
      pi=dacos(-1d0)
      xnsolid=0d0
      xdissolid=5E-4  ! we set the miniumum distance between the solid NaCl shells to 5 um
c     If the NaCl dendrite structure is larger or smaller than 5 um, one can modify this value here.

      
c     ------------------------------------------------------------------------------------------------------------------------------------------------      
c     Read input data
c     ------------------------------------------------------------------------------------------------------------------------------------------------      
       
      open(1,file='species.dat')
             read(1,'(A)') text 

      DO I=1,NP
             read(1,'(A)') text 
             do j=1,30
             if (text(J:J) .eq. ',') goto 591
             enddo
          print*, 
     & 'After the name of the species, a comma should follow!'
             stop

 591         filename=text(1:J-1) 
              text= text(j+1:200)

          read(text,*) ii , MM(i),
     & mv(i), dl_factor2(1,i) ,izc(i)
          dl_factor2(2,i)=dl_factor2(1,i) ! 

          write(6,'(A20,I5, 3F15.4,I5)') filename,i,MM(i),
     & mv(i), dl_factor2(1,i) ,izc(i)

             enddo
             close(1)
             
      open(1,file='input.dat')
      read(1,'(A)') filename
      read(1,*) r0
      read(1,*) press0
      read(1,*) NSmax,Deddy
     

      read(1,*) Imode_output       ! output mode , 0: log(time), 1: linear time
      iscenter=0
      if (imode_output.ge.10) then
         iscenter=1
         Imode_output=Imode_output-10
      endif
      
      imode_outputt=      imode_output  
      if (imode_output.ge.3)  imode_outputt=0
      
      ismode=0
      if (imode_output.eq.0) ismode=1
      if (imode_output.eq.1) ismode=1
      if (imode_output.eq.2) ismode=1
      if (imode_output.eq.4) ismode=1
      if (imode_output.eq.5) ismode=1
      if (ismode.eq.0) then
      print*, 'Input a valid output mode (0,1,2,4,5,10,11,12,14,15)'
         
         stop
      endif
      

      
       if(imode_output .eq.2 ) then

          open(11,file='input_eff.dat')
      read(11,*) xdissolid
          DO I=1,100
          read(11,*,end=26) xxx,timeeff(I),solidm(I)
      enddo
         endif
 26     continue

      read(1,*) Imode_ph

      if (imode_ph.lt.0 .or. imode_pH.gt.2) then
         print*, 'Input a valid pH mode (0,1,2) '
         stop
      endif


c     ------------------------------------------------------------------------------------------------------------------------------------------------     
c     For sensitivity calculations, one has the change the Idiff here
c     ------------------------------------------------------------------------------------------------------------------------------------------------

c     Idiff ist set here to 0 
      idiff=0
c line 6: 0
c         diffusion mode Idiff
c         0: activity gradient, with charge balance except (H+,OH-), SLF diffusivity
c         1: activity gradient,       SLF diffusivity, w.o. charge balance
c         2: concentration gradient,  SLF diffusivity, w.o. charge balance
c         3: concentration gradient,  sucrose diffusivity, with charge balance
c         4: concentration gradient,  citric acid diffusivity, w.o. charge balance.
c         5: Dl_ions = D_H2O_walker     activity with charge balance 
c         6: Dl_ions = D_H2O_thiswork   activity with charge balance


        Iexit=1 ! the index when the time=0, the exhaled air exits the nose or mouth
c     For the mixing, the air mixed air at time=0 with the ambient air.
        

        DO I=1,ntimes
         read(1,*,end=222) timetr0(I), rhtr0(I),ttr0(I),
     &  gasaatr0(I),gasamtr0(I),gasco2tr0(I),gashno3tr0(I),gashCLtr0(I)
         T=ttr0(I)
         gash2otr0(I) = vwater(t)*rhtr0(I)
         if (dabs(timetr0(I)).le.1D-10) Iexit=I
c         print*,i, timetr0(I), gash2otr0(I), RHtr0(I)

      enddo
 222     Ntr0=I-1
         print*, ntr0
         

         close(1)
       
       
c     test diffusion coefficients in the liquid for EDB simulation

cccccccccccccccccccccccccccccccccccccaccccccccc      

         if (imode_output.eq.2) then
               t=293.15
        DO aw =0.,1., 0.01
       call cal_dlaw(T,aw,dl)
       ddcl=dl_factor2(1,17)
       dd=dl_factor2(1,16)
       call cal_dlaw_walker(T,aw,dlw)
       call cal_dlaw_suc(T,aw,dlsuc)
       call cal_dlaw_citric(T,aw,dlcitric)
       call cal_dlaw_walker_mod(T,aw,dlm)
       write(52,'(15E15.6)')aw,dlm, dl*dd, dl*ddcl,dlcitric,DLSuc,dlw
      enddo
      DO xw=0d0,1.001,.01d0
       aw = aw_pinene(xw)
       call cal_dlaw_piene(T,aw,dl,xw)
       write(53,'(15E15.6)') aw,dl
       enddo
       endif

cccccccccccccccccccccccccccccccccccccccccccccc      

      

c     if time > timeeq, it is assumed that the water activity of the droplet is equal to the RH
c     no kinetic calculation for H2O

       
        timeeq=1d10
c     H2O Equilibrium for exhaled aerosol

        if (imode_output.ne.2) then 
        if (r0.le. .55D-4) then
           timeeq=0.5
           endif
        if (r0.le. .022D-4) then
           timeeq=0.3
           endif
           if (r0.ge. 0.3E-4)then
              timeeq=40* (r0*1D4)**2
              endif
              endif




c     define output times
       time11=1D-2
       if (imode_output.eq.5)        time11=1D-5
       if (imode_output.eq.5)      imode_output=0

       time22= timetr0(NTR0)
        Noutput=1001
       DtimeL = dlog(time22/time11)/(Noutput-1)
      
       


       
      if ( imode_output .eq. 0 .or. imode_output .eq. 4 )then

         if (TIMETR0(1).ge.0d0) then

         Noutput=1001
         do i=1,nOUTPUT
         outputtime(I)= time11* dexp((I-1)* dtimeL)
         enddo

         else
         Noutput=1101
         dtime= (TIMETR0(2)-TIMETR0(1))/99

         do i=1,100
          outputtime(I)= timetr0(1)+ dtime*(I-1)
c          print*, i, outputtime(I)
          enddo


         do i=1,nOUTPUT
         outputtime(I+100)= time11* dexp((I-1)* dtimeL)
         enddo
         endif

c         if (i.le.100)print*, I, outputtime(I)

         endif



      if ( imode_output .eq. 1   ) then

       DtimeL = (time22-time11)/(Noutput-1)
         do i =1,100
            outputtime(I)= TIMETR0(1) + (i)*DTIME
         ENDDO
         endif


            if (imode_output .eq. 2  ) then
c               timetr0(1)=timetr0(2) -1d0
               Noutput = Ntr0
               DO I=1,Ntr0
             outputtime(I)= timetr0(I)
               enddo
                           endif

             outputtime(Noutput+1)=  outputtime(Noutput)+100
                print*, noutput
          



c     ----------------------------------------------------------------------
c     Initialization  
c     ----------------------------------------------------------------------

        t0=ttr0(1)
        Ta=T0
        tdrop=T0
        T= t0
        rh=rhtr0(1)
        Amisch = apNaCL(t)
        press=press0

c     read the inital compostion (e.g., given from the SLF recipe or equilibrated with the given gas phase)
        
         open(1,file='./ML.dat')
         DO I=1,NP
         read(1,*) ii,ml(I)
         enddo
         close(1)
         M=ML

c     adjust the pH for mode 4  ( Adjust X+ and/or Y-) to achieve the required pH at a given gas phase, e.g. pH 6.6 of 2.5% CO2 and 135 ppb NH3)

         if (imode_output .eq. 4) then
c     calculate only one shell
            NSmax=1
            read(55,*) d19
            read(55,*) d20
            ML(19) = ml(19)*(1+d19)
            ML(20) = ml(20)*(1+d20)
            endif

         
c     find the efflorescence aw and Apeff (activity product at which efflorescence happens)
           call calHNew(T,m)
           call aw_back
     & (t,M,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

c     if the initial solution is dilute
           if (aw.ge.aweff ) then
                    ff=1.001
        DO xx = 1,1D10,1

           DO J= 2,NP
              M(J)=m(J)*ff
              enddo
           call calHNew(T,m)
           call aw_back
     & (t,M,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           app= m(16)*m(17)* gammacl*gammaNa
           
        if (aw.le.aweff) goto 44
        
        
      enddo
 44     print*, 'aw = ', aw
        else

c     if the initial solution is concentrated (e.g. at RH 50%)

                    ff=1.001
        DO xx = 1,1D10,1

           DO J= 2,NP
              M(J)=m(J)/ff
              enddo
           call calHNew(T,m)
           call aw_back
     & (t,M,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
           app= m(16)*m(17)* gammacl*gammaNa
           
        if (aw.ge.aweff) goto 144
        
        
      enddo
 144     print*, 'aw = ', aw

         DO I=1,NP
            print*, I, M(I)
            enddo

        endif

         APEFF= app  ! activity product at which efflorescence happens
         print*, 'ap_eff = ', apeff
         print*, 'ap_NaCl = ', apnacl(T)
         if (apeff.le. apnacl(T)) apeff=1.1* apnacl(I) ! assuming minimum S_NaCl of 10%



c     initializing xn: the moles of each species in each shell

c     set equal volume in each shell
         call set_shells_vol(T,RH,r0,ML,xn,x,NS,nsmax)

         stiter0=0
         DO J=1,Ns
            stiter0=stiter0+xn(j,21)
          enddo
            
         vv=4*pi/3*x(2)**3
         ctiter0=xn(1,21)/vv ! initial virus concentration


      partvap3=gasaatr0(1)*press0*1D-9 ! hPa
      partvap4=gasamtr0(1)*press0*1D-9 ! hPa
      partvapco2=gasco2tr0(1)*press0*1D-6 ! hPa
      partvapHNO3=gasHNO3tr0(1)*press0*1D-9 ! hPa
      partvapHCL=gasHcltr0(1)*press0*1D-9 ! hPa
      

       time=timetr0(1)
       print*, ' ns  =  ', ns,time



              NOUT=1
           print*, ' The first output may take several minutes to 1 h'


c     Here read the data for restart  
c     open output file

c     Molalities of all species in each shell
              open(3,file='output_shells_M.dat')

c     Moles of all species in each shell
              open(23,file='output_shells_n.dat')

c     Infectious number of virus etc
              open(33,file='output_virus.dat')
c     partial pressures
              open(15,file='output_partial.dat')

              open(14,file='output_vapour.dat')


c     continue the previous calculation if output_virus.dat exists
c     for a recalcualtion from the beginning, delete output_virus.dat

          INQUIRE (FILE='output_virus.dat', exist=ex)

           if (ex .and. imode_output.ne.4) then
                 DO Ii=1,5000
                  read(33,*,end=19)   
     &     time,sflu,ssars

                  read(14,*,end=19)timee
                  read(15,*,end=19)timee

                  enddo

 19                continue
                   Nout=Ii-1

                   print*,time, outputtime(nout), outputtime(nout+1)
                   
                   print*, 'NOUT ', nOut
               DO Ii=1,NOUT
              read(3,*,end=155) time,x(NS+1),ns,x(1)

              DO I = 1,ns
              read(3,*,end=155) jj,
     &        x(I+1)

              enddo
              

              DO I = 1,ns
              read(23,*,end=155) Iii,time11, (xn(i,j),j=1,NP)
c              print*, ii, I,ns
              enddo

            enddo
              goto 156
 155          Nout=II-1
              print*,'NOUT =', Nout

 156          continue

              
c     For a restart, one can reduce the number of shells, but not increase!
c     It is useful to reduce the number of shells when the equilibrium is reached and the calculation is slow due to high number of shells and one want to reach longer times (months - years)
              
              
              if (NSmax.gt. nS) then
                 print*,  ' The number of shells is higher than '
                 print*,  ' the previous calcaltion, stop! '
                 stop
              endif

c     group to 1 bin
              if(NSmax.eq.1 .and.NS.gt.1) then

                 DO I=2, NS
                    DO J=1,NP
                       xn(1,j) =xn(1,j) + xn(I,j)
                    enddo
                    enddo

                 NS=1
                 endif
     
c     reduces shells  to NSmax shells
                 

                 if (NSmax.le.NS/2) then
     
                    DO I=1,NSmax
                          DO J=1,NP
                             xnnew(i,j)=0d0
                          enddo
                          enddo

                    N = NS/NSMAX
                    Nrest= NS- N*NSMAx

                    N2=0
                    DO I= 1, NSMAx
                       NN=N
                       if (I.le.Nrest)  NN=N+1
                       N1= N2+1
                       N2= N2+NN
                       print*, N1,N2
                       DO II=N1,N2
                          DO J=1,NP
                             xnnew(I,j) =xnnew(I,j)+ xn(Ii,j) 
                          enddo
                       enddo
                    enddo

                    NS=NSMAx

                    DO I=1,NSmax
                          DO J=1,NP
                             xn(I,J)=xnnew(i,j)
                          enddo
                          print*, I, xn(I,1),xn(I,16)
                          enddo
                    endif


c     group outer  bins into one bin, when  2/NS < Nsmax < NS
                    
                    II= NS- NSMAX  
                    if ( II.lt.NSmax .and. II.ge.1) then

                       NG=II+1

                       DO I=NS-NG+2, Ns
                          II=NS-NG+1
                             DO J=1,NP
                                xn(II,J) =xn(II,J) +xn(I,J) 
                                enddo
                             enddo

                             endif

                          NS=NSMAx
                 
               DO I = 1,ns
c                 write(6,'(I6,501E14.5)') Iii, (xn(i,j),j=1,NP)
                 enddo

                 print*,' ns= ', ns
                 print*,' time= ', time
                 
               NOUT=NOUT+1

               
               DO I = 1,Noutput
                  if (outputtime(I).gt. time) goto 119
               enddo
 119           print*, I, nout

               nout=I


             endif


cccccccccccccccccccccccccccccccccccccccccccccccc
c     start of the time integral
cccccccccccccccccccccccccccccccccccccccccccccccc
             
 11           continue


          time1(1)=time

c     define the partial pressure
          
c     interpolates the T

c     if Eddy Diffusion coefficient is zero, then make a linear interpolation of input data

          if (Deddy.le.1D-10 ) then
c     T
                  time1(1)=time
             call intpl(timetr0, ttr0, ntr0, time1, y1, n1)          
         t=y1(1)
c     RH
         call intpl(timetr0, gash2otr0, ntr0, time1, y1, n1)          
         rh=y1(1)/vwater(T)
c   
         call intpl(timetr0, gasaatr0, ntr0, time1, y1, n1)          
         ppbace=y1(1)
         call intpl(timetr0, gasamtr0, ntr0, time1, y1, n1)          
                  ppbnh3=y1(1)
               call intpl(timetr0, gasco2tr0, ntr0, time1, y1, n1)          
                  partvapco2=y1(1)*1D-6*press0
                  call intpl(timetr0, gashno3tr0, ntr0, time1, y1, n1)          
                  partvaphno3=y1(1)*1D-9*press0  
                  call intpl(timetr0, gashcltr0, ntr0, time1, y1, n1)          
                  partvaphcl=y1(1)*1D-9*press0   
               else
c     take the eddy diffustion coefficient of time 0 
                  
                  A0= 2*pi*0.75**2
                  At = a0+ 4 * Deddy *time
                  if (time.le.0) At= A0
                  
            t = a0/At * ttr0(1) + (At-a0)/At*ttr0(Ntr0)
              pwater = a0/At * gash2otr0(1)
     & + (At-a0)/At*gash2otr0(ntr0)
              RH=pwater/vwater(T)

         ppbace = a0/At * gasaatr0(1) + (At-a0)/At*gasaatr0(Ntr0)
         ppbNH3 = a0/At * gasamtr0(1) + (At-a0)/At*gasamtr0(Ntr0)
         xxx = a0/At * gasco2tr0(1) + (At-a0)/At*gasco2tr0(Ntr0)
         partvapco2= xxx*1D-6*press0


         xxx = a0/At * gasHNO3tr0(1) + (At-a0)/At*gasHNO3tr0(Ntr0)
         partvaphno3= xxx*1D-9*press0         

         xxx = a0/At * gasHcltr0(1) + (At-a0)/At*gasHcltr0(Ntr0)
         partvaphcl= xxx*1D-9*press0         
         endif


       if (rh.le. .01d0) rh=.01d0 !not dryer than 1 %
         partvap3=ppbace*1D-9*press0
         partvap4=ppbnh3*1D-9*press0
              press=press0
              xnsolid=0d0
              DO I=1,NS
                 xnsolid=xn(I,np)+xnsolid
              enddo
c     end define partial pressures

              
c     merge shell with less molecules to the neighboring shells
          call merge(xn,x,NS)
              

          if ( xnsolid.le.1D-30) then

c  reset bins when one shell becomes too thin
c  do rebin when no solid is present

         call  reset_shells(NS,x,xn)


c     check efflorescence
                  xapmax=0
             xapmin=1D10
               DO JJ = 1,NS
               if (xapmax.le. sNaClshell(jj))xapmax=sNaClshell(jj)
               if (xapmin.ge. sNaClshell(jj))xapmin=sNaClshell(jj)
                enddo


                isnuc=0

c     nuc when ap > ap_nuc if activity produce of NaCl > AP_eff, NaCl nucleates ( isnuc=1) for non EDB

         if ( xapmax .ge. Apeff .and. 
     &  imode_outputt.eq.0)  isnuc=1


c     nuc at given time for EDB simulation
         if(time.ge.timeeff(Ieffcycle).and.imode_output .eq.2 ) isnuc=1
c                       print*, 'ieffcyle' , ieffcycle

                       if (isnuc.eq.1) then

            xnsolid=0d0
            
            ML(1)= 1000/MM(1)


            DO J =1,NS
               DO K=2,NP
                  ml(k)=ml(1)* xn(j,k)/xn(j,1)
               enddo
               M=ML 
        
        call calHNew(Tdrop,mL)
        call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa,
     & gammas1,gammas2)        

        DO K=1,NP
c           print*,k,ml(k)
        enddo
        ap= ML(16)*ml(17)*gammaCl*gammaNa
        print*, 'ap ', ap
        
c     do solid only if 
        if (ap.gt. apnacl(tdrop)) then
               
               call  effl_misch(Tdrop,mL,msalt)               

               print*, 'msalt', msalt
               
               
                  xnsolidJ = Msalt/M(16)*xn(j,16)
                  print*,j, xnsolidj
                  xnsolidJ = xnsolidj*(1+solidm(Ieffcycle))
                  print*, 'solidm = ',solidm(Ieffcycle)

                  print*, xnsolidj, solidm(Ieffcycle),Ieffcycle
                  print*, 'time' , time, timeeff(Ieffcycle)


 
              xnsolidjcl =xnsolidJ

              xn(j,16) =xn(j,16) -xnsolidj !mol
              xn(j,17) =xn(j,17) -xnsolidj  !mol


              xnsolid=xnsolid+xnsolidj

              ML(16) = ML(1) *xn(J,16)/xn(j,1)
              ML(17) = ML(1) *xn(J,17)/xn(j,1)

              
c     test activity products
        call calHNew(Tdrop,mL)
        call aw_back
     & (tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        

        ap= ML(16)*ml(17)*gammaCl*gammaNa


        endif

                    enddo
                       Ieffcycle=      Ieffcycle+1


                       print*, time, ' NaCl efforesed !',NS,xnsolid

                       if(iscenter.eq.1) then
                          xn(1,np) = xnsolid
                          goto 228
                       endif

             if(iscenter.eq.0) then
c     distribute the solid into shells
c     with a spacing of xdissolid
                
                xx=0d0
                
c     the first shell
           JJ=1
                ishell(JJ)=1
c     out the solid only in shell with s > 1
                
                I=2
                x11= x(2)
                ss=x11**2
 209            continue
                    dx =  (x(I+1)-x(I))
                write(6,'(I5,5E15.6)')
     & I ,(x(I+1)-x11)*1D4, dx*1D4,(xdissolid*jj - dx)*1D4


                if (x(I+1)- x11 .ge. xdissolid  ) then

                    if (snaclshell(I).ge. .8*apnacl(T))then
                    JJ=JJ+1

                    print*,' jj ', jj,I
                   x11=x(I+1)
                    
                    ishell(JJ)=I
                   ss=ss+ x(I+1)**2
                  
                endif
                   endif
                   I=I+1
                   if (I.LE.NS) goto 209

                   
                   
             print*,'solid and total bins ', jj,NS



                snn=0
             DO I=1, JJ
                II= ishell(I)
                xn(II,NP) =xnsolid/ss*x(II+1)**2
                snn=1/ss*x(II+1)**2
                enddo

c     distributes the molecules weighted with the area
                print*, ' normalize factoer =1 ', snn
                

                    DO I=1,NS
        write(6,'(I5,6E15.6)')i, xn(I,np),ml(1)*xn(I,16)/xn(I,1)
     & ,ml(1)*xn(I,17)/xn(i,1)
                       enddo
                       endif !iscneter=0
                       
                       
 228                   continue
        endif
        endif
ccccccccccccccccccccccccccc
c     End efflorescence
cccccccccccccccccccccccccc
        loop=loop+1

c     calculate the time step and fluxes
          call cal_flux(time,xn,x,flux,dtime,NS)
          Tdrop=ta              ! the droplet T is calculated in cal_flux
          



        xmw=0

        xms=0d0
        
        DO I=1,NS
           xmw=xmw+xn(I,1)*mm(1)
           xms=xms+ xn(I,2)* MM(2)
        DO J=6,Np
           xms=xms+ xn(I,J)* MM(J)
        enddo
        enddo
c     NaCl solid in bin NP
        
        if (loop.eq.5000) then
c     diagnostic output finds the time limiting shell and species
              write(26,'(2I7,10E15.6)')nmin, jmin, time,dtime
              loop=0
              endif



      rad=x(NS+1)



              
       if (time+dtime.gt. outputtime(Nout)) 
     & dtime= outputtime(NoUT)-time+1D-12

       time=time+dtime
       times=times+dtime

          time1(1)=time



c     calculates the infectious titer of viruses

          DO I=1,ns
             vv=4*pi/3 *(x(I+1)**3-x(I)**3)

c             pH= -dlog( ML6shell(I))/dlog(10d0) !pH
             pHa= -dlog( xhshell(I)*ML6shell(I))/dlog(10d0) !pH

             aww=awshell(I)

             if (time.ge.0d0) then
c     mass ratio of liquid
              xm0=xn(I,NP)* MM(NP)+xn(I,2)* MM(2)
                DO kk =6,20
                xm0=xm0+ xn(I,kk)* MM(kk)
                enddo

                xm0= xm0/(xm0+ xn(I,1)*mm(1))

                taua=tau_ivea(pha,xm0)

                if (time.le.0) taua=1D10
              xn(I,21)=xn(I,21)  * dexp(-dtime/taua)

             call cal_tau_sars(T,aww,pHa,tauSa)
                if (time.le.0) tausa=1D10
             xn(I,22)=xn(I,22)  * dexp(-dtime/tausa)
          endif
          
          enddo
                
          
          




c     calculates the new composition after the timestep dtime
c     for dissociation species, the liquid phase diffusion is set to the sum of the individual species.


		
 
          DO J=1, Np-1
             if (J.eq.5 .or. J.eq.13 .or. J.eq.15) then
                fluxs(J)=flux(ns,j)*dtime+ fluxs(J) ! for CO2 take the flux from shell ns-1 to ns
                else
              if (dl_factor2(1,j).gt.0) then
                fluxs(J)=flux(ns+1,j)*dtime+ fluxs(J)
                endif

             endif

              if (dl_factor2(1,j).gt.0) then
                DO I=1,ns
                  df=(flux(I+1,j)-flux(I,j))
                  if (J.eq.5 .and. I.eq.NS) df=0d0

                  xn(I,j)=xn(I,j)+ df*dtime
                 if (xn(I,j).lt.1D-40) xn(I,j)=0d0
              enddo
          endif
       enddo


c     consider the Na and Cl flux to  the crystal 

       DO I=1,ns

          if (xn(I,NP) .gt.0d0) then
c          write(6,'(1I5,5E15.6)'), I, dtime,xn(I,np), flux(I,np),
c     & flux(i,np+1), dtime

          xnss=xn(I,NP)+(flux(I,np+1)-flux(I,np))*dtime

                 if (xnss.lt.0d0) then
                    fff= -xn(I,NP) /dtime
               flux(I,np+1) = fff 
               flux(I,np) = 0
                    endif


c     subtract f2 the Na Cl from liquid layer 
c     

c     reduce flux  that xn16 
          xss=xn(I, 16) - flux(I,np+1)*dtime
          if (xss.lt.0) then
             flux(I,np+1) =  xn(I, 16)/dtime
             endif

          xss=xn(I, 17) - flux(I,np+1)*dtime
          if (xss.lt.0) then
             flux(I,np+1) =  xn(I, 17)/dtime
             endif

c     add Na Cl to the liquid shells
        xn(I, 16) =xn(I, 16) - flux(I,np+1)*dtime
        xn(I, 17) =xn(I, 17) - flux(I,np+1)*dtime


c     subtract  f1 the Na Cl from liquid layer below
        if (I.ge.2) then
c     not negative
              xss=xn(I-1, 16) - flux(I,np)*dtime
              if (xss.lt.0) then
              flux(I,np) =  xn(I-1, 16)/dtime
             endif

c     not negative
               xss=xn(I-1, 17) - flux(I,np)*dtime
             if (xss.lt.0) then
             flux(I,np) =  xn(I-1, 17)/dtime
             endif

              xn(I-1, 16) =xn(I-1, 16) - flux(I,np)*dtime
              xn(I-1, 17) =xn(I-1, 17) - flux(I,np)*dtime

              endif


          xn(i,np)=xn(I,NP)+(flux(I,np+1)+flux(I,np))*dtime
          if (xn(i,np).lt.0d0) xn(i,np)=0d0

           endif !xn(I,np) > 0

       enddo

c     consider the Na and Cl flux to the crystal 
          
      

c     CO2  equilbrium in shell the outermost shell

       if (partvapco2.ge.1D-40) then
         xkelvin = dexp( 2* 85d0 * MV(5) /(8.314E7*T*x(NS+1)) )
           pp=partvapco2/xkelvin

c     calculates the CO2, HCO3-, CO2- concentration at given pCO2 partial pressure in the outermost shell
       DO J=2,nP
          ML(j)=xn(NS,J)/xn(NS,1)* ml(1)
       enddo
           call calhnewco2(Tdrop,ML,pp)
           call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2)


        xn(NS,5) =  ML(5)/ML(1)* xn(ns,1)
        xn(NS,13) =  ML(13)/ML(1)* xn(ns,1)
        xn(NS,14) =  ML(14)/ML(1)* xn(ns,1)
        xn(NS,15) =  ML(15)/ML(1)* xn(ns,1)
       else
        xn(NS,5) =0d0
        xn(NS,13) =0d0
        xn(NS,14) =  0d0
        xn(NS,15) =  0d0
        endif

       pco2eq=pco2 ! for output

c     begin output
      if (time.ge.outputtime(Nout) ) then

       xhs=0

c     write titer concentration


        write(3,'(1F17.7,1E15.6,I6,1E16.6,I5)') time,x(NS+1),ns,x(1),NP

          write(6,'(I10,1E15.6,1E15.6,1F10.4,1E15.7)') 
     *         Nout,TIME,dtime,rh, awshell(NS)
          
          

          
             DO J=1,NP
             xns(J)=0
             enddo
             sna=0
             sco2=0
             sx=0
             scl=0
             sy=0
             xn11=0
             sflu=0
             ssars=0
             
             DO I=1,NS

                sna=sna+xn(I,16)
                scl=scl+xn(I,17)
                sco2=sco2+xn(I,5)
                sx=sx+xn(I,19)
                sy=sy+xn(I,20)


                xn11=xn11+xn(I,11)
                                sflu =sflu +xn(I,21)
                                ssars=ssars+xn(I,22)

                ML(1)=1000d0/18d0


             DO J=1,NP
                ML(J)=xn(I,j)* ML(1)/xn(I,1)
                xns(J)=xns(J)+ xn(I,J)
             enddo
c                      call ml2m(ML)
                      M=ML

       call vapnew(Tdrop,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2)
              call calHNew(Tdrop,mL)
      call aw_back(tdrop,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
       Smisch=ML(16)*ML(17)*gammaNa*gammaCl/Apnacl(tdrop)
       phh = dlog(ML(6)*gammah)/dlog(.1d0)

       DO J=6,NP
          xn(I,j) = ML(J)/ML(1)*xn(I,1)
       enddo
       
       ppb=pacetic/press * 1D9
       vv=4*pi/3 * (x(I+1)**3-x(I)**3)


              flu=xn(i,21)/vv/ctiter0

              if (flu.le.1D-20) flu=1D-20
              sarsa=xn(i,22)/vv/ctiter0
              if (sarsa.le.1D-20) sarsa=1D-20
              
              

              M=ML
              if (m(18).le.1D-20) m(18)=1D-30



             xis=(x(I)**3+ mv(NP)*xn(I,np)/4/pi*3)**(1/3d0) ! radius fo shell I with solid only
c     X(I+1) is the size with solid and liquid of shell I
c     xis is used for plot program
c     
c     output molalites,aw, radius and pH etc. in each shell 
             phh=dlog(xhshell(i)* ml6shell(i))/dlog(.1d0)
              write(3,'(I6,501E14.5)') I,
     &        x(I+1),aw,(M(j),J=1,20),flu,sarsa,(M(j),J=23,NP),smisch, 
     & xhshell(i),pHH,xis

c     output the moles in each shells 
              write(23,'(I5,501E19.11)') i,time,
     &        (xn(i,j),J=1,np)
             
 
       enddo


        if (imode_output.eq.4) print*, 'phh  =' , phh

c     Titer conversion to relative concentration
      sflu= sflu/stiter0
      if (sflu.le.1D-30) sflu=1D-30
      ssars= ssars/stiter0
      if (ssars.le.1D-30) ssars=1D-30

      write(33,'(1F17.7,2E16.4,1E16.8,14F15.5)')
     &     time,sflu,ssars,xnsolid,tdrop,t
      write(32,'(1F17.7,2E16.4,1E16.8,14F15.5)')
     &     time, flux(1,np),   flux(1,NP+1)

c     output the partial pressures
          WRITE(14,'(1F17.7,10E15.6)') TIME, RH, pacetic/press*1D9
     &   ,1D9*pnh3/press,pHNO3*1D9/press,PHCL*1D9/press,PCO2/press*1D6
     &     ,tdrop

c     mass solutes without crystal
             xmsL = xms -xnsolid*Mm(NP)
          
c     xms : mass of all solutes inclusive solid
c     xmw: mass of H2O
             
          rhoo=(xms+xmw)/4/pi* 3 / x(NS+1)**3

c     output the vapour pressures
          WRITE(15,'(1F17.7,10E15.6)') TIME, aw, partvap3/press*1D9
     &         ,1D9*partvap4/press,partvapHNO3*1D9/press,
     &         PartvapHcl*1D9/press,PartvapCO2/press*1D6, xms/(xms+xmw)
     & ,      xmsl/(xmw+xmsL),awin, rhoo

                    Nout=nout+1
          
          endif
c     end output

          

c     goto time loop 11
          
          IF (TIME .LE. timetr0(ntr0) ) GOTO 11

          

 200    continue

        open(2,file='output_ML.dat')
c     output molalities of outermost shell
c     it serves for the input ML.dat

        if (imode_output.eq.4) then
        DO J=1, NP
        write(2,*) J,ml(1)*xn(ns,J)/xn(ns,1)
        enddo
            write(2,'(A,2E15.6)') 'aw = ', aw
            write(2,'(A,2E15.6)') 'pH = ', phh
            write(2,'(A,2E15.6)') 'pCO2 = ', pco2*1E3 !ppm
            write(2,'(A,2E15.6)') 'pAA = ', pacetic*1E6 !ppb
            write(2,'(A,2E15.6)') 'pNH3 = ', pNH3*1E6 !ppb
            write(2,'(A,2E15.6)') 'pHNO3 = ', pHNO3*1E6 !ppb
            write(6,'(A,2E15.6)') 'RH = ', aw
            print*, 'ML(5)/ML(16)',ML(5)/ML(16)
        endif

        close(2)

        

        
       stop
       end

      
c     ---------------------------------------------------
c     Functions and subroutines
c     ---------------------------------------------------


c     -------------------------------------------------
c
c     LOOKUP TABLE for water activity and activity coefffients
c     The model is aw_back_model
c      
c     -------------------------------------------------

      subroutine aw_back
     & (t,MLa,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

      implicit real*8 (a-h,m,o-z)
      integer NP
      parameter (np=29)
      real*8 ML(Np),mla(*),ml0(NP)


      real*8 MM(NP) ! molar mass
      real*8 Mv(NP) ! molar volume
      common /M/ MM,mv          
      common /awin/ awin,aws,xvol,xmi
      data key /-1/
      

      parameter (IT=5001,iA=101,IN=101)
       real*4 data(6,IT,IA,IN)
       save data,key

       real*8 p2(6,2,2)
       real*8 p1(6,2)
       real*8 p0(6)
       real*8 x1(1),y1(1),x2(2),y2(2)

      common/solid/xnsolid
      common /output/ imode_output,idiff,imode_pH
      
       n2=2
       n1=1
c     no look up take for imode_output 4

       xms=mla(28)+mla(27)
       if (imode_output.eq.4 .or. xms.gt.0d02) then
      call aw_back_model
     & (t,MLa,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      return
          endif
       if (key.ne.1) then
          key=1
       ML=0
       ML(1)=1000d0/MM(1)
      DO I=1,IT
         if (I/100*1d0.eq. I/100d0) print*,i

         xtotal = (I-1d0)/50d0
       if (I.le.1) xtotal =1D-5
    
      DO j=1,IA
        fa = (J-1)/100d0
         ML(16)=xtotal *(1-fa)
         ML(12)=xtotal *(fa)
         if (ML(12).ge.47) ML(12)=47

      DO K=1,IN
        fN = (K-1)/100d0
         ML(17)=xtotal*(1-fn)
         ML(18)=xtotal*fn

      call aw_back_model
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
      xmi=2*xtotal
      
        data (1, i,j,k) =aw
        data (2, i,j,k) =gammah
        data (3, i,j,k) =gammanh4
        data (4, i,j,k) =gammaNa
        data (5, i,j,k) =gammaCl
        data (6, i,j,k) =gammaNo3

        enddo
        enddo
        enddo
        endif
        
cc
 11     continue

        DO I=1,np
           ml0(I)=mla(I)
        enddo
        
        xmi = ML0(6)+ml0(12)+ml0(16)+ ML0(17)+ml0(18)+ml0(19)+ml0(20)
        xmi= xmi+ mL0(25)+mL0(26)
        xt=xmi/2


        
        if (xt.lt.0d0) xt=0d0
        if (xt.ge.99.9) xt=99.9
        
        i1= xt*50d0 + 1
        i2= I1 + 1

        if (I2.ge.it) then
           i2=it
           i1=it-1
        endif
        
        
        j1= ML0(12)/xt*100 + 1
        k1= ML0(18)/xt*100 + 1

        if (j1.le.1) J1=1
        J2=J1+1

        if (j2.gt.IA) then
           J2=IA
           J1=IA-1
        endif
        
        if (k1.le.1)  k1=1

        k2=k1+1
        if (k2.gt.IN) then
           k2=IN
           k1=In-1
        endif

          X2(1) = (I1-1d0)/50d0
          X2(2)=  (I2-1d0)/50d0
          x1(1)= xt
          
c          write(6,'(3F10.3,1E15.6)') x2(1),x1(1),x2(2)

          II1=1
          II6=6


          DO k=K1,K2
          DO j=J1,J2


             DO II=ii1,ii6

             if (I2.le.1) then
              p2(II,j-j1+1,k-k1+1) =data (II, 1,j,k)
               goto 13
               endif

             if (I1.ge.IT) then
               p2(ii,j-j1+1,k-k1+1)  = data (II, IT,j,k)
               goto 13
               endif

             y2(1)=data (II, I1,j,k)
             y2(2)=data (II, I2,j,k)
             call intpl(x2,y2,N2,x1,y1,n1)
             p2(ii,j-j1+1,k-k1+1)  = y1(1)
                
 13            continue
c               print*,j-j1+1,k-k1+1, p2(ii,j-j1+1,k-k1+1) 
               enddo

              enddo
              enddo


c     interpolate NH4

          X2(1)= xt *(j1-1)/100d0
          X2(2)= xt *(j2-1)/100d0
          x1(1)= ml0(12)
          if (x1(1).gt.xt) x1(1)=xt
          
c          print*, 'NH4 '
c          write(6,'(3F10.3,1E15.6)') x2(1),x1(1),x2(2)


          DO k=k1,k2

             DO II=ii1,ii6

             if (j2.le.1) then
               p1(II,k-k1+1) = p2(II,1,k-k1+1)
               goto 14
               endif




             y2(1)=p2(II,1,k-k1+1)
             y2(2)=p2(II,2,k-k1+1)
             call intpl(x2,y2,N2,x1,y1,n1)

               p1(ii,k-k1+1) = y1(1)

                
 14            continue
c               print*,k-k1+1, p1(ii,k-k1+1) 
               enddo
            enddo

c     interpolate NO3-
          X2(1)= xt *(k1-1)/100d0
          X2(2)= xt *(k2-1)/100d0
          x1(1)= ml0(18)
          if (x1(1).gt.xt) x1(1)=xt

c          print*, 'NO3- '
c          write(6,'(3F10.3,1E15.6)') x2(1),x1(1),x2(2)
            
             DO II=ii1,ii6

             if (k2.le.1) then
               p0(II) = p1(II,1)
               goto 15
               endif



             y2(1)=p1(II,1)
             y2(2)=p1(II,2)
             call intpl(x2,y2,N2,x1,y1,n1)
               p0(ii) = y1(1)

 15            continue

         ML(6)=ML(18)+ml(17)-ml(12)-ml(16)

c      call aw_back_model
c     & (t,ML0,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

         p0(II)=(p0(II))
         
c      write(6,'(3F10.3,11E15.6)') xt,ML0(12),ml0(18),p0(ii)
      enddo


         awin = p0(1)
         gammah=p0(2)
         gammaNH4=p0(3)
         gammaNA=p0(4)
         gammaCL=p0(5)
         gammaNO3=p0(6)

                xvH2O = 1000/MM(1)*MV(1)

        xmH2O=1000/mm(1)

        xvol = xmH2O/(XmH2O
     &       + ML0(8)+ ML0(9)+ ML0(5) + ML0(10) ) ! acetic acid ammonium acetate

        aws = xmh2o/(ML0(11)+ml0(2)+xmh2o) ! organics Raoult's law

        aw = awin*aws*xvol
c     makes the empirical correction
        if (imode_output.eq.2) then
        aww = aw
        aw1=.75
        aw2=.55
        
        if (aww.ge.aw1) aww = aw1
        if (aww.le.aw2) aww = aw2
                ff = 1 -.11 * ( aw1- aww )/ (aw1-aw2)

                 if (xnsolid.gt.1d-30) then
        aww = aw
        aw1=.75
        aw2=.55
                if (aww.ge.aw1) aww = aw1
        if (aww.le.aw2) aww = aw2
                ff = 1! -.015 * (  aww -aw2)/ (aw1-aw2)
                endif
                

        aws = aws*ff ! organics Raoult's law with correction
        aw = awin*aws*xvol
      endif
      
      return
      end

c     -----------------------------------------------------------------------------------------------------
c     This subroutine calcualtes the flux f of all species and timestep dtime including solid NaCl in all shells.
c
c     x
c     f(i,np) is the flux of the inner liquid shell to the NaCl crystal
c     f(i,np+1) the flux from the liquid to NaCl solid in the shell i
c     f(i,j): is the flux of species j from shell i to i-1.
c     f(NS+1,NP) is the flux from gas phase to the particle
     
c     output x: NS+1 array, x(i) and x(i+1) the lower and higher boundary (radius) of first shell
c
c     -----------------------------------------------------------------------------------------------------
c     the time step is dynamically calculated to ensure a stable compostion and pH
c      


      subroutine cal_flux(time,xn,x,f,dtime,NS)

       implicit real*8 (a-h,m, o-z)
       parameter (np=29,nsmm=1000)
       real*8 mm(NP),mv(NP),ml(NP)
       
       real*8 vshell(nsmm),
     & awshell(nsmm),xhshell(NSMM),ml6shell(NSMM),sNaClshell(NSMM)
       
       
      real*8 dl_factor2(2,np),dl_factor(NP)
      common /DL/ DL_factor2 ,deltaxgas,jmin,nmin
      common /timeeq/ timeeq

       

       common/flux/T,Ta,press,rh,partvap3,partvap4,partvapco2,
     & partvaphno3,
     &      partvaphcl
                    
      real*8 x(*),xn(NSMM,np),f(NSMM,np+1)
      real*8 w1(NSMM),vap(np),gammaclshell(NSMM)

      real*8 g(NSMM,np),c(NSMM,np),a2(NSMM,np)

      common /pi/ pi
      integer izc(np)

      common /M/ MM,mv,izc

      common /awshell/ awshell,xhshell,ml6shell,sNaClshell,gammaclshell

      common /solid/xnsolid
      common /output/imode_output,idiff,imode_pH
      common /enhance/radius,r1
      common /Ienhance/Ienh,iscenter
      

       n1=1
       diffmin=1D10
       
       dtime=1D10

       nmin=-1
       jmin=-1
       jgas=-1


c     calculate the thickness and concentration in each shell

       xvsolid=0d0
       vols=xvsolid ! cm3

       DO I=1,Ns
         vol= mv(1)*xn(I,1)+mv(2)*xn(I,2)

         DO J =6,20
            vol=vol+xn(I,j)*mv(j)
         enddo
c     exclude Titers
         DO j=23,np
            vol=vol+xn(I,j)*mv(j)
         enddo

         vols=vols+vol
         ff= 1000d0 /(MM(1)*xn(I,1))
         ML(1)=1000d0/MM(1)
         DO J=2,NP
            ML(J)= ff* xn(I,j)
         enddo


       	call calHNew(Ta,mL)
        call aw_back(ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

c     check aw
         sigma = 72
          if (time.ge.timeeq) then
     
c     check aw
         xkelvin = dexp( 2* sigma * 18d0 /(8.314E7*T*x(NS+1)) )
         rhh= rh/xkelvin
        if (dabs(Aw-rhh).ge..3D-3) then
           write(28,'(A,1F15.4)')'time = ', time
           print*, time,xkelvin*aw, rh
           write(28,'(A,2F10.4)')'before ',  xkelvin*aw, rh

  	   call cal_ml(T,rhh,ML)
           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

c           print*, 'after  cal_ml'
c           print*, xkelvin*aw, rh
c           write(28,'(A,1E15.6)') 'after'
           write(28,'(A,2F10.4)')'after ',  xkelvin*aw, rh



           xn(I,1) = (xn(I,2)+xn(I,16))*ML(1)/(ML(2)+ml(16))

c     repeat check
         ff= 1000d0 /(MM(1)*xn(I,1))

         DO J=2,NP
            ML(J)= ff* xn(I,j)
         enddo

           call aw_back
     & (ta,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

          print*, 'repeat check', xkelvin*aw, rh
c           write(28,'(A,2E15.6)') 'repeat check', xkelvin*aw, rh


        endif
                endif

        a2(I,1)= aw
        DO j=2,np
        a2(i,j)= ML(J)
        enddo
        xmplus= ml(6)+ml(16)+ ml(12)
        gammaOH =get_gammaOH(xMplus)

        a2(I,6)= a2(I,6)*gammaH
        a2(I,12)= a2(I,12)*gammaNH4
        a2(I,16)= a2(I,16)*gammaNa
        a2(I,17)= a2(I,17)*gammaCl
        a2(I,18)= a2(I,18)*gammaNO3
        a2(I,19)= a2(I,19)*gammaNa
        a2(I,20)= a2(I,20)*gammaCl
        a2(I,7)= a2(I,7)*gammaOH

         awshell(I)=aw
         Amisch1= ML(16)*ml(17)*gammaCl*gammaNa
         sNaClshell(I)=Amisch1
         xhshell(I)=gammah
         ml6shell(I)=ML(6)
         gammaclshell(I)        =gammacl

         if (I.eq.NS) then
         call vapnew(Ta,ML,aw,pacetic,pnh3,pHNO3,PHCL,PCO2)

        vap(1)=aw*vwater(T) ! take the air T head considered via Pruppacher

         vap(3) = pacetic
         vap(4) = pNH3
         vap(5) = pCO2
         vap(18) = pHNO3
          vap(17) = pHCL

          endif
       
         ff= xn(I,1)/ML(1)

c     recalculates equilibrium species
         DO j=6,15
            xn(I,j) = ff*ML(J)
         enddo

         DO J=23,np
            xn(I,j) = ff*ML(j)
         enddo

         vshell(I)= vol -mv(np)*xn(I,np)
c     liquid volume
         
         xx = 3d0/4d0/pi* vols
              x(I+1) = (xx)**(1/3d0)
              DO J=1,NP
              c(I,j)= xn(I,j)/vshell(I) ! mol / cm3 
              enddo              
           enddo

c     important needed for the enhancement factor 
           radius=x(NS+1)
           r1= (MV(NP)*xn(1,np)/4d0/pi*3d0)**(1/3d0)

          DTIMEMIN= 1
          xnsolid=0d0

          DO I=1,NS
              if (xn(I,np).gt.0d0 ) then
                 xnsolid=xnsolid+ xn(I,np)
                 Ienh=1
                 call cal_Misch_flux (NS, I,ta, x, xn, ff1,ff2)
                 f(I,NP+1) = ff2  !flux from outer shell to solid
                 f(I,NP) = ff1  !flux from solid to inner shell

                if (xn(I,np).ge. 1D-23) then
                   dtime1= 1D-2/(dabs( ff1+ff2)/xn(i,np))
                   if(dtime1 .le.dtime) then
                   dtime=dtime1
                   NMIN=I
                   jmin= NP
                 endif
                endif

              endif

              enddo

        Ienh=0
        DO J=NP-1,1,-1
c      set the flux at radius = 0 to zero when no solid

           f(1,J)=0
           
           if (dl_factor2(1,J).gt.0) then

c     calculates aw
c     Special treatment for solid NaCl for f(1,16) Na+ and f(1,17) Na+ and Cl-
c     f(1,16) = f(1,17) Na+ and Cl- has the same flux
c     D16 * (C1 -c0) / dx
              Ienh=0
                             
       DO I=2,NS

       thickness = (x(I+1) -  x(I-1))/2
c     calculates the gradient in concentration
          
          g(i,J)= (c(I,J)-c(I-1,J))/thickness


          cmm= (c(I,j)+c(I-1,j))



          if (a2(I,j) + a2(I-1,j) .gt.1D-30 .and. 
     & (Idiff.le.1 .or. idiff.ge.5) ) then

c     takes the activity/aw as the driving force
          g(i,J)= cmM/thickness *
     & (a2(I,J)-a2(I-1,j))/(a2(I,j)+a2(I-1,j))
          endif

c     calculates the diffusion coefficients
c     takes the mean value of shell I and I-1
       aw=(awshell(I)+awshell(I-1))/2 ! take the mean value of the neighboring shells

c     aw0=awshell(J)
         dl_factor(J)=dl_factor2(1,J)  

         ienh=0
c     enhanced diff only when solid is in the center
         if (j.eq.16 .and. xnsolid.gt.0d0 .and. iscenter .eq.1) ienh=1
         if (j.eq.17 .and. xnsolid.gt.0d0 .and. iscenter .eq.1) ienh=1
 
       call cal_dlaw(T,aw,dl)  !diffusion coefficient of ions

c      if (idiff.eq.3)   call cal_dlaw_suc(T,aw,dl) 
      if (idiff.eq.4)   call cal_dlaw_citric(T,aw,dl)
      if (idiff.eq.5)   call cal_dlaw_walker(T,aw,dl)
      if (idiff.eq.6)   call cal_dlaw_walker_mod(T,aw,dl)
       
c     neutral species take  Dl  of H2O
      if (j.eq. 9 .or. j.eq. 10 .or. j.eq. 23 .or. 
     & j .eq. 14 .or.j.eq. 13 ) then
       call cal_dlaw_walker_mod(T,aw,dl)
      if (idiff.eq.5)   call cal_dlaw_walker(T,aw,dl)
       endif

             d= dl*dl_factor(J)       ! scaling facotr for species J, 

c     for water takes modified r data 

       if (j.eq.1) then
         if (idiff.le.2)    then
        call cal_dlaw_walker_mod(T,aw,dl)  ! present work
         d=dl*dl_factor(J)
c     take the walker dîffusivity  for H2O
      if (idiff.eq.5)   call cal_dlaw_walker(T,aw,dl)
        endif
       endif

       
c     calculate the flux rate from shell i+1 to i : the direction

          f(I,j)= g(I,j) *d* 4*pi* x(I)**2 ! the missing minus sign changes the diffusion direction

c     for acetic acid take the sum of 8,9 10
                if (dl_factor2(1,8).gt.0 .and. j.eq.3) then
 		f(I,3)=f(I,8)+f(I,9)+f(I,10) !total acetic acid flux
                endif

c     takes the sum of NH4CH3COO, NH4+, NH3
                if (dl_factor2(1,12).gt.0 .and. j.eq.4) then
 		f(I,4)=f(I,10)+f(I,12) +f(I,23) ! total NH3 flux
                endif

                if (dl_factor2(1,15).gt.0 .and. j.eq.5) then
 		f(I,5)=f(I,13)+f(I,14) +f(I,15) ! total Co2 flux
                endif
                

       enddo              !I=2,Ns


c     the gas phase deposition has to be calculated extra.
c     calculates the gas phase deposition! no gas phase depletion is considered, constant gas phase is assumed here
       
      DH2O0=0.211*1013.d0/PRESS*(T/273.15d0)**1.94 ! Pruppacher + Klett

      rad= x(Ns+1)
      alpha=1d0
       velocity=(8.d0*8.314D7*T/PI/MM(1))**0.5
         DH2O=DH2O0/(1d0 +
     +        (4.d0*DH2O0/RAD/alpha/velocity))
         DGas =DH2O* dsqrt(MM(1)/MM(J))
c     H2O deposition
         xkelvin = dexp( 2* sigma *b 18d0 /(8.314E7*T*rad) )

         f(NS+1,j)=0d0

         
         if (J.eq.1) then
         xngaspartial = rh*vwater(T)*1D-4 /8.314/T ! mol/cm3air
         ph2ovap = vap(1) *xkelvin
         xngasvap = ph2ovap*1D-4 /8.314/Ta ! mol/cm3air
         delta = ( xngaspartial-xngasvap) 
cccccccccccccccccccccc

      vapheat=597.3*(273.15/t)**(0.167+3.67d-4*t) !cal/g
      thermcond0=((5.69+0.017*(t-273.15))*1.d-5) ! cal/s/cm/K
      velocity=(8.d0*8.314D7*T/PI/28.97d0)**0.5
       cair =1 /4.183 ! cal/K/g

       rhoair = press/1000  * 1E-1 /8.314d0/T*28.97d0 !  g/cm3

c 
      thermcond=thermcond0/(1.d0+4.d0*thermcond0/
     +                     (RAD*velocity*Rhoair*Cair))
      
      fheat=18d0*vapheat/ thermcond /T *
     &(vapheat*18d0/1.987d0/t-1)*
     &  vwater(t) *100/8.314/T*1E-6 ! mol/cm3
     & *dgas 

      f(NS+1,j) = 4*pi*rad* DGAS  *delta/(1+ fheat)

      Ta=t + vapheat*mm(1) * f(NS+1,1) /4/pi/rad/thermcond
ccccccccccccccccccccc

         
         
c     print*, dgas, ph2o0, rad, xngasvap,t,xkelvin
      endif

c     acetic acid deposition
         IF (J.eq.3 .and.  imode_ph.le.1) then
            
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvap3*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(3)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air
         dgastang= 94d0/760 *1013d0/press *(T/296)**1.75
         dgas=dgastang
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)
        endif

c NH3
      IF (J.eq.4 .and.  imode_ph.le.1) then
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )
         xngaspartial = partvap4*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(4)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air
         dgastang= 176d0/760 *1013d0/press *(T/296)**1.75
         dgas=dgastang
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)


      endif

       IF (J.eq.5 .and.  imode_ph.le.1) then
c     for CO2 steady state
c     not used
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

          dx=x(NS+1)-x(NS)

       xngaspartial = partvapco2*1D-4 /8.314/T ! mol/cm3air
       xngasvap = vap(5)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air
         
         dgastang= 0.16 *1013d0/press *(T/293.15)**1.75 ! engineering tool box
c         print*, j,dgas,dgastang
         dgas=dgastang
        f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)

        
      endif

       IF (J.eq.18 .and.  imode_pH.le.1) then !HNO3
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvapHNO3*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(18)*1D-4 /8.314/Ta*xkelvin ! mol/cm3air

         dgastang= 87d0/760 *1013d0/press *(T/296)**1.75
c         print*, j,dgas,dgastang
         dgas=dgastang
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)         

       endif

      
      IF (J.eq.17 .and.  imode_pH.le.1) then         !HCL
         xkelvin = dexp( 2* sigma * MV(J) /(8.314E7*T*x(NS+1)) )

         xngaspartial = partvapHcl*1D-4 /8.314/T ! mol/cm3air
         xngasvap = vap(17)*1D-4 /8.314/T*xkelvin ! mol/cm3air
          dgastang= 118d0/760 *1013d0/press *(T/296)**1.75
          dgas=dgastang
          
         f(NS+1,j)=4*pi*rad* dgas*( xngaspartial-xngasvap)
         if (xngaspartial .lt. -1D-40) f(NS+1,j)=0

      endif

      
      DO I = 1, Ns
                 deltax=1D-2

                 df=dabs(f(I+1,j)-f(I,j) )

                  if (I.eq.NS) then
c                     if (j.eq.3 .or. j.eq.4 .or. j.eq.17 .or. j.eq.18 )
c     & deltax=deltaxgas
                  endif

                  
                  
            if ( DL_factor2(1,J) .gt.0 .and. df.ge.1D-30 ) then

              isdis=0

        if (j.eq.8 .or. J.eq.10 .or. j.eq.9 .or. j.eq.12 .or. j.eq.23)
     & isdis=1

        if (I.eq.NS .and. J.eq.5) isdis=1
        if (j.eq.6 .or. J.eq.7) isdis=1
        if (j.eq.26 .or. J.eq.25) isdis=1
        if (j.eq.13 .or. J.eq.14 .or. j.eq.15)
     & isdis=1


        
                xm= 55.51/xn(I,1)*xn(I,j)
                

c     CO2 at NS shell  and disocciated species no check
c                write(30,*) J, isdis,df

c     
             if (  isdis.ne.1 .and. Idiff .gt. 0 .and. Idiff .le.4) then

                    if (xm.ge.1D-7) then
                       deltax = deltax /xm  !normalized  change step
                       if (deltax.ge.0.01) deltax=0.01 ! maximal 1% change 

c     for EDB fixed deltax
                  if (imode_pH.ge.2 .or. time.le.0) deltax=1d-3

                    cc=dabs( f(I+1,j)-f(I,j))/xn(I,j)  
                    dtime1=1D11
c     the total for phosphoric acid

                     xnp= xn(i,24)+xn(i,25)+xn(i,26)
                  if (J.eq.24 ) then
                  if ( xnp.ge.1D-30) then
                     df = f(I+1,24)-f(I,24)
                     df = df +f(I+1,25)-f(I,25)
                     df = df +f(I+1,26)-f(I,26)
                     cc=dabs(df)/xnp
                     else
                        cc=0d0
                     endif
                     endif



                 if (cc.ge.1D-30) then
                 dtime1=deltax/cc
                 endif

c     no H2O time constant for time > timeeq


                 if (j.eq.1 .and. time.ge.timeeq) then
                    f(I+1,1) =0
                    f(I,1) =0
                    dtime1=10d0
                    endif

                 

                  if (dtime1.lt.dtime) then
                   dtime=dtime1
                  jmin=J
                  nmin=I
                  endif

c     check HCl evaporation less than 10% of H+

c     calculates the balance of cations and anions

            if (j.eq.1 .and. imode_pH .le.1  .and. NS.eq.I ) then
                    xnn=xn(I,6) ! h+
                    if (xn(I,7).gt.xnn) xnn=xn(I,7)

                    xmh = 55.51 /xn(I,1) * xnn ! max(H+,OH-)

                    deltaH = 1D-3 * xmh  !* xhshell(I)
                    if (deltaH .lt. 1D-8)deltaH = 1D-8

c     considers anion flux except CO2 (equilibrium)
                    ff= f(I+1,17)-f(I,17)
                    ff=ff +f(I+1,18)-f(I,18)
                    ff=ff +f(I+1,20)-f(I,20)
                    ff=ff +f(I+1,25)-f(I,25)
                    ff=ff +2*f(I+1,26)-2*f(I,26)
                    ff=ff +f(I+1,3)-f(I,3)

c     consider cation flux except CO2 (equilibrium)
                    ffp= f(I+1,16)-f(I,16)
                    ffp=ffp +f(I+1,19)-f(I,19)
                    ffp=ffp +f(I+1,4)-f(I,4)


                    ffm= ff-ffp
                    ffm = 55.51 /xn(I,1) * dabs(ffm) ! dMH/dt 

                   dtime1=1
c     in case of buffer, allow larger H+ change
                   xmbuffer =ml(5)+ ml(24)+ml(25)+ml(26)
                   if (xmbuffer.ge.0.01) deltaH= 100* deltaH*xmbuffer

                   if (ffm.ge.1D-25) dtime1= deltaH/ffm
                   

                   if (dtime1.lt.dtime .and. imode_pH.eq.0 ) then
                  dtime=dtime1
                  jmin=1000
                  nmin=I
                  endif
                  endif


c     special treatment for solid diffusion

                                else

                  xnmin= 1D-7 * xn(I,1)/55.51*.3
c     ignores dissociated species, consider only the total amount
                  dtime1=1d10
                  IF ( DABS(DF).GE.1d-23 )dtime1=xnmin/df

                  if (dtime1.lt.dtime) then
c                   dtime=dtime1
c                   jmin=j
c                   nmin=i+100
                   endif
                        endif
                        
c     minimize deposition

         DD= DABS(F(NS+1,1))/XN(NS,1) * Dtime ! change of  H2O
         jgas=1

         if (xn(ns,3).ge.1D-18) then
         DD3= DABS(F(NS+1,3))/XN(NS,3) * Dtime

         if (dd3.ge.dd) then 
           dd=dd3
          jgas=3
          endif
                  endif

         if (xn(ns,4).ge.1D-19) then
         DD4= DABS(F(NS+1,4))/XN(NS,4) * Dtime

         if (dd4.ge.dd) then
           dd=dd4
         jgas=4
         endif
                  endif

         if (xn(ns,18).ge.1D-20) then
         DD4= DABS(F(NS+1,18))/XN(NS,18) * Dtime
         if (dd4.ge.dd) then
           dd=dd4
           jgas=18
           endif
          endif


         if (xn(ns,17).ge.1D-20) then
         DD4= DABS(F(NS+1,17))/XN(NS,17) * Dtime
         if (dd4.ge.dd) then
          dd=dd4
            jgas=17
         endif
                  endif


                  gasfac=deltax/dd
c                  print*,deltax


                     dtime0=dtime

                  if (dd.gt.deltax) then
                     dtime=dtime*deltax/dd
c                    gasfac= deltax/dd
                   endif


                     endif
             endif

            ENDDO

         endif
         
      ENDDO !J 


c     flux adjustment for charge balance of liquid phase diffusion
c     f+ + c*n+ = f- - c*n-
c     takes only the species with DL_factor > 0 
c     Idiff =0 or Idiff=5 and Idiff =6

      if (Idiff.gt.0 .and. Idiff .le.4 ) goto 234




      DO I=2,Ns

         ff=0
         czz=0
 
        DO J=1,20
        if (izc(j).ne.0) then
          ff=ff+ izc(J)* f(i,j)
          czz=czz+(c(i-1,j)+c(i,j))/2d0* izc(J)*izc(J)
          endif
         enddo

        DO J=23,NP
             if (izc(j).ne.0) then
          ff=ff+ izc(J)* f(i,j)
          czz=czz+(c(i-1,j)+c(i,j))/2d0* izc(J)*izc(J)
          endif

         enddo

c         write(6,'(A,13E15.6)')'pre ',dtime, fplus,fminus,fplus-fminus
         vel = -ff/czz
         
c     reduces the residue

         ff=0d0
         Vel12 = 0d0! vel

C     c         print*,vel
         DO J=1,20
          if (izc(j).ne.0) then
            f(i,j)=f(i,j) + izc(J)*vel*(c(i-1,j)+c(i,j))/2d0
            endif

            enddo

         DO J=23,NP
          if (izc(j).ne.0) then
            f(i,j)=f(i,j) + izc(J)*vel*(c(i-1,j)+c(i,j))/2d0
            endif
            enddo



         




c     for acetic acid take the sum of 8,9 10
                if (dl_factor2(1,8).gt.0d0 ) then
 		f(I,3)=f(I,8)+f(I,9)+f(I,10) !total acetic acid flux
                endif

c     takes the sum of NH4CH3COO, NH4+, NH3
                if (dl_factor2(1,12).gt.0d0 ) then
 		f(I,4)=f(I,10)+f(I,12) +f(I,23) ! total NH3 flux
                endif

c     takes the sum of CO2, H2CO3, HCO3-
                if (dl_factor2(1,15).gt.0d0 ) then
 		f(I,5)=f(I,13)+f(I,14) +f(I,15) ! total CO2 flux
                endif

      enddo


c     finish charge balance
c     recalculates dtime
      dtime=10
      
      DO j=1,NP
      
            DO I = 1, Ns
                 deltax=1D-2

                 df=dabs(f(I+1,j)-f(I,j) )

                  
            if ( DL_factor2(1,J) .gt.0 .and. df.ge.1D-40 ) then

              isdis=0

        if (j.eq.8 .or. J.eq.10 .or. j.eq.9 .or. j.eq.12 .or. j.eq.23)
     &             isdis=1
        
        if (j.eq.13 .or. J.eq.14 .or. j.eq.15)
     & isdis=1
        if (j.eq.6 .or. J.eq.7) isdis=1
        if (I.eq.NS .and. J.eq.5) isdis=1
        if (j.eq.26 .or. J.eq.25) isdis=1

                        xm= 55.51/xn(I,1)*xn(I,j)
                
c     CO2 at NS shell and disocciated species no check
c                write(30,*) J, isdis,df

                if (  isdis.ne.1) then

                    if (xm.ge.1D-7) then
                       deltax = deltax /xm  !normalized  change step
                       if (deltax.ge.0.01) deltax=0.01 ! maximal 1% change 

c     for EDB fixed deltax
                  if (imode_pH .eq. 2 .or. time.le.0) deltax=1d-3

                    cc=dabs( f(I+1,j)-f(I,j))/xn(I,j)  
                    dtime1=1D11

c     the total for phosphoric acid
                     xnp= xn(i,24)+xn(i,25)+xn(i,26)
                  if (J.eq.24) then
                  if ( xnp.ge.1D-30) then
                     df = f(I+1,24)-f(I,24)
                     df = df +f(I+1,25)-f(I,25)
                     df = df +f(I+1,26)-f(I,26)
                     cc=dabs(df)/xnp
                     else
                        cc=0d0
                     endif
                     endif


                 if (cc.ge.1D-30)dtime1=deltax/cc
                 if (j.eq.1 .and. time.ge.timeeq) then
                    f(I+1,1) =0
                    f(I,1) =0
                    dtime1=10d0
                    endif



                  if (dtime1.lt.dtime) then
                   dtime=dtime1
                  jmin=J
                  nmin=I
                  endif


c     calculates the gas phase

         if (j.eq.1 .and. imode_pH .le.1 .and. time.ge.0d0
     * .and. I.eq.NS) then

c
c                    xnn= xn(I,17)+xn(I,18)+xn(I,20)+xn(I,3) !concentration of anions
c                    xnn=xnn-xn(I,16)-xn(I,19)-xn(I,4)
                    xnn=xn(I,6)!* xhshell(I) ! h+
                    if (xn(I,7).gt.xnn) xnn=xn(I,7)

                    xmh = 55.51 /xn(I,1) * xnn ! max(H+,OH-)

                       deltaH = 1D-2 * xmh !* xhshell(I)

                    if (deltaH .lt. 1D-8)deltaH = 1D-8

c     considers the pH of outer most shell

c     anions without CO2
                    ff= f(I+1,17)-f(I,17)         !CL-
                    ff=ff +f(I+1,18)-f(I,18)      !NO3-
                    ff=ff +f(I+1,20)-f(I,20)      !Y-
                    ff=ff +f(I+1,25)-f(I,25)      !H2PO4-
                    ff=ff +2*f(I+1,26)-2*f(I,26)     !HPO4-2
                    ff=ff +(f(I+1,3)-f(I,3)) !total acetic acid
c     cations without H+
                    ffp= f(I+1,16)-f(I,16)         !Na+
                    ffp=ffp +f(I+1,19)-f(I,19)      !x+
                    ffp=ffp +f(I+1,4)-f(I,4) !total ammonium


c        net flux
                    ff=ff-ffp

                    ffm = 55.51 /xn(I,1) * dabs(ff) ! dMH/dt 
 

c     in case of buffer, allow larger H+ change
                   xmbuffer =ml(5)+ ml(24)+ml(25)+ml(26)
                   if (xmbuffer.ge.0.01) deltaH= 100* deltaH*xmbuffer

                    dtime1=1
                  if (ffm.ge.1D-25) dtime1= deltaH/ffm
c     xm


                   if (dtime1.lt.dtime .and. imode_pH.eq.0 ) then
                    dtime=dtime1
                  jmin=1001
                  nmin=I
                  endif
                  endif




                                else

                  xnmin= 1D-7 * xn(I,1)/55.51*.3
c     ignores dissociated species, consider only the total
                  dtime1=1d10
                  IF ( DABS(DF).GE.1d-40 )dtime1=xnmin/df
                  
c     
                  if (dtime1.lt.dtime) then
c                   dtime=dtime1
c                   jmin=j
c                   nmin=i+100
                   endif
                        endif
                        
c     minimize deposition




              endif
             endif

            ENDDO
            ENDDO

c     takes the gas phase flux 
            if (time.le.timeeq) then
         DD= DABS(F(NS+1,1))/XN(NS,1) * Dtime ! change of  H2O
         jgas=1
         else
            jgas=-1
            dd=-1
            endif
            
         if (xn(ns,3).ge.1D-18) then
         DD3= DABS(F(NS+1,3))/XN(NS,3) * Dtime

         if (dd3.ge.dd) then 
           dd=dd3
         jgas=3
          endif
                  endif

         if (xn(ns,4).ge.1D-19) then
         DD4= DABS(F(NS+1,4))/XN(NS,4) * Dtime


         if (dd4.ge.dd) then
           dd=dd4
         jgas=4
         endif
                  endif

         if (xn(ns,18).ge.1D-20) then
         DD4= DABS(F(NS+1,18))/XN(NS,18) * Dtime
         if (dd4.ge.dd) then
           dd=dd4
         jgas=18
           endif
          endif


         if (xn(ns,17).ge.1D-20) then
         DD4= DABS(F(NS+1,17))/XN(NS,17) * Dtime
         if (dd4.ge.dd) then
          dd=dd4
         jgas=17
         endif
         endif

         if (Jgas.eq. 4 .or. jgas.eq.3) then
         deltax = .01d0
                  if ( xn(NS,5) .ge. 5*xn(NS,jgas)) deltax=.1d0
                  if ( xn(NS,5) .ge. 100*xn(NS,jgas)) deltax=dd
            endif
c         print*, deltax,jgas


                     dtime0=dtime
                        if (dd.gt.deltax) then
                            dtime=dtime*deltax/dd
                     jmin=Jgas
                     Nmin=NS+1
                    endif
                   

      
 234  continue
      

c      xm6= ML(1) * xn(ns,6)/xn(ns,1)
c      write(50,'(2I5,14E15.6)') Jmin,nmin,time, dtime, xm6,x(NS+1)

      if (Jmin.le.0)  dtime=5d0
      if (dtime.ge.5d0) dtime=5d0

         return
       end

c     ------------------------------------------------------------------------------------
c
c     calculates the initial number of moles of each species in all shells xn at given molalities and size r0
c     each shell has the same volume 
c
c     ------------------------------------------------------------------------------------

      subroutine set_shells_vol(t,RH,r0,ML,xn,x,NS,NSmax)

       implicit real*8 (a-h,k,m, o-z)
       parameter (np=29,nsmm=1000)
       real*8 MM(NP),MV(NP),dl_factor(2,NP),ntiter(NP)
       common /DL/ DL_factor2 ,deltaxgas
       
       common /M/ MM,mv          
       real*8 x(*),xn(NSMM,np),ml(*)

         xnn=(10D-7)**3/(0.2D-7)**3
         XNTITER0=xnn/6.023D23


 
         
       pi=dacos(-1d0)

        xv=ML(1)*MV(1)+ML(2)*Mv(2)
       DO I=6, 20
          xv=xv+ ml(I)*mv(I)
       enddo
          xv=xv+ ml(NP)*mv(NP)
       xvmol =XV
       print*, xvmol

       

       
c     divide into shells with distances of about  0.3 - 2  nm.
c
c     calculates the diffusion length
         
       

c       dx = r0/Ns
       call cal_dlaw(T,rh,dl)
       dmin= dsqrt(dl*2)


       v11 = ( r0**3- (r0-dmin)**3)
       
       NS = ( r0**3)/ V11
        ns=nsmax


        if (NS.gt.NSmax) NS=NSmax
       if (NS.le.1) NS=1


       print*, NS

       
       v1 = 4d0*pi/3d0*r0**3/NS
       xntiter=xntiter0/ns
       
       if (ns.gt.9999) then
          print*, '    increase the dimension for the size bins '
          print*, ' or increase the minimum thickness!  '
          stop
          endif
        

c     calculates the moles of H2O and organics in each shell
          DO I=1,NS
             ntiter(I)=xntiter
             
             v= v1*I
          x(I+1)= (v/4d0/pi*3d0)**(1d0/3d0)
          DO J=1,nP
          xn(I,j)= ML(J)*v1/xvmol
c          print*,j,mm(J),mv(J)
       enddo
          v2=xn(I,1)*mv(1)+xn(I,2)*mv(2)
          DO J=6,NP
             v2=v2+ MV(J)*XN(I,j)
             
          enddo
          print*,I, v1,v2
       enddo
       print*, xn(NS,NP), ntiter(NS)

       return
      end



c     ----------------------------------------------------------------------------------------
c
c     calculates the diffusion coefficient of ions as a function of aw
c
c     ----------------------------------------------------------------------------------------

      

c     aw: according to He (is not used)
c     dl: diffusion coefficient in cm2/s


      subroutine cal_dlaw(T,aw0,dl)
      implicit real*8 (a-h,o-z)
      T0=293.15
      call cal_ions(T0,aw0,dl)

c     takes the T dependence of citric acid
      call  cal_dlaw_citric(t,aw0,dlt)
      call  cal_dlaw_citric(t0,aw0,dlt0)

      dl= dl* dlt/dlt0
      return
      end


c     ----------------------------------------------------------------------------------------
c
c     calculates the diffusion coeffcieint of H2O as a function of aw obatined from Figure 5 of Walker et al., 2021
c
c     ----------------------------------------------------------------------------------------

      subroutine cal_dlaw_walker(T,aw0,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1)
      common /enhance/radius,r1,enh_factor
      common /Ienhance/Ienh

      t0=293.15
      aw1=1
            call cal_dlaw_suc(T0,aw1,dl0)
      enhance = 1


      ff= enhance
            fac(1)=7E-9*ff
            fac(2)=8E-9*ff   !10E-9
           fac(3)= 8D-9*ff !11D-9
           fac(4)= 9E-9*ff

       fac(5)= 1.8E-8
       fac(6)= 6E-8
       fac(7)= dl0/10
       fac(8)= dl0/2
       fac(9)= dl0

      rh(1)= 0
      rh(2)=.5
       rh(3)=.65
       rh(4)=.75

       rh(5)=0.85
       rh(6)=0.90
       rh(7)=0.96
       rh(8)=0.99
       
       rh(9)=1

      fac=dlog(fac)
       N1=1
       N5=9
       x1(1)=aw0
       call intpl(rh,fac,n5,x1,y1,n1)
       dl=dexp(y1(1))

       call  cal_dlaw_citric(t,aw0,dlt)
       call  cal_dlaw_citric(t0,aw0,dlt0)
       dl= dl* dlt/dlt0
       
      

       return
      end



c     ----------------------------------------------------------------------------------------
c
c     calculates the diffusion coeffcieint of H2O as a function of aw obatined from from present EDB study used in Luo et al., 2022
c
c     ----------------------------------------------------------------------------------------
   
      subroutine cal_dlaw_walker_mod(T,aw0,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1)

      aw1=1

      t0=298.15
      call cal_dlaw_citric(T0,aw1,dl298)
      t293=293.15
      call cal_dlaw_citric(T,aw1,dl293)
c      print*, dl0
      
      Dl0= dl293 * 2.44D-5/dl298  ! scale d0 to 2.44E-5 cm2/s at 298d0 and aw=1
           
            fac(1)=7E-9*5
            fac(2)=8E-9*5   !10E-9
           fac(3)= 8D-9*5 !11D-9

           fac(4)= 9E-9*5

       fac(5)= 1.8E-8*5
       fac(6)= 6E-8*5
       fac(7)= dl0/4
       fac(8)= dl0*.8
       fac(9)= dl0

      rh(1)= 0
      rh(2)=.5
       rh(3)=.65
       rh(4)=.75

       rh(5)=0.85
       rh(6)=0.90
       rh(7)=0.96
       rh(8)=0.99
       
       rh(9)=1

      fac=dlog(fac)
       N1=1
       N5=9
       x1(1)=aw0
       call intpl(rh,fac,n5,x1,y1,n1)
       dl=dexp(y1(1))

       t0=293.15
       call  cal_dlaw_citric(t,aw0,dlt)
       call  cal_dlaw_citric(t0,aw0,dlt0)
       dl= dl* dlt/dlt0
       
      

       return
      end

c     ----------------------------------------------------------------------------------------
c
c     calculates the diffusion coeffcieint of H2O of citric acid Lienhardt et al., 2014
c
c     ----------------------------------------------------------------------------------------
      
      subroutine cal_dlaw_citric(t,aw,dl)
      implicit real*8 (a-h,o-z)

       a1=0.61477E+00
       a2=0.42200E+00
       a3=0.90000E-02
       b1=0.13800E+02
       b2=0.15693E+00
       b3= -0.90000E-02

c    7    -0.78000E+01
c    8   -0.20000E+02


         DH2O_1=10**(-6.514-387.4/(T-118.0))
         DH2O_0=10**(-15.0-175.0/(T-208.0))
         tc=t-273.15
         if (tc.ge.-7.8) tc=-7.8
         A=a1+a2*Tc+a3*Tc**2
         if (tc.ge.-20) tc=-20
         b=b1+b2*Tc+b3*Tc**2
         alpha=dexp((1-aw)**2*(A+aw*B))
!         print*, a,b
         dl=DH2O_1**(alpha*aw)*(DH2O_0**(1-aw*alpha)) ! m2/s
         dl=dl*1D4 ! cm2/s
         return
         end

      
c     ----------------------------------------------------------------------------------------
c
c     calculates the diffusion coeffcieint of H2O in sucrose Zobrist et al., 2011
c
c     ----------------------------------------------------------------------------------------

      subroutine cal_dlaw_suc(T,aw,dl)
         implicit real*8 (a-h,o-z)

       real*8 x(9)


      data x / 0.175,  -46.46,1.7,262.867, 10.53,-0.3,127.9,
     & 0.4514,-0.5/

c       omega1= 1-w

        

c      DO Aw= 0,1.01,.2
      x2= 1-aw
      if (x2.lt.0d0) x2=0d0


       if (x2.lt.0) x2=0

        a =  x(1)*(1+x(2)*x2)
        b =  x(4)*(1+x(5)*x2+x(6)*x2**2)
         t0 = x(7)*(1+x(8)*x2+x(9)*x2**x(3))
c         print*, aw,x2, t0
c         enddo





        dmin=1d-30

        if (t.le. t0+0.1d0) then
           dl=dmin
           else

       xx= a+ b/(T-t0)

       fcal = 10d0**(-xx)

       fcal=fcal*1d-7 ! m2/s

       dl =fcal
       endif
       if (dl.lt.dmin) dl=dmin
c      print*, aw, t0, dl
       dl=dl*1D4 !cm2/s
c      enddo

c
       return
       end




c     ----------------------------------------------------------------------------------------
c
c     calculates aw of sucrose: Zobrist et al. 2011
c     
c     ----------------------------------------------------------------------------------------

       subroutine calaw_beni(T,omega1,aw)
       implicit real*8 (a-h,o-z)
       data a /-1/
       data b /-0.99721/
       data c /0.13599/
       data d /0.001688/
       data e /-0.005151/
       data f /0.009607/
       data g /-0.006142/
 

       w2=1-omega1
       T0 =298.15d0
       aw = (1+a*w2)/(1+b*w2+c*w2**2)
       aw =aw+ (T-t0)*(d*w2+ e*w2**2+f*w2**3+g*w2**4)

      
       return
       end

c     ------------------------------------------------------------------------
c   henry's law coefficient of 
c      H+(aq) + NH3(gas)    --> NH4+(aq)  
c      k =( a_H+ * pNH3 ) / a_NH4+
c     in bar-1 
c     Rennard 2004
      
      function xknh3(T)
      implicit real*8 (a-z)

      xkw =  -0.61205E+01+ 0.44820E+04/t+0.17055E-01*t
      xkw=10d0**(-xkw)
      xkb = dexp( 16.9732 - 4411.025/T -0.044*t)
      xH =  dexp( -8.09694 + 3917.507/T -0.00314*t)
c      print*, xkw, xkb, xh
c      Renard 2004


      xknh3 = Xh*xkb/xkw
      return
      end

c     ------------------------------------------------------------------------
c   dissociationw coefficient of 
c      H+(aq) + NH3(aq)    --> NH4+(aq)  
c      k =( a_H+ * aNH3 ) / a_NH4+
c     in M
c     Renard 2004

      
            function xknh4b(T)
      implicit real*8 (a-z)

      xkw =  -0.61205E+01+ 0.44820E+04/t+0.17055E-01*t
      xkw=10d0**(-xkw)
      xkb = dexp( 16.9732 - 4411.025/T -0.044*t)
c      xH =  dexp( -8.09694 + 3917.507/T -0.00314*t)
c      print*, xkw, xkb, xh
c      Renard 2004


      xknh4b = xkb/xkw
      
      return
      end


c     --------------------------------------------
c
c     required by calHnew , calculated the H+ and OH- concentration to ensure the charge balance
c     
c     ------------------------------------------
      function fcnnew(x)

        IMPLICIT REAL*8 (A-H,O-Z)
       integer NP

       parameter (np=29)
       real*8 ML(NP), NL(NP)
       real*8 mm(NP),mv(NP)
        integer izc(NP)
       common /M/ mm,mv,izc
       common/suls/ T, ML, ppartco2
       common/gammas/ gammas1,gammas2
       
       common /isco2/iseqco2

       HCO2 = 0.034 * dexp(2300 *(1/T-1/298.15))

      
c     new ML(23): molecular nh3: RG bates Pinching 1949
c     KNH4 (NH4+ + H2O - NH3 + H3O+
c     

c     NH3 + H+ --> NH4

c      xKNH4bb = 10d0**(9.4d0)

      xKNH4bb=xknh4b(T)



      ML(6)=x
      
        ml(6)=x
        xsur=ml(2)
        XNH4= ml(4) !total NH4+ + NH4CH3COO
        xace= ML(3) !total acetic acid
        xco2=ML(5)
        xh= x
        xsul=0d0 ! no sulfat
        xno3= ML(18)

        
c     takes the activity coefficients of H+, NH4+ into account
c
        gammaH=1d0
        gammaNH4=1d0
        gammaAA=1d0
        gammaHCO3=1d0
        gammaCl=1d0
        gammas1=1d0
        gammas2=1d0

        DO II=1,2
         xkw=  dexp( -0.92644d1-0.68727E+04/t) !dissociation of H2O
c         xkw=xkw/gammaH



c         xka=1.74D-5            ! dissociation constant of acetic acid !
c     Harned and Echlers 1933


       xx= -1500.65/T-6.50923 * dlog(T)/dlog(10d0)-0.0076792* T
     & +18.67257d0
         xka= 10.d0**(xx)

        xka=xka/gammaH/gammaAA  ! takes the activity coefficient of H+
c        xksalt=10d0**(-9.53d0)    !  k for NH4CH3COO

        xksalt=4.E-5    !   Jaffe 1991 H =110 KJ/mol  and boiling point 117.1
c total =1.327E-2 hP
        xksalt=3.44E-3

        

        xksalt=xksalt/gammaNH4/gammaAA

        xkNH4bb=xknh4b(T)*gammaH/gammaNH4

         av = xksalt *( 1 +1/(x*xkNH4bb))

        bb=-xace+XNH4 + av*(1+x/xka) 
 
        aminus = 1d0/2d0/(1d0+x/xka) *(-bb+
     & dsqrt( 4d0 * (1d0+x/xka)*xace*av + bb**2d0))
 
        HA= aminus*x /xka
        xsalz = xace - ha-aminus

c     dissociation NH3
        
c        XNH4plus = XNH4 - xsalz
        XNH4A = XNH4 - xsalz
        XNH4plus = XNH4A /(1 +  1/(x *xkNH4bb)) 

c     set to ML
        ML(10)=xsalz
        ML(9)= HA              !acetic acid
        ML(8) = aminus         !CH3COO-
        ML(12)= xNH4plus
        ML(23)= xnh4A-xNH4plus

        ml(6) =x
        

        call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
        aOH = aw *xkw/(x*gammaH)
        xmplus=ML(6)+ML(12)+ml(16)
        gammaAA = 1d0! gammaCl

c     takes the water activity of NaOH for OH-
        gammaOH =get_gammaOH(xMplus)
        
        ML(7)= aOH/gammaOH
        
        
c        print*,'H+= ', x
c        print*,'HA= ', HA ! HA, A-
c        print*, 'A-  OH- ', aminus, xoh
c        print*, 'Asalt = ', xsalz, xnh4plus*aminus/xksalt
c        print*, 'NH4+= ',  xnh4plus
c        print*, 'xace, xnh4= ', xace, xnh4
c        Xk3 = 1.7D-3            ! CO2 + H2O --> H2CO3 at 298.15K
c        Xk4 = 2.5D-5            ! H2CO3 -> H+ HCO3- 
        gammaHCO3= 1d0! GammaCl

c     disscociation CO2(aq) + H2O --> H+ HCO3-

        xk34 = 4.448E-7*dexp(-2133*(1/T-1/298.15)) ! https://www.sciencedirect.com/science/article/pii/S0070457108703303
        xk34=xk34/gammaH/gammaHCO3
c     disscociation HCO3- --> H+ + CO3-2

         xk2 = 1.0855E-9*dexp(3347.3*(1/t-1/298.15d0))
        xk2 = XK2/gammaH
 
        xm15 =xk34*aw/(x)
        xm14= xm15* xk2/x

        if (iseqco2.eq.0) then

        ML(13) = ML(5)/(1+xm15+xm14)   !CO2
        if (ml(13).le.0)ml(13)=0d0
        ML(15)= ML(13)*xm15! HCO3-
        ML(14) = ML(13)*xm14    !CO3-2
        else
c           print*,' ppartco2 ', ppartco2,hco2
           ML(13)=HCO2*ppartco2/1013.5d0
           ML(15) =xm15*ML(13)
           ML(14) =xm14*ML(13)
           ML(5)=ML(13)+ml(14)+ml(15) 
           endif

c        write(6,'(8E15.6)') x,ML(3),ml(4),ml(8),ml(9),ml(10),
c     &  ML(12), ML(23)

           xkP1= 6.9E-3/gammaH
           xkp2= 6.2E-8/gammaH

           xm25 = xkp1/x
           xm26 = xkp2/x * xm25

           xptot= ml(24)+ml(25)+ml(26)
           ML(24) = xptot/(1+xm25+xm26)
           mL(25) = ml(24)* xm25
           mL(26) = ml(24)* xm26

           xms = ml(27) + ml(28)
	xlk=-4.556380021818660+dh0*(1/298.15-1/T)-275./8.314*
     &  dlog(t/298.15d0)
	xks=exp(xlk)
        xks = xks /gammah/gammas2*gammas1

        ML(27) = xms /(1 + xks/ML(6))
        ML(28) = ML(27) *xks/ML(6)


           fcnNEW=0
        do j=1,20
           FCNNEW=fcnnew+ IZC(J)* ML(J)
           ENDDO

        do j=23,NP
           FCNNEW=fcnnew+ IZC(J)* ML(J)
           ENDDO
 
        fcnnew=-fcnnew


           enddo

 	return
	end

      
c     -----------------------------------------------------------------
c
c     Calculates the vapore pressures in hPa
c
c     -----------------------------------------------------------------
      
             subroutine vapnew(T0,M0,aw,pacetic,pnh3,pHNO3,PHCL,PCO2)
        IMPLICIT REAL*8 (A-H,O-Z)
        parameter (np=29)         ! number of species
        real*8 m(NP)
        real*8 m0(NP)

        common/suls/ T,m
        
        M=M0
        T=t0
       	call calHNew(T0,m0)
        M=M0
        
        call aw_back(t0,M0,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

        xh=gammaH
        
         pNH3 = m0(12)/m0(6)/ xknh3(T)*1013.5*gammaNH4/gammaH 

         Xha =4000 * dexp(6200*(1/t-1/298.15d0)) ! Henrys law acetic acid M/bar
c     Sander 2015
         pacetic = 1013.5 * M0(9)/XhA

         xk= xkhx(T)*1000.d0**2/18.**2
c     Luo 1996
         phno3=m0(6)*m0(18)/xk*1013.5*gammaH*gammaNO3
c     Sanders 2015

         XHCO2 = 0.034 * dexp(2300 *(1/T-1/298.15))
         pco2 = m0(13) * 1013.5 /XHCO2
         	xk=xhenry(T)

c     Luo 1996
         phcl=m0(6)*m0(17)/xk*gammaH*gammacl*1013d0
        
        return
        end
      


c     --------------------------------------------
c
c     Calculates the H+ and OH- concentration and dissociation  to ensure to charge balance 
c     
c     ------------------------------------------

      	subroutine calHNew(Ta,ma)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=29)
        real*8 fcnnew
	external  fcnnew
        real*8 MA(NP),m(NP)
	common/suls/ T, M

       common /isco2/iseqco2

        common /kout/xx2

      common /output/imode_output,idiff,imode_pH

        
       iseqco2=0

	t=ta
        M= MA
        

c     H+ concentrations
        if (imode_output.eq.2 ) then
        mA(6)=1D-7
        mA(7)=1D-7
        return
        endif

 	xmin=1D-15              !pH 14
        
        if(xmin.le.0) xmin=5.E-14
 	xmax=10        ! 100% dissociation of acetic acid
        
	erabs=0.d0
	errel=0d0
	ITmax=2000

	xx1= fcnnew(xmin)
	xx2= fcnnew(xmax)


c	write(6,'(4E14.6)') xmin,xmax,xx1,xx2
	if(xx1.le.0.) then
	xmax=xmin
	goto 2
	endif
	if(xx2.ge.0.) 	goto 2
	call dzbren(fcnnew, erabs,errel,xmin,xmax,ITMAX)
 2      continue
c        print*, xmax
        M(6)=xmax
        xx2= fcnnew(xmax)
        MA=M
        
        
 999            return
	end  

c     -----------------------------------------------------------------------------
c
c     calculates the H+ and OH and discociations, CO2 is in equilibrium with the gas phase
c     defined by iseqCO2=1
c     
c     -----------------------------------------------------------------------------
      subroutine calHNewco2(Ta,ma,ppartco2)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=29)
        real*8 fcnnew
	external  fcnnew
        real*8 MA(NP),m(NP)
	common/suls/ T, M,pp
        common /isco2/iseqco2

        common /kout/xx2

      common /output/imode_output,idiff,imode_pH
        
       iseqco2=1
       pp=ppartco2


	t=ta
        M= MA
        
c        print*, 'T ' ,t

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	

c     H+ concentrations
        if (imode_output.eq.2 ) then
        mA(6)=1D-7
        mA(7)=1D-7
        return
        endif

 	xmin=1D-15              !PH 14
        
        if(xmin.le.0) xmin=5.E-14
 	xmax=10        ! 100% dissociation of acetic acid
        
	erabs=0.d0
	errel=0d0
	ITmax=2000

	xx1= fcnnew(xmin)
	xx2= fcnnew(xmax)

c       print*, xx1
c        print*, xx2
        


c	write(6,'(4E14.6)') xmin,xmax,xx1,xx2
	if(xx1.le.0.) then
	xmax=xmin
	goto 2
	endif
	if(xx2.ge.0.) 	goto 2
	call dzbren(fcnnew, erabs,errel,xmin,xmax,ITMAX)
 2      continue
c        print*, xmax
        M(6)=xmax
        xx2= fcnnew(xmax)
        MA=M
        
        
 999            return
	end  



c     -------------------------------------------------------------
c
c     Calculate inactivaton times for SARS-CoV-2
c     
c     -------------------------------------------------------------


      subroutine cal_tau_sars(T,RH,pH,tau)
      implicit real*8 (a-h,o-z)
      parameter (N=3)
      real*8 pharr(N), tauarr(N), k(N),x(4),y(1),klog(N)

        x(1)=0.50000E+01 
        x(2)=0.54638E+01 
        x(3)=0.24943E+01 
        x(4)=0.21971E+01 


        phh=ph-x(4) 
       tau_sars = x(1)   + x(2)*datan(x(3)   *phh )
       tau_sars =  dexp(tau_sars)
       tau = tau_sars/dlog(100d0)

      return
            end




c     -------------------------------------------------------------
c
c     ions diffusion coefficent, obtained from EBD data, scaled to 293 K 
c     
c     -------------------------------------------------------------

      subroutine cal_ions(T,aw0,dl)
      implicit real*8 (a-h,o-z)

      real*8 rh(17),fac(17),x1(1),y1(1)
      
      common /enhance/radius,r1,enh_factor
      common /Ienhance/Ienh,iscenter



      aw1=1

      t0=298.15
      call cal_dlaw_citric(T0,aw1,dl298)
      t293=293.15
      call cal_dlaw_citric(T293,aw1,dl293)
c      print*, dl0
      
      Dl0= dl293/dl298 * 2.44E-5  ! scale d0 to 2.44E-5 cm2/s at T0 and aw=1

      enhance = 1


      if (r1.ge.1D-12 .and. Ienh.eq.1) then

         xlc = 340d-4/20d0**3*(radius*1D4)**3

         enhance = xlc/2/r1
c         print*, enhance, radius, xlc,r1

         
      if (enhance.lt.1d0) enhance=1d0
         endif


         ff=0.54*(enhance)
            ff1=1/1.46
            fac(1)=.06D-9*ff*ff1
            fac(2)=0.07E-9*ff*ff1
            fac(3)=0.15E-9*ff*ff1
            fac(4)= .4D-9*ff *.9*ff1
            fac(5)= 1.3D-9*ff*1.1*ff1
       fac(6)= 1.8E-8*ff**(1/2d0)*ff1**(0.5d0)
       fac(7)= 6E-8 
       fac(8)= dl0/10
       fac(9)= dl0/2
       fac(10)= dl0
       
      rh(1)= .0
      rh(2)=.5
      rh(3)=.575
      rh(4)=.65
      rh(5)=.72
      
       rh(6)=0.85
       rh(7)=0.90
       rh(8)=0.96
       rh(9)=0.99
       rh(10)=1

      fac=dlog(fac)
       N1=1
       N5=10
       
       x1(1)=aw0
       call intpl(rh,fac,n5,x1,y1,n1)
       dl=dexp(y1(1))

c     no reduction of diffusion coefficient for crystal in the center
       if (iscenter.eq.1) return
c     DO reductions for co-shells crystals
       aww=aw0
       if( aww.le.0.5d0) Aww=0.5
       ff=dlog(0.4d0)*(1-aww)/0.5d0
       dl = dexp(ff)* dl        ! multiply by 0.4 at aw <= 0.5
       
       
       return
       end



c     --------------------------------------------------
c
c    activity product crystal = a_Na * a_Cl
c
c     --------------------------------------------------

      
      function apNaCL(t)
      IMPLICIT REAL*8 (A-H,O-Z)
          apnacl=      0.13852E+02 -0.21891D-1 *(t-298.15)

      return
      end
      


c     --------------------------------------------------
c
c     calculates the flux f1 (from the shell to to crystal) and f2 (the flux liquid to crystal of the same shell)
c
c     NS: the total number of shells
c      II: the location of the crystal 
c     
c     --------------------------------------------------

      subroutine  cal_Misch_flux (NS, iI,Tdrop, x, xn,f1,f2)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      parameter (np=29,nsmm=1000)
      real*8  x(*),ml(NP),xa(NSMM),ml0(np),ml1(NSMM),mlm1(NSMM)
      
       
      real*8 dl_factor2(2,np),dl_factor(NP)
      common /DL/ DL_factor2

      real*8 xn(NSMM,NP),x1(NP)
      

      common /output/imode_output,idiff,imode_pH
      
     

      common /ienhance/Ienh,iscenter
      common /pi/pi
      


      real*8 MM(NP) ! molar Mass
      real*8 Mv(NP) ! molar volume
      common /M/ MM,mv          
      common /Nacl/T,ML,ML0
      external fcn_ap


      I=II

      f1=0d0
      xvsolid=mv(Np) *xn(I,np)
      vol= 4*pi/3d0*x(I)**3 + xvsolid
      xvliquid =4*pi/3d0*( x(I+1)**3-x(I)**3) -xvsolid
      r1=(vol/4d0/pi*3)**(1/3d0)
      r2=x(I+1)

      if (I.lt.NS) then
         if (r2.le. x(I+2)/2) r2 = x(I+2)/2
      endif
      
      v1= xvliquid
      
            
      ML(1) = 1000d0/mm(1)

      t=Tdrop
c     take the mean composition of shell I and I+1
      DO J=1,NP
         ML1(J)=ML(1)* xn(I,J)/xn(I,1)
                  if (I.lt.NS) then
         ML(j)= (ML(1)* xn(I,J)/xn(I,1)+ ML(1)* xn(I+1,J)/xn(I+1,1))/2
         
        else
         ML(j)= ML(1)* xn(I,J)/xn(I,1)
         endif
c         print*, i,j, ML(j) 
      enddo

      ML0=ML
      xmax=1000d0

      xmin=1D-3
       xmax =xmax*.999999


       erabs=0.
	errel=0d0
        itmax=100
c        pi=dacos(-1d0)
      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
      Ienh=0
      if (iscenter.eq.1) ienh=1
      call cal_dlaw(T,aw,dl)
c      if (idiff.eq.3)   call cal_dlaw_suc(T,aw,dl)
      if (idiff.eq.4)   call cal_dlaw_citric(T,aw,dl)
      if (idiff.eq.5)   call cal_dlaw_walker(T,aw,dl)
      if (idiff.eq.6)   call cal_dlaw_walker_mod(T,aw,dl)

       dl_factor(16)=dl_factor2(1,16)
       dl_factor(17)=dl_factor2(1,17)
       dl_factor(11)=dl_factor2(1,11)
       dlfNa=dl_factor(16)
       dlfCl=dl_factor(17)

c       print*, 'before ap ', ml0(16),ml0(17)
        call dzbrens(fcn_ap,erabs,errel,xmin,xmax,ITMAX)
c       print*, 'after ap ', ml0(16),ml0(17)
c        call calHNew(T,mL0)
       xx2=fcn_ap(xmax)

       dlna= dl* dl_factor(16)
       DlCl= dl* dl_factor(17)
       Dl11= dl* dl_factor(11)

       c1 = xn(I,1)/v1* ML(16)/ML(1)
       c0 = xn(I,1)/v1* ML0(16)/ML0(1)
       
       fNa = dlna  * (c1-c0)*4 *pi*r1*(r2+r1)/(r2-r1)
       fna=fna +  vel12*(c0)*dl_factor(16)
       
       c1 = xn(I,1)/v1* ML(17)/ML(1)
       c0 = xn(I,1)/v1* ML0(17)/ML0(1)
       fcl= dlcl  * (c1-c0)*4*pi*r1*(r2+r1)/(r2-r1)
       fcl=fcl - vel12*c0*dl_factor(17)


c       print*, fcl, fna
       ff = 1 !vliquid/VOL
       if (I.eq.1) ff=1d0
       f2=fcl*ff
c 

       call aw_back(t,ML1,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
       Amisch1= ML1(16)*ml1(17)*gammaCl*gammaNa 

       if (Amisch1.le.apnacl(T) .and. f2 .gt.0d0) f2=0d0 !  no growth when the I shell is not saturated
       if (Amisch1.ge.apnacl(T) .and. f2 .lt.0d0) f2=0d0 !  no dissolving when NaCl is super-saturated 
       
       
       if (I.ge.2) then
          
c     calculate f1
c     the flux from inner shell to crystal  !!!!!
c     c0-c1
      r1=x(I)  ! crystal inner radius
      r2=x(I-1) ! inner radius of shell I-1

      vliq= xn(I-1,1) *mv(1)+ xn(I-1,2)*mv(2)
      DO J= 6,NP-1
         vliq= vliq + xn(I-1,J)*mv(J)
      enddo
      
         
c       ff =vsolid/(vliq+vsolid)
c     take the mean compostion of the shell I and I+1
      DO J=1,NP
         I1= I-1
         ML1(j)= ML(1)* xn(I1,J)/xn(I1,1)
         ML(j)= (ML(1)* xn(I,J)/xn(I,1)+ ML(1)* xn(I1,J)/xn(I1,1))/2
      enddo

      ML0=ML



      xmax=1000d0

      xmin=1D-3
       xmax =xmax*.999999

       xx1=fcn_ap(xmin)



       erabs=0.
	errel=0d0
        itmax=100
      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa,
     & gammas1,gammas2)        

      isenh=0
      
      call cal_dlaw(T,aw,dl)

      if (idiff.eq.4)   call cal_dlaw_citric(T,aw,dl)
      if (idiff.eq.5)   call cal_dlaw_walker(T,aw,dl)
      if (idiff.eq.6)   call cal_dlaw_walker_mod(T,aw,dl)

       dlfNa=dl_factor2(1,16)
       dlfCl=dl_factor2(1,17)


        call dzbrens(fcn_ap,erabs,errel,xmin,xmax,ITMAX)
       xx2=fcn_ap(xmax)

       dlna= dl* dl_factor2(1,16)
       DlCl= dl* dl_factor2(1,17)

       c1 = xn(1,1)/vliq* ML(16)/ML(1)
       c0 = xn(1,1)/vliq* ML0(16)/ML0(1)
       
       fNa = dlna  * (c1-c0)*4 *pi*r1*(r1+r2)/(r1-r2)

       c1 = xn(1,1)/vliq* ML(17)/ML(1)
       c0 = xn(1,1)/vliq* ML0(17)/ML0(1)
       fcl= dlcl  * (c1-c0)*4*pi*r1*(r2+r1)/(r1-r2)
c     oppsite sign as F2

       f1=(fcl+fna)/2 ! *vsolid/(v1+vsolid)/2


       call aw_back(t,ML1,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
       Amisch1= ML1(16)*ml1(17)*gammaCl*gammaNa 
c     f1 and f2 has oppsote sign
       
       if (Amisch1.le.apnacl(T) .and. f1 .gt.0d0) f1=0d0 !  no growth when the I shell is not saturated
       if (Amisch1.ge.apnacl(T) .and. f1 .lt.0d0) f1=0d0 !  no desolving when super saturated

          

       endif
            
        return
        end
      


c     --------------------------------------------------
c
c     Required by cal_misch_flux, to calculate the equilibrium Na+ and Cl- concentration over the NaCl crystal
c      
c     --------------------------------------------------

      function fcn_ap(xmcl)
      IMPLICIT REAL*8 (A-H,m,O-Z)
      parameter (np=29)
      real*8 ml(NP),ml0(NP)

     

       
      real*8 dl_factor2(2,np),dl_factor(NP)
      common /DL/ DL_factor2 ,deltaxgas

      common /Nacl/T,ML,ML0


      common /output/imode_output,idiff,imode_pH
      
  


      
c      call aw_back(t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        


       ML0(17)=xmcl
       d17=dl_factor2(1,17)
       if (d17 .le.1D-30) d17=1d0
       d16=dl_factor2(1,16)
       if (d16 .le.1D-30) d16=1d0

       ML0(16)= d17/d16 * (ML0(17)-ML(17))+ ML(16)
      
       if (ML0(16).le.1D-10 )  ml0(16)=1D-10


       call aw_back(t,ML0,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)        
 
c       print*, 'ml0 ',ml0(16),ml0(17) 
       Amisch1= ML0(16)*ml0(17)*gammaCl*gammaNa 
       fcn_ap= Amisch1-Apnacl(T)
       
       return
      end
      

c     ------------------------------------------------
c
c     calculates the H2O diffusion coefficient in Pinene solution Lienhard et al., 2015
c
c     ------------------------------------------------

      subroutine       cal_dlaw_piene(T,aw,dl,xw)
      implicit real*8 (a-z)
      aw1=1d0
      call  cal_dlaw_citric(t,aw1,dw1)
      dw0 = 7E-11* dexp ( -65500/8.314*(1/T-1/300d0))
      alpha=1
      ta=T
      if (T.ge.273d0) ta =273d0
      A= -18.31+0.063*TA
      B= -10.65 + .039*TA
      alpha= (1-xw)**2 *( A +3*b-4*b*(1-xw))
c      print*,xw,dw1,alpha
      alpha=dexp(alpha)
      
      Dl = dw0** (1-xw*alpha) * dw1**(xw*alpha)
      return
      end
      

c     -------------------------------------------------------
c
c     calculates aw Pinene solution Lienhard et al., 2015
c
c     -------------------------------------------------------
      
      function       aw_pinene(xw)
      implicit real*8 (a-z)
      xmfs  = (1-xw)* 150 / (xw*18 + (1-xw)* 150)
      aw_pinene = (1-xmfs)/( 1 - 0.85848*xmfs-0.09026*xmfs**2)

      return
      end

c     ---------------------------------
c
c     calculates the OH- activity coefficient in a NaOH solution
c
c     ---------------------------------
      
      function get_gammaOH(xMplus)
      implicit real*8 (a-h,o-z)
      parameter (N=240,n1=1)
      real*8  x(N),y(N),x1(N1),y1(N1)
      data key /0/
      data x /
     &    0.0000,     0.0010,     0.0020,     0.0030,     0.0040, 
     &    0.0050,     0.0060,     0.0070,     0.0080,     0.0090, 
     &    0.0100,     0.0110,     0.0120,     0.0130,     0.0140, 
     &    0.0150,     0.0160,     0.0170,     0.0180,     0.0190, 
     &    0.0200,     0.0210,     0.0220,     0.0230,     0.0240, 
     &    0.0250,     0.0260,     0.0270,     0.0280,     0.0290, 
     &    0.0300,     0.0310,     0.0320,     0.0330,     0.0340, 
     &    0.0350,     0.0360,     0.0370,     0.0380,     0.0390, 
     &    0.0400,     0.0410,     0.0420,     0.0430,     0.0440, 
     &    0.0450,     0.0460,     0.0470,     0.0480,     0.0490, 
     &    0.0500,     0.0510,     0.0520,     0.0530,     0.0540, 
     &    0.0550,     0.0560,     0.0570,     0.0580,     0.0590, 
     &    0.0600,     0.0610,     0.0620,     0.0630,     0.0640, 
     &    0.0650,     0.0660,     0.0670,     0.0680,     0.0690, 
     &    0.0700,     0.0710,     0.0720,     0.0730,     0.0740, 
     &    0.0750,     0.0760,     0.0770,     0.0780,     0.0790, 
     &    0.0800,     0.0810,     0.0820,     0.0830,     0.0840, 
     &    0.0850,     0.0860,     0.0870,     0.0880,     0.0890, 
     &    0.0900,     0.0910,     0.0920,     0.0930,     0.0940, 
     &    0.0950,     0.0960,     0.0970,     0.0980,     0.0990, 
     &    0.1000,     0.1100,     0.1200,     0.1300,     0.1400, 
     &    0.1500,     0.1600,     0.1700,     0.1800,     0.1900, 
     &    0.2000,     0.2100,     0.2200,     0.2300,     0.2400, 
     &    0.2500,     0.2600,     0.2700,     0.2800,     0.2900, 
     &    0.3000,     0.3100,     0.3200,     0.3300,     0.3400, 
     &    0.3500,     0.3600,     0.3700,     0.3800,     0.3900, 
     &    0.4000,     0.4100,     0.4200,     0.4300,     0.4400, 
     &    0.4500,     0.4600,     0.4700,     0.4800,     0.4900, 
     &    0.5000,     0.5100,     0.5200,     0.5300,     0.5400, 
     &    0.5500,     0.5600,     0.5700,     0.5800,     0.5900, 
     &    0.6000,     0.6100,     0.6200,     0.6300,     0.6400, 
     &    0.6500,     0.6600,     0.6700,     0.6800,     0.6900, 
     &    0.7000,     0.7100,     0.7200,     0.7300,     0.7400, 
     &    0.7500,     0.7600,     0.7700,     0.7800,     0.7900, 
     &    0.8000,     0.8100,     0.8200,     0.8300,     0.8400, 
     &    0.8500,     0.8600,     0.8700,     0.8800,     0.8900, 
     &    0.9000,     0.9100,     0.9200,     0.9300,     0.9400, 
     &    0.9500,     0.9600,     0.9700,     0.9800,     0.9900, 
     &    1.0000,     1.2000,     2.2000,     3.2000,     4.2000, 
     &    5.2000,     6.2000,     7.2000,     8.2000,     9.2000, 
     &   10.2000,    11.2000,    12.2000,    13.2000,    14.2000, 
     &   15.2000,    16.2000,    17.2000,    18.2000,    19.2000, 
     &   20.2000,    21.2000,    22.2000,    23.2000,    24.2000, 
     &   25.2000,    26.2000,    27.2000,    28.2000,    29.2000, 
     &   30.2000,    31.2000,    32.2000,    33.2000,    34.2000, 
     &   35.2000,    36.2000,    37.2000,    38.2000,    39.2000, 
     &   40.2000,    41.2000,    42.2000,    43.2000,    44.2000, 
     &   45.2000,    46.2000,    47.2000,    48.2000,    49.2000/
      data y/
     &    1.0000,     0.9593,     0.9407,     0.9259,     0.9131, 
     &    0.9016,     0.8911,     0.8813,     0.8721,     0.8635, 
     &    0.8553,     0.8474,     0.8399,     0.8327,     0.8257, 
     &    0.8190,     0.8125,     0.8061,     0.8000,     0.7941, 
     &    0.7883,     0.7826,     0.7771,     0.7718,     0.7665, 
     &    0.7614,     0.7564,     0.7515,     0.7467,     0.7420, 
     &    0.7374,     0.7328,     0.7284,     0.7240,     0.7198, 
     &    0.7156,     0.7114,     0.7074,     0.7034,     0.6995, 
     &    0.6956,     0.6918,     0.6881,     0.6844,     0.6807, 
     &    0.6772,     0.6737,     0.6702,     0.6668,     0.6634, 
     &    0.6601,     0.6568,     0.6535,     0.6503,     0.6472, 
     &    0.6441,     0.6410,     0.6380,     0.6350,     0.6320, 
     &    0.6291,     0.6262,     0.6234,     0.6206,     0.6178, 
     &    0.6150,     0.6123,     0.6096,     0.6070,     0.6043, 
     &    0.6017,     0.5992,     0.5966,     0.5941,     0.5916, 
     &    0.5892,     0.5867,     0.5843,     0.5819,     0.5796, 
     &    0.5772,     0.5749,     0.5726,     0.5704,     0.5681, 
     &    0.5659,     0.5637,     0.5615,     0.5594,     0.5573, 
     &    0.5551,     0.5531,     0.5510,     0.5489,     0.5469, 
     &    0.5449,     0.5429,     0.5409,     0.5389,     0.5370, 
     &    0.5351,     0.5167,     0.4998,     0.4842,     0.4698, 
     &    0.4563,     0.4438,     0.4321,     0.4212,     0.4108, 
     &    0.4011,     0.3920,     0.3834,     0.3752,     0.3674, 
     &    0.3601,     0.3531,     0.3465,     0.3402,     0.3341, 
     &    0.3284,     0.3229,     0.3176,     0.3126,     0.3078, 
     &    0.3032,     0.2988,     0.2945,     0.2904,     0.2865, 
     &    0.2827,     0.2791,     0.2756,     0.2722,     0.2690, 
     &    0.2659,     0.2628,     0.2599,     0.2571,     0.2544, 
     &    0.2518,     0.2492,     0.2468,     0.2444,     0.2421, 
     &    0.2399,     0.2377,     0.2356,     0.2336,     0.2316, 
     &    0.2297,     0.2279,     0.2261,     0.2243,     0.2227, 
     &    0.2210,     0.2194,     0.2179,     0.2164,     0.2149, 
     &    0.2135,     0.2121,     0.2108,     0.2095,     0.2082, 
     &    0.2069,     0.2057,     0.2046,     0.2034,     0.2023, 
     &    0.2013,     0.2002,     0.1992,     0.1982,     0.1972, 
     &    0.1963,     0.1954,     0.1945,     0.1936,     0.1927, 
     &    0.1919,     0.1911,     0.1903,     0.1895,     0.1888, 
     &    0.1881,     0.1874,     0.1867,     0.1860,     0.1853, 
     &    0.1847,     0.1749,     0.1684,     0.1909,     0.2288, 
     &    0.2803,     0.3468,     0.4305,     0.5349,     0.6639, 
     &    0.8224,     1.0159,     1.2508,     1.5341,     1.8740, 
     &    2.2790,     2.7586,     3.3226,     3.9815,     4.7456, 
     &    5.6254,     6.6308,     7.7710,     9.0538,    10.4853, 
     &   12.0693,    13.8072,    15.6969,    17.7328,    19.9053, 
     &   22.2005,    24.6002,    27.0816,    29.6178,    32.1779, 
     &   34.7276,    37.2297,    39.6452,    41.9340,    44.0561, 
     &   45.9729,    47.6480,    49.0486,    50.1465,    50.9192, 
     &   51.3502,    51.4303,    51.1570,    50.5356,    49.5783/

      save x ,y

      

      x1(1)= xmplus
      call intpl(x,y,n,x1,y1,n1)
      get_gammaOH=y1(1)
      return
      end

      
c     --------------------------------------------
c      
c     calcalates inactivation time of IAV virus in seconds
c
c     --------------------------------------------
      
      function tau_ivea(ph,xm0)
      implicit real*8 (a-h,o-z)
      
      real*8 ph, tau_ivea, tau_ivea_1,xx,dd,xm0,xm,x(16), phh
      integer i,ii,key
      data key /0/
      save key,x
       

      if(key.eq.0) then
         key=1
        x(1) =   0.50131E+01 
           x(2)=    0.62214E+01   
          x( 3)=         0.12972E+01   
          x(4)=    0.53802E+01   
          x(5)=    0.14134E+01   
          x(6)=   -0.47140E+00   
          x(7)=   -0.20450E+01   
          x(8)=    0.46916E+00   

         DO I=1,8
c          read(25,*) ii, x(i)
       enddo
       close(25)
      endif
      
       xm=xm0-0.027
       if (xm.le.0d0) xm=0d0

         phh=ph-x(4) + xm*x(5)
         if (phh.ge.7.5d0) phh= 7.5
         tau_ivea = x(1)   + x(2)*datan(x(3)*(1+xm*x(6))   *phh )

         taul = x(7)+ x(8)* ph
         tau_ivea = ( dexp(tau_ivea) + dexp(taul))

       
       
             return
      end 


c     ---------------------------------------------------------
c
c     calculate the molality of the solution at given aw (in equlibrium) and T
c
c     ---------------------------------------------------------

      subroutine cal_ml(T,rh,ML)

         
	 implicit real*8 (a-h, o-z)
	 external  fml

	parameter (NP=29)
	real*8 mm(NP),mv(NP),mL(NP),mlA(NP),ml0(NP) ! mole mass
	common /m/ mm, mv
	common /mLA/ rh1,t1, xmfs0, ml0,MLA
        

         
         ML0=ML
c         pwrint*, 'x', x, MV(1)

         xms=MM(NP)*ML(NP)+MM(2)*ML(2)
         DO I=6, np
            xms=xms+ MM(I)*ML(I)
            enddo
            xmfs0=xms/1000


	 t1=t
	 rh1=rh

c     x = solute /H2O ratio
	 xmin=.5
	 xmax=1.1
         

	 erabs=0.d0
	 errel=0d0
	 ITmax=1000

         xmin = 1D-10
 	 ymin=fml(xmin)
	 ymax=fml(xmax)
c	 ymax1=fml(15d0)
 

         if (ymin.ge.0) goto 12
         if (ymax.le.0) goto 12
	 call dzbrens(fml,erabs,errel,xmin,xmax,ITMAX)
	 ymax=fml(xmax)

  12	 continue
         ML=MLA

      call aw_back
     & (t,MLA,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
 
  
	 return

	 end

c     ---------------------------------------------------------
c
c     required by cal_ml, to calculate the molalities at given aw
c      
c     ---------------------------------------------------------
      
	 function fml(X) 
	 implicit real*8 (a-h,o-z)
	parameter (NP=29)
	 real*8 mm(NP),mv(NP),mL0(NP),mla(NP) ! mole mass

	common /m/ mm, mv
	common /mLA/ rh,t,xmfs0, ml0,MLA

c        print*, x,xmfs0
	
        MLA(1)= 1000d0/mm(1)

           DO I=2,NP
            MLA(I)= ML0(I) *x
           enddo

      call aw_back
     & (t,MLA,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)

 	   Fml=(RH-aw)

c           print*, aw



 10 	return
	 end
 
c     ---------------------------------------------------------------------------------------
c      
c     calculates the water activity and acitivity coeffcients of ions using the Pitzer ion interaction model for ions
c      
c     for neutral species, The Raoult's law is used
c
c     ---------------------------------------------------------------------------------------
      subroutine aw_back_model
     & (t,ML,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
        implicit real*8 (a-h,m,o-z)
      integer NP
      parameter (np=29)
      real*8 ML(Np)

      real*8 MM(NP) ! molar mass
      real*8 Mv(NP) ! molar volume
      common /M/ MM,mv          
      common /awin/ awin     ,aws   ,xvol
      
	Parameter( NMAx=5 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax),C0(nmax,nmax)
 	real*8 C1(nmax,nmax),omega(nmax,nmax),mc0(nmax),ma0(nmax)
	real*8 MC(nmax),MA(nmax),ZC(nmax),ZA(nmax),xs(100),xfit(10)
        
	integer Iflag(NMAX,NMAX)
	data ZC /5*1./,ZA/5*1./
        common /a/ ah2so4

      common /output/imode_output,idiff,imode_pH
      
      
        
       common/gammas/ gammas1,gammas2
        data xfit /
     1    0.65165E-02, 
     2    0.92664E-01, 
     3    0.19114E-04, 
     4    0.78199E-01 ,
     5    0.12843E+01 ,
     6    0.79342E-02, 
     7   -0.72711E-01, 
     8    0.69578E-04 ,
     9    0.22108E-01 ,
     1    0.50531E+00 /

        za(4)=2

        iflag=0 
      	Iflag(1,1)=1 !H+ NO3-
	Iflag(1,2)=2 !H CL

	Iflag(2,1)=7            !NH4 NO3
	Iflag(2,2)=8  ! NH4 Cl
	NA=4
	NC=2

	Iflag(2,3)=5 ! NH4-HSO4
	Iflag(2,4)=6 ! NH4-SO4
	Iflag(1,3)=3 !H-HSO4
	Iflag(1,4)=4  ! H- SO4
        
        
	call calpar(T,NC,NA,b0,b1,C0,C1,omega,xs,Iflag)
	NC=3

c NaNO3

        b0(NC,1)= -0.34873E-02  
        b1(NC,1)=0
        c0(NC,1)= 0.31888E-05   
        c1(NC,1)=0.66301E-02   
        omega(NC,1)= 0.84726E+00   


        
c     replace NH4Cl
        b0(2,2)=0.29535E-02
        b1(2,2)=0.31808E+00  
        c0(2,2)=-0.14979E-04  
        c1(2,2)=0.49032E-01  
        omega(2,2)= 0.12116E+01 


        
c     NaCl
        b0(NC,2)=-0.11243d0
        b1(NC,2)=0
        c0(NC,2)= 0.92338E-03
        c1(NC,2)= 0.12494E+00
        omega(NC,2)= 0.84788E+00

c     with Ulrich Kriegers Data
        b0(NC,2)=-0.16602E-01
        b1(NC,2)= -1d0
        c0(NC,2)= 0.11126E-03
        c1(NC,2)=  0.10771
        omega(NC,2)=  0.9829d0
        

        b0(NC,3)=xfit(1)
        b1(NC,3)=xfit(2)
        c0(NC,3)=xfit(3)
        c1(NC,3)=xfit(4)
        omega(NC,3)=xfit(5)


c     NA - SO4
        b0(NC,4)=xfit(6)
        b1(NC,4)=xfit(7)
        c0(NC,4)=xfit(8)
        c1(NC,4)=xfit(9)
        omega(NC,4)=xfit(10)



        
        
        MC(1)= ml(6)
        MC(2)= ml(12)
c     treat  ML(19) as Na+  halfing
        MC(NC)= ml(16) +ML(19)*1.0
        
        mc0=mc
c     set the Na molality to maximal 40

        Mmax=100d0
        if (MC(2).ge.mmax) MC(2)=Mmax
        if (MC(NC).ge.mmax) MC(NC)=Mmax

        Ma(1)= ml(18) +ML(25)+ml(26) !  treat phosphoric acid as nitric acid
c     Treat ML(20) as Cl- halfling
        Ma(2)= ml(17)+ML(20)*1d0 

        ma0=ma
        
        if (Ma(2).ge.Mmax) Ma(2)=Mmax
        if (Ma(1).ge.Mmax) Ma(1)=Mmax

c     set the Cl molality to maximal 40
        ma(3)=ML(25)
        ma(3)=ML(26)
        
        
      	awin=gammasn(T,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega) ! Pitzer model ions

        xmm=0
        xmm0=0
        DO I=1,2
        xmm=xmm+ma(I)
        xmm0=xmm0+ma0(I)
        enddo

        DO I=1,Nc
        xmm=xmm+mc(I)
        xmm0=xmm0+mc0(I)
        enddo


        awin=awin* (1000d0/MM(1)+ xmm)/(1000d0/MM(1)+xmm0)

        
c     calculates the aw of acetic acid, NH4CH4COO, and CO2

        xvH2O = 1000/MM(1)*MV(1)

c     vmol fraction of acetic acid and Co2 in ammonium acetate
c     and Phosphoric acid
        xmH2O=1000/mm(1)

c     all inorganic species not treated in Pitzer model

        
        xvol = xmH2O/(XmH2O
     &       + ML(8)+ ML(9)+ ML(5) + ML(10) +ml (23) ) ! acetic acid ammonium acetate

c     all inorganic species not treated in Pitzer model
        aws = xmh2o/(ML(11)+ml(2)+xmh2o) ! organics Raoult’s law



c     do empircial correction only for EDB
c        if (imode_output.ne.0) then
c        awin= aw_corr(awin)
c        endif

        aw=awin*aws*xvol
c        print*, ' Aaws-model ', awin, aws, xwol
        
c     print*, awin,aws,xvol
c        if  (imode_output.eq.0) then !ony valid for labo experiments
        N1=1
	gammah=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        N1=2
	gammaNH4=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        n2=NC+1
        gammaNO3=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

        n2=NC+3
        gammas1=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

        n2=NC+4
        gammas4=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        
c     else
c        gammah=1
c        gammaNO3=1
c        gammaNH4=1
c        endif


        n2=NC+2
	gammacl=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

        N1=NC
	gammaNa=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
        ap = gammaNA*gammaCl*ml(16)*ml(17)
        snacl= ap/apnacl(T)

          
        return
        end

c     ---------------------------------------------
c
c     The Henrys law coefficent of HCl in M²/atm, required by vapnew
c
c     HCl (gas )  --> H+(aq) + Cl-(aq)
c      
c     ---------------------------------------------

      function xhenry(T)
	
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 xhenry
	real*8 K0,K
	k0=2.05E6

	R=8.314
	T0=298.15	
	dh=-74852.
	dcp=-165.52
	k=log(k0)+(dh-dcp*t0)/R*(1/T0-1/T)+ dcp/R*log(T/T0)
	k=exp(K)
	xhenry=k	
	return
	end

c     ---------------------------------------------
c
c     The Henrys law coefficent of HNO3 in 1/atm, required by vapnew
c
c     HNO3 (gas )  --> H+(aq) + NO3-(aq)
c      
c     ---------------------------------------------
      
      
	function xkhx(T)

      IMPLICIT REAL*8(A-H,O-Z)
	real*8 xkhx
	x=385.9722-3020.3522/T-71.002*log(T)+.131442311*T
     * -.420928363E-4*T**2
	xkhx=exp(x)
	return
	end

c     -----------------------------------------------
c
c     calcualte the activity of ions N (1 <= N <= Nc+Na) using
c     Pitzer ion interactions model. 
c     
c     -----------------------------------------------
        function gammann(t,N,NC,NA,mC,mA,zC,zA,b0,b1,C0,C1,omega)
C	N: the index of the ion for which the activity is to be calculated
C	IF ( for cations: N <= Nc; else for anions N-NC )
C	zC, ZA: the ionic charge of type C or A.
C	b0,b1,c: Pitzer  coefficients B0(Nc, Na)
C	Nc : type of cations
C	Na : type of anions

C
C
C	This only a simple version, which only take the unsymmetrical 
C	factor e-theta and E-theta ' into account.
C	
C

	implicit real*8 (a-h,o-z)
	real*8 gammann
	Parameter( NMAx=5 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax)
 	real*8 b(nmax,nmax),Bs(nmax,nmax)
 	real*8 c0(nmax,nmax),C1(nmax,nmax)
 	real*8 c(nmax,nmax),Cs(nmax,nmax)
	real*8 MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 omega(nmax,nmax)
	real*8 I, I2

C	Calculate I
	alpha=2.
	 I=0
	DO IC=1,NC
	I=I+mC(IC)*zC(IC)**2
	enddo
	DO IA=1,NA
	I=I+mA(IA)*zA(IA)**2
	enddo
	I=I/2.
	I2=sqrt(I)

	z=0	
        gam=0.


	DO IC=1,NC
	Z=Z+ZC(IC)*MC(IC)
	enddo
	DO IA=1,NA
	Z=Z+ZA(IA)*MA(IA)
	enddo



CCC	Calculate B,Bs, C, Cs
	x=sqrt(I)*alpha
	gg=g(x)
	ggs=gs(X)

	DO Ic = 1, NC
	DO Ia = 1, Na
	B(IC,IA) = b0(Ic,Ia) + gg* b1(Ic,Ia)
	BS(IC,IA) = ggs*b1(Ic,Ia)/I
	omega1=omega(Ic,IA)
        xo=omega1*sqrt(I)
	xhx=1./xo**4*(6.-exp(-xo)*(6+6*xo+3*xo**2+xo**3))
	xhxs=exp(-xo)/2-2*xhx

	C(Ic,Ia)=(C0(Ic,Ia)+4*C1(Ic,Ia)*xhx)
	Cs(Ic,Ia)= C1(Ic,Ia)/I*xhxs
	enddo
	enddo

	Aphi=.377+4.684E-4*(T-273.15)+3.74e-6*(T-273.15)**2

	F1= -Aphi*(I2/(1.+1.2*I2) +2./1.2*log(1+1.2*I2))
	F2=0

	DO Ic = 1, NC
	DO Ia = 1, Na
	F2=F2+mc(Ic)*Ma(Ia)*(BS(IC,Ia)+2*Z*CS(Ic,IA) )
	enddo		
	enddo

	xm=mc(1)
	J1=1
	J2=2
	call EFUNC(J1,J2,Aphi,I,E,ED)
	
	F3= 0

	DO Ic1 = 1, NC
	DO Ic2 = IC1+1, NC	
	z1=ZC(IC1)
	z2=ZC(IC2)
	IF(z1.eq.z2) goto 21
	F3=F3+ED*MC(Ic1)*MC(IC2)
21	continue
	enddo
	enddo
	f4=0.


	DO IA1 = 1, NA
	DO IA2 = IA1+1, NA	
	z1=ZA(IA1)
	z2=ZA(IA2)
	IF(z1.eq.z2) goto 22
	F4=F4+ED*MA(IA1)*MA(IA2)
22	continue
	enddo
	enddo
	

	F=F1+F2+F3+f4


	IF (N.le.NC) then

	F3=0
	DO Ic1 = 1, NC
	z1=ZC(N)
	z2=ZC(IC1)
	IF(z1.eq.z2) goto 121
	F3=F3+E*MC(Ic1)
121	continue
	enddo

	a1=zc(N)**2 *F
	a2=0
	DO Ia=1,NA
	a2=a2+MA(IA)*( 2*B(N,IA)+Z*C(N,IA) )
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
	a3=a3+MA(IA)*MC(IC)*C(IC,IA)
	enddo
	enddo

	a3=Zc(N)*A3
	gam=a1+a2+a3 +F3


	else 

	N1=N-Nc
	a1=zA(N1)**2 *F

	f4=0.
	DO IA1 = 1, NA
	z1=ZA(N1)
	z2=ZA(IA1)
	IF(z1.eq.z2) goto 122
	F4=F4+E*MA(IA1)
122	continue
	enddo

	
	a2=0

	DO IC=1,NC
	a2=a2+MC(IC)*( 2*B(IC,N1)+Z*C(IC,N1) )
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
	a3=a3+MA(IA)*MC(IC)*C(IC,IA)
	enddo
	enddo

	a3=ZA(N1)*A3
	gam=a1+a2+a3+f4
	endif
	gammann= exp(gam)
	return
	end
        

c     -----------------------------------------------
c
c     calcualtes water activity of ionic solution using
c     Pitzer ion interactions model. 
c     
c     -----------------------------------------------
      

	function gammasn(T,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

	implicit real*8 (a-h,o-z)

	Parameter( NMAx=5 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax),C0(nmax,nmax),c1(nmax,nmax)
	real*8 Bphi(nmax,nmax),C(nmax,nmax)
	real*8 I, I2, MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 omega(nmax,nmax)
	

C	Calculate I
	alpha=2.
	 I=0
	DO IC=1,NC
	I=I+mC(IC)*zC(IC)**2
	enddo
	DO IA=1,NA
	I=I+mA(IA)*zA(IA)**2
	enddo
	I=I/2.
	I2=sqrt(I)

CCC	Calculate Bphi, C
	x=sqrt(I)*alpha

	DO Ic = 1, NC
	DO Ia = 1, Na
	Bphi(IC,IA) = b0(Ic,Ia) +exp(-x)*b1(Ic,Ia)

	omega1=omega(Ic,IA)
        xo=omega1*I2
	C(Ic,Ia) = C0(Ic,Ia)+c1(Ic,IA)*exp(-xo)
	enddo
	enddo

	Aphi=.377+4.684E-4*(T-273.15)+3.74e-6*(T-273.15)**2
C	calculate sum mi and Z

	z=0	
	xmi=0
	DO IC=1,NC
	Z=Z+ZC(IC)*MC(IC)
	xmi=xmi+MC(IC)
	enddo
	DO IA=1,NA
	Z=Z+ZA(IA)*MA(IA)
	xmi=xmi+MA(IA)
	enddo


	fphi1=  -Aphi*I**(3./2.)/(1+1.2*I2)

	xs=0
	DO Ic=1,NC
	DO Ia=1,Na
	xs=xs+MA(IA)*MC(IC)*(Z*C(IC,IA)+Bphi(IC,IA))
	enddo
	enddo

	J1=1
	J2=2
	call EFUNC(J1,J2,Aphi,I,E,ED)
	pp=e+I*ed	
	F3 = 0

	DO Ic1 = 1, NC
	DO Ic2 = IC1+1, NC	
	z1=ZC(IC1)
	z2=ZC(IC2)
	IF(z1.eq.z2) goto 21
	F3=F3+pp*MC(Ic1)*MC(IC2)
21	continue
	enddo
	enddo
	f4=0.


	DO IA1 = 1, NA
	DO IA2 = IA1+1, NA	
	z1=ZA(IA1)
	z2=ZA(IA2)
	IF(z1.eq.z2) goto 22
	F4=F4+pp*MA(IA1)*MA(IA2)
22	continue
	enddo
	enddo
	xmh0=0
	DO IA=3,NA
	xmh0=xmh0+mA(IA)
	enddo

	phix=  fphi1 + xs  +f3+f4
	phi=1+ phix*2/xmi
	as=(-phi)*18/1000.*xmi
        gammasn=exp(as)
	return
	end

      
c     ----------------------------------------------------------
c      required by the Pitzer model gammann and gammasn
c     ----------------------------------------------------------
      function gs(X)
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 gs
	gs=2*(-1.+(1+x+x**2/2.)*exp(-x) )/x**2
	return
	end


c     ----------------------------------------------------------
c      Linear interpolation subroutine
c     ----------------------------------------------------------

      subroutine intpl(x1,y1,n1,x2,y2,n2)
           implicit real*8 (a-h,o-z)
           REAL*8 X1(*),X2(*),y1(*),y2(*)

           DO I=1,n2

              if(x2(i).le.x1(1)) then
                 y2(i)=y1(1)
                 goto 20
              endif

              if(x2(i).ge.x1(n1)) then
                 y2(i)=y1(n1)
                 goto 20
              endif


              DO J=2,n1
                  if(x1(J).ge.x2(I)) then
                     ff= (x2(I)-X1(J-1))/(x1(J)-x1(J-1))
                  yy=y1(J-1)+ff*(y1(J)-y1(J-1))
                  y2(I)=yy
                  goto 20
                  ENDIF

             enddo

 20           continue

           enddo



           return
           end

c     ----------------------------------------------------------
c      subroutine to find the zero place of function func
c     ----------------------------------------------------------

	subroutine dzbren(func, erabs,tol,x1,x2,ITMAX)
	implicit real*8 (a-h,o-z)
	INTEGER ITMAX
	INTEGER iter
	REAL*8 tol, x1, x2, func,EPS
	EXTERNAL func
	PARAMETER (EPS=1.D-14)
	save fa,fb,a,b,c,fc,xm,d,e
	
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)	
	if((fa.gt.0. .and.fb.gt.0.).or. (fa.lt.0. .and.fb.lt.0.))
     * pause 'root muss be bracketed for zbrent'
	c=b
	fc=fb	
	do  iter=1,ITMAX
	if((fb.gt.0. .and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))	then
	c=a	!Rename a, b, c and adjust bounding interval d.
	fc=fa
	d=b-a
	e=d
	endif
	if(dabs(fc) .lt.dabs(fb)) then	
	a=b
	b=c
	c=a
	fa=fb
	fb=fc
	fc=fa
	endif
	toli=2.*EPS*dabs(b)+0.5*tol !     Convergence check.
	xm= .5*(c-b)
	if (dabs (xm) .le.toli .or. fb.eq.0.)then
	x2=b
	return
	endif
	if(dabs(e) .ge.toli .and. dabs(fa) .gt.dabs(fb)) then
	s=fb/fa	 !Attempt inverse quadratic interpolation.
	if(a.eq.c) then
	p=2.*xm*s
	q=1D0-s
	else
	q=fa/fc
	r=fb/fc
	p=s*(2. *xm*q* (q-r)-(b-a)*(r-1.))
	q=(q-1.)*(r-1.)*(s-1.)
	endif
	if(p.gt.0.) q=-q!	Check whether in bounds.
	p=dabs (p)
	if(2.*p .lt. min(3.D0*xm*q-dabs(toli*q),dabs(e*q))) then
	e=d
	d=p/q !Accept interpolation.
	else
	d=xm !Interpolation failed. use bisection.b
	e=d
	endif
	else !Bounds decreasing too slowly, use bisection.
	d=xm
	e=d
	endif
	a=b	!Move last best guess to a.
	fa=fb
	if(dabs(d) .gt. toli) then 
	   b=b+d
	else
	   b=b+dsign(toli,xm)
	endif
	fb=func(b)
	enddo 
	pause 'dzbrent exceeding maximum iterations'
	x2=b
	return
	END



c     ----------------------------------------------------------
c      subroutine to find the zero place of function func, for the case when they used twice
c     ----------------------------------------------------------
      
	subroutine dzbrens(func, erabs,tol,x1,x2,ITMAX)
c	FUNCTION dzbren(func,x1 ,x2,tol)
	implicit real*8 (a-h,o-z)
	INTEGER ITMAX
	INTEGER iter
	REAL*8 tol, x1, x2, func,EPS
	EXTERNAL func
	PARAMETER (EPS=1.D-14)
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)	
	if((fa.gt.0. .and.fb.gt.0.).or. (fa.lt.0. .and.fb.lt.0.))
     * pause 'root muss be bracketed for zbrent'
	c=b
	fc=fb	
	do  iter=1,ITMAX
	if((fb.gt.0. .and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))	then
	c=a	!Rename a, b, c and adjust bounding interval d.
	fc=fa
	d=b-a
	e=d
	endif
	if(dabs(fc) .lt.dabs(fb)) then	
	a=b
	b=c
	c=a
	fa=fb
	fb=fc
	fc=fa
	endif
	toli=2.*EPS*dabs(b)+0.5*tol !     Convergence check.
	xm= .5*(c-b)
	if (dabs (xm) .le.toli .or. fb.eq.0.)then
	x2=b
	return
	endif
	if(dabs(e) .ge.toli .and. dabs(fa) .gt.dabs(fb)) then
	s=fb/fa	! Attempt inverse quadratic interpolation.
	if(a.eq.c) then
	p=2.*xm*s
	q=1D0-s
	else
	q=fa/fc
	r=fb/fc
	p=s*(2. *xm*q* (q-r)-(b-a)*(r-1.))
	q=(q-1.)*(r-1.)*(s-1.)
	endif
	if(p.gt.0.) q=-q!	Check whether in bounds.
	p=dabs (p)
	if(2.*p .lt. min(3.D0*xm*q-dabs(toli*q),dabs(e*q))) then
	e=d
	d=p/q !	Accept interpolation.
	else
	d=xm !Interpolation failed. use bisection.
	e=d
	endif
	else !Bounds decreasing too slowly, use bisection.
	d=xm
	e=d
	endif
	a=b	!Move last best guess to a.
	fa=fb
	if(dabs(d) .gt. toli) then 
	   b=b+d
	else
	   b=b+dsign(toli,xm)
	endif
	fb=func(b)
	enddo 
	pause 'dzbrent exceeding maximum iterations'
	x2=b
	return
	END
	
c     ---------------------------------------------
c      
c     Generate the parameters of ion pairs for Pitzer ion interaction model
c
c     
c     --------------------------------------------
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this program is valid only for the interaction, where at least C
C     one ion is single charged.                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



  
C     Iflag              Molality range          T-range
C     1: H-NO3            0-16                   190-320
C     2: H-Cl             0-60                   190-320
C     3: H-HSO4           0-40                   190-320 
C     4: H-SO4            0-40                   190-320
C     5: NH4-HSO4         0-30                   220-320
C     6: NH4-SO4          0-30                   220-320
C     7: NH4-NO3          0-10                   220-320
C     8: NH4-Cl           0-7                    220-320 
C     9: Na-SO4		  0-4                    220-320
  


C	xs is the mixture parameter
C	1: psi(H,HSO4,NO3)	
C	2: theta(SO4,NO3)	
C 	3: psi(H,SO4,NO3)
C	4: psi(H,HSO4,Cl)	
C	5: theta(HSO4,Cl)	
C 	6: psi(H,SO4,Cl)
CC	7: psi-NH4-SO4-HSO4
CC	8: psi-NH4-H-SO4
CC	9: psi-NH4-H-HSO4
CC     10: Phi-NH4-H
CC     11: psi NH4-NO3-SO4 
CC    


   	subroutine calpar(T,NC,NA,b0,b1,C0,C1,omega,xs,Iflag)
  
	Parameter( NMAx=5 )
        IMPLICIT REAL*8(A-H,O-Z)
  	real*8 b0(nmax,nmax),B1(nmax,nmax)
 	real*8 c0(nmax,nmax),C1(nmax,nmax)
	integer Iflag(NMAX,NMAX)
 	real*8 bb(50),b2(50),b3(50),b8(50)
	real*8 b5(50),b6(50),b7(50),omega(5,5),xs(100),b9(50)

CCCCCCCCCCCC data for H-Cl CCCCCCCCCCCCCCCCCCCC
        data (b2(I),i=1,21) /
     > 0.23378,-7.21238E-02,-1.7335667E-02,5.760665E-03,-8.29279E-03,
     >0.2897,7.575434E-02,-1.1474E-03,0.38038,-0.309442,-2.794885E-03,
     >2.309349E-04,9.322982E-04,-2.398E-04,2.85959E-04,-0.21154,
     > 0.101481,5.945618E-02,-0.107864, 8.81749E-02,  1.9916/
       
CCCCCCCCCCCCC data for H-NO3 CCCCCCCCCCCCCCCCC
        data (bb(I),I=1,27)/
     * 3.895835E-03,-1.55571E-02,1.703729E-02,
     * -5.6173712E-03,  5.732047E-03,  0.91622,  0.613523,
     * -0.68489, 0.3038,  -0.32888, 7.6086113E-07, 7.2714678E-05,
     * -1.0037E-04,3.475E-05,-3.62927E-05,5.380465E-02,-2.2163E-02,
     * -1.0166E-02, 6.5423E-03,  -8.80248E-03, 0.907342,
     *         -6.78428E-4,9.576E-4,2*0D0,  7.769E-3, -5.819E-4 /

c      !mixturedata
   
CCCCCCCCCCCC data for H-HSO4, H-SO4 CCCCCCCCCCCCCCCCCCCC
        data (b3(I),i=1,42) /  0.148843, -7.769E-2 ,  2.8062E-2 ,
     *  4.7903E-4,  7.25E-4 ,  0.17843,  0.678,  8.7381E-2, 
     *  -0.57881,  7.58E-2 , -9.878E-4 ,  5.447651E-4, -2.58798E-4,
     *   1.8466527E-5 ,  1.23457E-5 ,  0.37138, -9.24874E-2 ,
     *  -9.21372E-3 , -1.065158E-2 ,  5.4987733E-2 ,  0.2726312,
     * -1.34824E-3 , -0.24711,  1.25978E-2 ,  0.11919,  0.7397,
     *  -3.01755,  -4.5305,  -3.1072, -0.8555842,  9.2223E-4 ,
     * -4.1694532E-3 ,  7.141266E-3 ,  2.32984E-3 , -6.98191E-4 ,
     *  -2.242,  0.71925,   2.52, -0.7391,  -1.548503, 1.5452, 2./


CCCCCCCCCCCCCCCC NH4-HSO4 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     with out chan's data

c	data (b5(I),I=1,32)/
c     * -8.746E-4,  -2.3125, -9.56785E-6,   2.58238,   2.38, -3.1314E-4,
c     *  1.6896E-2, -0.7351,  0.6883,  1.813E-3, -0.1012515 ,
c     * -2.66E-2, -2.86617E-3,  0.22925,  0.438188,  2.522E-4,
c     * -2.90117E-5,  0.9014,  0.41774,  -1035.9,  0.0,  -299.69,
c     *  0.0, -4.9687E-4,  0.0,  1.21485E-2,  0.0, -1.0334E-3,
c     *  0.0,  8.48374E-2,  0.0,0d0/

	data (b5(I),I=1,32)/
     &    -0.11671E-01,  
     2   -0.49047E+00 ,  
     3    0.10378E-03 ,  
     4    0.10170E+00 ,  
     5    0.11722E+01 ,  
     6   -0.30795E-03 ,  
     7   -0.20000E-03 ,  
     8   -0.35877E+01 ,  
     9    0.23663E+01 ,  
     1   -0.10000E-03 ,  
     1    0.50000E-04 ,
     &	  0.169702108786968   
     *, -0.387640333038779    
     *,  0.0
     *,  0.0
     *,  2.601561402128945E-004
     *, -2.337658723484350E-003
     *,  0.0
     *,  0.0
     *,  -612.096813408415 
     *,  0.00
     *,  -306.048406704207     
     *,  0.00
     *,  5.320349188422027E-004
     *,  0.0
     *,  0.0
     *,  0.0
     *,  3.532241943713784E-003 
     *,  0.0
     *, -4.534300563622422E-002 
     *,  0.0000
     *,  0.0000/
	data ikk /0/
	save ikk
	
	ikk=1
	    
	     if (ikk.eq.0) then
	   ikk=1
	   open(333,file='/home/bluo/lib/b.mix')
	   b5=0

	   DO I=1,11
	      read(333,*) iii, b5(I)
	   enddo

	   DO I=12,32
c	      read(333,*) iii, b5(I)
	   enddo

         endif


      
C      ! with chan's data
c	data (b5(I),I=1,31)/
c     * -7.8224E-3,  -1.722,  7.5882E-5,  0.962,   1.8914, -4.285E-4,
c     *  3.99E-4, -0.76151,  0.53133,  1.16726E-3,
c     * -9.45855E-2, -2.007E-2,  6.0482E-3, -0.1392,
c     *   2.9544,  2.45133E-4,   -5.492E-5,  0.8358,
c     * -0.3791,  -980.81,  0.,   265.7,  0., -1.7568E-3,
c     * -7.72E-4,  3.437E-4,   -1.6456E-2, -1.2034E-3,
c     * -2.2943E-3,  7.9281E-2,  3.772E-2/

CCCCCCCCCCCCCCCC NH4-SO4 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c	data (b6(I),I=1,13) /-1.2058223E-2 ,1.1043,4.79018E-5
c     * ,2.14346E-2 , 0.58,  -2.9146E-2,  1.9631E-4  ,   1.1378,
c     *  0.9283,1.28548E-4,1.684E-5,2.6267E-2, -2.6E-4/   ! wt=1
 
C	data (b6(I),I=1,13)/
C     * -1.2058223E-2,   1.1043,  4.79018E-5,2.14346E-2 ,  0.58,
C     * -0.1188,  8.5E-2,   2.10514,  0.5942,  7.888E-4, -5.503E-4,
C     *  5.815E-2 , -3.766E-2 / !wt=5d-5

	data (b6(I),I=1,13)/
     *  -4.327681689677379E-003 
     *      ,  0.953787255317688      
     *      ,  1.113167041454729E-005 
     *      ,  2.101340928698026E-002 
     *      ,  0.614050552546192      
     *      , -2.606487385557435E-003 
     *      ,  5.672030237892362E-004 
     *      ,   1.75633200902706      
     *      ,   1.10907828739749      
     *      ,  1.177526981905975E-005 
     *      ,  6.713875789826932E-007 
     *      ,  5.973238697090055E-003 
     *      , -2.812698651271791E-003/ 
	 

	IKK1=1
	if (ikk1.eq.0) then

	   ikk1=1
	   open(333,file='/usr/users/luo/amsulf/bo.nh42so4t')
	   DO I=1,13
	      read(333,*) iii, b6(I)
c	      print*, I, b6(I)
	   enddo
	endif


CCCCCCCCCCCCC data for NH4-NO3 CCCCCCCCCCCC
	data (b7(I),I=1,13) /-2.3275E-2, 0.15, 1.1634E-4,1.62E-3,
     * 0.43, ! 0.107264,  0.221, -3.8442E-4, -1.3872E-2, 4*0D0/
     *  8.78E-2 ,  0.2753645, -3.349E-4 , -1.093E-2 ,
     * -4.769E-2 ,  0.1776,  1.25E-4 ,  6.9751E-3 /

C
C     *  2.15438E-2,0.67073,  2*0.0, -1.3662E-2,3.4747E-002,2*0.0/
 

CCCCCCCCCCCCCCCC NH4-Cl CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	data (b8(I),I=1,9) / -6.333E-4 , -3.99546E-4 , 0.3155,
     *  0.1414, -3.837E-5, 1.08331E-4, 5.2436E-2,1.6827E-2,1.19/

CCCCCCCCCCCCC data for Na-SO4 CCCCCCCCCCCC
	data (b9(I),I=1,17) /1.63E-03 , 2.092E-03,3.484156E-02 ,
     * -1.057E-02,1.0775,  0.9, 0.8206,-9.6425E-02,2.7492E-03,
     *-9.2838E-04,-6.8268E-04,4.8126E-05,7.182E-02,
     *  0.7586, -0.30291, 0.2311,1.7/

	xs(4)=0
	xs(5)=0
	xs(6)=0
	DT=(T-298.15)/100D0

	DO I=1,NC
	DO J=1,NA
	if (Iflag(I,J) .eq. 2) then

	omega(i,j)=b2(21)
        B0(i,j)=b2(1)+dt*b2(2)+dt**2*b2(3)+dt**3*b2(4)+dt**4*b2(5)
        B1(i,j)=b2(6)+dt*b2(7)+dt**2*b2(8)+dt**3*b2(9)+dt**4*b2(10)
        c0(i,j)=b2(11)+dt*b2(12)+dt**2*b2(13)+dt**3*b2(14)+dt**4*b2(15)
        c1(i,j)=b2(16)+dt*b2(17)+dt**2*b2(18)+dt**3*b2(19)+dt**4*b2(20)
        endif

	 IF (Iflag(I,j).eq.7) then
	b0(I,J)=b7(1)+b7(6)*dt+b7(10)*dt*dt
	b1(I,J)=b7(2)+b7(7)*dt+b7(11)*dt*dt
	c0(I,J)=b7(3)+b7(8)*dt+b7(12)*dt*dt
	c1(I,J)=b7(4)+b7(9)*dt+b7(13)*dt*dt
	omega(I,j)=b7(5)

	xs(1)=bb(22)+bb(23)*dt
	xs(2)=bb(24)+bb(25)*dt
	xs(3)=bb(26)+bb(27)*dt
        xs(11)=4.75458E-4-4.0577E-003*dt
c4.633E-4-4.093E-3*dt
        endif

	IF (Iflag(I,J).eq.9) then

        B0(i,j)=b9(1)+dt*b9(2)+dt**2*b9(3)+dt**3*b9(4)
        B1(i,j)=b9( 5)+dt*b9(6)+dt**2*b9(7)+dt**3*b9(8)
        c0(i,j)=b9( 9)+dt*b9(10)+dt**2*b9(11)+dt**3*b9(12)
        c1(i,j)=b9(13)+dt*b9(14)+dt**2*b9(15)+dt**3*b9(16)
	omega(I,j)=b9(17)
        endif

	IF (Iflag(I,j).eq.8) then

	omega(I,j)=b8(9)
	b0(I,J)=b8(1)+b8(2)*dt
	b1(I,J)=b8(3)+b8(4)*dt
	c0(I,J)=b8(5)+b8(6)*dt
	c1(I,J)=b8(7)+b8(8)*dt
        endif

	 IF (Iflag(I,j).eq.3) then
	dt=(t-298.15)/100.
	omega(I,J)=b3(41)
	omega(I,J+1)=b3(42)
        B0(i,j)=b3(1)+dt*b3(2)+dt**2*b3(3)+dt**3*b3(4)+dt**4*b3(5)
        B1(i,j)=b3(6)+dt*b3(7)+dt**2*b3(8)+dt**3*b3(9)+dt**4*b3(10)
        c0(i,j)=b3(11)+dt*b3(12)+dt**2*b3(13)+dt**3*b3(14)+dt**4*b3(15)
        c1(i,j)=b3(16)+dt*b3(17)+dt**2*b3(18)+dt**3*b3(19)+dt**4*b3(20)
       B0(i,j+1)=b3(21)+dt*b3(22)+dt**2*b3(23)+dt**3*b3(24)+dt**4*b3(25)
       B1(i,j+1)=b3(26)+dt*b3(27)+dt**2*b3(28)+dt**3*b3(29)+dt**4*b3(30)
       c0(i,j+1)=b3(31)+dt*b3(32)+dt**2*b3(33)+dt**3*b3(34)+dt**4*b3(35)
       c1(i,j+1)=b3(36)+dt*b3(37)+dt**2*b3(38)+dt**3*b3(39)+dt**4*b3(40)
        endif

	if ( Iflag(I,J) .eq. 1) then
	omega(i,j)=bb(21)
	dt=(T-298.15)/100.
        B0(i,j)=bb(1)+dt*bb(2)+dt**2*bb(3)+dt**3*bb(4)+dt**4*bb(5)
        B1(i,j)=bb(6)+dt*bb(7)+dt**2*bb(8)+dt**3*bb(9)+dt**4*bb(10)
        c0(i,j)=bb(11)+dt*bb(12)+dt**2*bb(13)+dt**3*bb(14)+dt**4*bb(15)
        c1(i,j)=bb(16)+dt*bb(17)+dt**2*bb(18)+dt**3*bb(19)+dt**4*bb(20)


	xs(1)=bb(22)+bb(23)*dt
	xs(2)=bb(24)+bb(25)*dt
	xs(3)=bb(26)+bb(27)*dt
          xs(11)=4.75458E-4-4.0577E-003*dt
c4.633E-4-4.093E-3*dt
        endif

	if ( Iflag(I,J) .eq. 5) then

	B0(i,j)=b5(1)+b5(11+1)*(dt)+b5(11+2)*dlog(t/298.15)
	B1(i,j)=b5(2)+b5(11+3)*(dt)+b5(11+4)*dlog(t/298.15)
	c0(i,j)=b5(3)+b5(11+5)*(dt)+b5(11+6)*dlog(t/298.15)
	c1(i,j)=b5(4)+b5(11+7)*(dt)+b5(11+8)*dlog(t/298.15)
	omega(i,j)=b5(5)+b5(32)*dt

	xs(7)=b5(6)  +b5(11+13)*(dt)+b5(11+14) *dlog(t/298.15)
	xs(8)=b5(7)  +b5(11+15)*(dt)+b5(11+16) *dlog(t/298.15)
	xs(9)=b5(10) +b5(11+17)*(dt)+b5(11+18) *dlog(t/298.15)
	xs(10)=b5(11)+b5(11+19)*(dt)+b5(11+20) *dlog(t/298.15)
c	xs(10)=0d0

	T0=298.15

        sld1 = -86.+2791.9/T+13.482*dlog(T)
	sld2=b5(8)+(1/t-1/298.15)*b5(11+9)+b5(11+10)*log(t/t0)
	sld3=b5(9)+(1/t-1/298.15)*b5(11+11)+b5(11+12)*log(t/t0)
        endif

	if ( Iflag(I,J) .eq. 6) then

	omega(I,J)=b6(5)
	b0(I,J)=b6(1)+b6(6)*dt+b6(7)*dt*dt
	b1(I,J)=b6(2)+b6(8)*dt+b6(9)*dt*dt
	C0(I,J)=b6(3)+b6(10)*dt+b6(11)*dt*dt
	C1(I,J)=b6(4)+b6(12)*dt+b6(13)*dt*dt
                   dd=00

	endif


	enddo
	enddo

	xs(11)=0d0
	
	return
	end

c     ---------------------------------
c
c     vapor pressure of pure water Koop and Murphy
c
c     -----------------------------------
      
      function vwater(temp)
      implicit real*8 (a-h,o-z)
      vwater= dexp(54.842763 - 6763.22/temp - 4.210*dlog(temp) +
     + 0.000367*temp + dtanh(0.0415*(temp - 218.8))*(53.878 -
     + 1331.22/temp - 9.44523*dlog(temp) + 0.014025*temp))
	vwater=vwater/100d0

	return
        End

c     ---------------------------------
c
c     required by Pitzer model gammann and gammasn
c
c     -----------------------------------

	function g(x)
      IMPLICIT REAL*8(A-H,O-Z)
      		real*8 g
      	g=2*(1.-(1.+x)*exp(-x))/X**2
	return
	end

c     ---------------------------------
c
c     required by Pitzer model gammann and gammasn
c
c     -----------------------------------

	SUBROUTINE EFUNC(Icharg,JCHARG,A,XI,E,ED)
       IMPLICIT REAL*8(A-H,O-Z)
		REAL*8	X(3),J0(3),J1(3),DUM

	IF((ICHARG.EQ.JCHARG) .OR. (XI .LE. 1.D-30) )THEN
	E=0
	ED=0
		ELSE

	X(1)=6*ICHARG*JCHARG*SQRT(XI)*A
	X(2)=6*ICHARG*ICHARG*SQRT(XI)*A
	X(3)=6*JCHARG*JCHARG*SQRT(XI)*A

	DO I=1,3		
	DUM=-1.2D-2*X(I)**.528
	J0(I)=X(I)/(4.+4.581*X(I)**(-.7238)*EXP(DUM))
	J1(I)=(4.+4.581*X(I)**(-.7238)*EXP(DUM)*(1.7238-DUM*.528
     * ))/(4.+4.581*X(I)**(-.7238)*EXP(DUM))**2
	ENDDO
	
	XE=ICHARG*JCHARG/4./XI*(J0(1)-.5*J0(2)-.5*J0(3))
	XED=ICHARG*JCHARG/8./XI**2*(X(1)*J1(1)-.5*X(2)*J1(2)-.5*
     * X(3)*J1(3))-XE/XI
	E=XE
	ED=XED
	ENDIF
	
	RETURN
	END


c     --------------------------------------------------------
c
c      Reset shells to equal volume when one shell becomes too thin or too thick
c
c     -------------------------------------------------------
      subroutine reset_shells(NS,x,xn)
       implicit real*8 (a-h,k,m, o-z)
      parameter (NSMM=1000)
      parameter (np=29)         ! number of species
      real*8 xn(NSMM,np),xn0(NSMM,np),vshell(NSMM)

      real*8 xn2(NSMM,np),ml(np)
      real*8 x(*), x0(NSmm+1)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      common /M/ MM,mv,izc
c      common /vshell/ vshell
      common /pi/ pi

      
     

c     set equidstance shells to equil volume

c      print*, 'diff tau = ', difftau
      volmin=1D8
      vommax=0
      
      DO I=1,NS
         vv = mv(1)*xn0(I,1)+mv(2)*xn0(I,2)

          DO j=6,np
          if (j.lt.21 .or. j.gt. 22)  
     & vv = vv +mv(j)*xn0(I,j)
             enddo
             
         if (vv.gt. volmax) volmax= vv
         if (vv.lt. volmin) volmin= vv
         
         enddo
         

         if (volmax.ge. 3d0* volmin) then
            print*, ' Rebin '

            write(28,'(A)') ' rebin'


       DO I=1,ns
          DO J=1,NP
          xn0(I,j)=xn(i,j)
          xn(I,j)=0d0
          enddo
       enddo

c          x0(ns+1)=x(ns+1)

            vol=0
       DO I=1,ns
          vol = vol+mv(1)*xn0(I,1)+mv(2)*xn0(I,2)
          vv = mv(1)*xn0(I,1)+mv(2)*xn0(I,2)

          DO j=6,np
          if (j.lt.21 .or. j.gt. 22)  
     & vol = vol +mv(j)*xn0(I,j)
          if (j.lt.21 .or. j.gt. 22)  
     & vv = vv +mv(j)*xn0(I,j)

             enddo
             x0(I+1) = (vol/4d0/pi*3d0 )**(1/3d0)
c             print*,i,x0(I+1), x(I+1)
c             print*,vshell(I), vv

             enddo

             vol1=vol/NS

             print*,' vol1 = ', vol1

         X(1)=0d0
     

       DO I=1,NS
C     CALCULATE XN1,XN2
          X(I+1)= (vol1*I/4d0/pi*3d0)**(1/3d0)
          
          RI = X(I)
          RA = X(I+1)

          DO J=1 ,ns
             if (x0(J+1).ge. RI) goto 11
          enddo
 11       NI= J

          DO J=1 ,ns
             if (x0(J+1).ge. RA) goto 12
          enddo
 12       NA= J
           if (NA.gt.NS)NA=NS

c          print*,I, NI,NA,x(I+1)



          if (NA.eq.NI)  then
          vnew = (RA**3  - RI**3)
          vS=( x0(NI+1)**3 -x0(NI)**3 )
          DO J=1,NP
          xn(I,j)= vnew/vs* xn0(NI,j)
          enddo
          goto 13
          endif


c     take the RI to x(NI+1)
          vnew= x0(NI+1)**3 -rI**3 
          vS= x0(NI+1)**3 -x0(NI)**3 

          DO j=1,Np
          xn(I,J)=xn0(NI,j)*vnew/vs
          enddo

c     take x(NA) to RA
          vnew=( RA**3 -x0(NA)**3 )
          vS=( x0(Na+1)**3 -x0(Na)**3 )

          DO j=1,Np
          xn(I,j)=xn(I,j)+ xn0(NA,j)*vnew/vs
          enddo

C     TAKE THE SHELLS NI+1 TO NA-1
          DO jJ = ni+1,NA-1
             DO j=1,np
          xn(I,j) =xn(I,j)+xn0(jj,j)
          enddo
          enddo
13       continue
       ENDDO



       DO J=1,NP
          xns0=0d0
          xns=0d0
         DO I=1,ns
            xns=xns+ xn(i,j)
            xns0=xns0+ xn0(i,j)
          enddo
          if (xns0.gt.0d0) print*,j, xns, (xns-xns0)/xns0
c          if (xns0.gt.0d0) write(28,*) j, xns, (xns-xns0)/xns0
          enddo
          print*, 'finish rebin'          
c          stop

         vol=0d0
         vol0=0d0

          DO I=1,NS

          enddo

            vol=0
            xn0=xn

            T=298.15
            ML(1)=1000D0/mm(1)
            DO I=1, NS
               DO J=1,np
                  ML(J)= xn(I,j)*ml(1)/xn(I,1)
                  enddo

           call calHNew(T,ml)
                      call aw_back
     & (t,Ml,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa,
     & gammas1,gammas2)        
                      phh= dlog(mL(6)* gammah)/dlog(.1d0)
                      print*, 'phh = ', i,phh

               enddo
          endif

       return
       end
 

c     --------------------------------------------------------
c
c      Merge shell when it becomes very thin
c
c     -------------------------------------------------------

      subroutine merge(xn,x,NS)
       IMPLICIT REAL*8(A-H,O-Z)
       parameter (np=29)
      parameter (NSMM=1000)
      real*8 MM(NP) ! mol Mass
      real*8 Mv(NP) ! mole volume
      real*8 xn(NSMM, np),xna(NSMM, np)
      real*8 x(NSMM+1),vol(NSMM)
      common /M/ MM,mv
      common /pi/pi

c     not possible to merge
      if (NS.eq.1) return

      
      DO I=1,NS
         vol(I)= 4d0*pi/3d0*(x(I+1)**3-x(I)**3)
       enddo

         DO I=1,NS
c     compare shell below 
      ismerge =0

            if (I.gt.1) then
               if (vol(I)/vol(I-1).le. 0.2d0) ismerge=1
               endif

            if (I.lt.NS) then
               if (vol(I)/vol(I+1).le. 0.2d0) ismerge=1
               endif
                     enddo

                     if (ismerge.eq.0) return

                     
         DO I=1,NS
c     compare shell below 
      ismerge =0

            if (I.gt.1) then
               if (vol(I)/vol(I-1).le. 0.2d0) ismerge=1
               endif

            if (I.lt.NS) then
               if (vol(I)/vol(I+1).le. 0.2d0) ismerge=1
               endif

               if (ismerge.eq.1) then

                  if (I.eq.NS) then
                  DO J=1,1,NP
                     xn(NS-1,j)=xn(NS-1,j)+xn(NS,j)
                  enddo
                  NS=NS-1
                  x(NS+1)=x(NS+2)
                  print*, ' The shell'  , I, 'is merged to, shell',  I-1
                  return
                  else

                  DO J=1,1,NP
                     xn(I,j)=xn(I,j)+xn(I+1,j)
                  enddo
                  NS=NS-1
                  
                  DO J=I+1,NS+1
                     x(J)= x(J+1)
                     enddo
                  return
                     endif

                   endif

                     enddo
                     


          return
          end



c     --------------------------------------------------------
c
c     Calcualtes the equilibrium composition of solution with NaCl crystal
c     Ma: the inital compostion as input and output is the composition with NaCl crystal
c     Msalt: is the molalty of NaCl crystal
c      
c     -------------------------------------------------------

      
      subroutine effl_misch(Ta,ma,msalt)

        IMPLICIT REAL*8(A-H,m,O-Z)
        parameter (np=29)
        real*8 fcnmisch
	external  fcnmisch
        real*8 MA(NP),m(NP),msalt
	common/suls/ T, M
        common /kout/xx2

        
        
        t=ta
        M= MA

 	xmin=1D-14              
          	 xmax=m(16)-1D-20        
         xmax1=m(17)-1D-20        
        if (xmax1.le.xmax) xmax= xmax1
        
	erabs=0.d0
	errel=0d0
	ITmax=2000

	xx1= fcnmisch(xmin)
	xx2= fcnmisch(xmax)

c	write(6,'(4E14.6)') xmin,xmax,xx1,xx2
	if(xx1.le.0.) then
	xmax=xmin
	goto 2
      endif
      
	if(xx2.ge.0.) 	goto 2
	call dzbren(fcnmisch, erabs,errel,xmin,xmax,ITMAX)
 2      continue

        Msalt= xmax
        MA(16)=Ma(16)-xmax
        MA(17)=Ma(17)-xmax

      call aw_back(t,Ma,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNa)
              Ap= Ma(16)*ma(17)*gammaCl*gammaNa


              print*, ap
              


        
 999            return
	end  

c     ----------------------------------
c    Required by eff_misch
c
c     ---------------------------------
           function fcnmisch(x)

        IMPLICIT REAL*8 (A-H,m,O-Z)
      integer NP
      parameter (np=29)
      real*8 M(NP), NL(NP),msalt,ML(NP)
      common/suls/ T, M

      
      
      msalt= x
      ML=M
      
      ML(16)=ML(16)-Msalt
      if (ML(16).le.1D-20) ML(16)=1D-20
      
      ml(17)=ML(17)-Msalt
      if (ML(17).le.1D-20) ML(17)=1D-20



      
      call aw_back(t,Ml,aw,gammaH,gammaNO3,gammaNH4,gammaCl,gammaNA)

       Aproduct= ML(16)*ml(17)*gammaCl*gammaNa
      
       fcnmisch= Aproduct -Apnacl(T)
       

      
         
        return
	end

