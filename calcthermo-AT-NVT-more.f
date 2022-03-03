      program thermoanly
      Implicit none
      integer ::k,  Nst,i,step,s
      integer, parameter :: iprint = 2
      double precision, parameter :: nchain = 300 ! number of polymers in simulations
      double precision, parameter :: kB = 1.38065e-23 ! J/K
      double precision, parameter :: NA = 6.022141e23 ! mol-1
      double precision, allocatable::temp(:),ke(:),pe(:),te(:),press(:)
      double precision, allocatable::rg(:)
      double precision, allocatable::epair(:),ebond(:),eang(:),edih(:)
      double precision, allocatable::pang(:),ppair(:),pbond(:),pke(:)
      double precision, allocatable::pvir(:)
      double precision :: sumtemp,sumke,sumpe,sumte,sumpress
      double precision :: avgtemp,avgke,avgpe,avgte,avgpress
      double precision :: sumepair,sumebond,sumeang,sumedih
      double precision :: avgepair,avgebond,avgeang,avgedih
      double precision :: sumpang,sumppair,sumpbond,sumpke
      double precision :: sumpvir
      double precision :: avgpang,avgppair,avgpbond,avgpke
      double precision :: avgpvir
      double precision :: sumrg, avgrg
      double precision :: sumrgsq, avgrgsq
      double precision :: instte,instpress
      double precision :: instrg, avginstrg
      double precision :: instrgsq, avginstrgsq
      double precision :: instepair,instebond,insteang,instedih
      double precision :: instpang,instppair,instpbond,instpke
      double precision :: instpvir
      double precision :: insttesq, sigmaE, Cv
      double precision :: avginstte,avginsttesq
      double precision :: dummy1,t,r
      double precision :: a,b,c,d
      double precision :: e,f,g,h
      
      open(unit=11,file='average-thermo.dat',status='unknown')
      open(unit=12,file='evol-totE.dat',status='unknown')
      open(unit=13,file='evol-avg-pressure.dat',status='unknown')
      open(unit=16,file='evol-inst-pressure.dat',status='unknown')
      open(unit=14,file='more_thermo.dat',status='old')
      open(unit=15,file='evol-Cv.dat',status='unknown')
      open(unit=17,file='evol-rg.dat',status='unknown')
      open(unit=18,file='evol-rgsq.dat',status='unknown')
      open(unit=19,file='evol-avgrg.dat',status='unknown')
      open(unit=20,file='evol-avgrgsq.dat',status='unknown')

      open(unit=30,file='evol-epair.dat',status='unknown')
      open(unit=31,file='evol-ebond.dat',status='unknown')
      open(unit=32,file='evol-eang.dat',status='unknown')
      open(unit=33,file='evol-edih.dat',status='unknown')
      open(unit=34,file='evol-pang.dat',status='unknown')
      open(unit=35,file='evol-ppair.dat',status='unknown')
      open(unit=36,file='evol-pbond.dat',status='unknown')
      open(unit=37,file='evol-pke.dat',status='unknown')
      open(unit=38,file='evol-pvir.dat',status='unknown')
      

      k=0
 10   k=k+1
      read(14,*,end=11) dummy1
      go to 10
 11   k=k-1
      close(14)
      Nst = k

      allocate(temp(Nst))
      allocate(ke(Nst))
      allocate(pe(Nst))
      allocate(te(Nst))
      allocate(press(Nst))
      allocate(rg(Nst))
      allocate(epair(Nst))
      allocate(ebond(Nst))
      allocate(eang(Nst))
      allocate(edih(Nst))
      allocate(pang(Nst))
      allocate(ppair(Nst))
      allocate(pbond(Nst))
      allocate(pke(Nst))
      allocate(pvir(Nst))

      open(unit=14,file='more_thermo.dat',status='old')
      do i=1,Nst
      Read(14,*)step,temp(i),press(i),ke(i),pe(i),te(i),rg(i) !,epair(i),ebond(i),eang(i),edih(i),pang(i),ppair(i),pbond(i),pke(i)
      End do
      close(14)

      open(unit=14,file='more_thermo.dat',status='old')
      do i=1,Nst
      Read(14,*)step,temp(i),press(i),ke(i),pe(i),te(i),rg(i),
     *epair(i),ebond(i),eang(i),edih(i) !,pang(i),ppair(i),pbond(i),pke(i)
      End do
      close(14)

      open(unit=14,file='more_thermo.dat',status='old')
      do i=1,Nst
      Read(14,*)s,a,b,c,d,e,f,g,h,t,r,a,b,c,d,e,f,a,b,c,d,pke(i),pvir(i)
      End do
      close(14)

      open(unit=14,file='more_thermo.dat',status='old')
      do i=1,Nst
      Read(14,*)s,a,b,c,d,e,f,g,h,t,r,a,b,c,d,e,f,a,b,c,pbond(i)
      End do
      close(14)


c       !Calculated average quantities
      sumtemp=0
      sumke=0
      sumpe=0
      sumte=0
      sumpress=0
      sumrg=0
      sumrgsq=0
      sumepair=0
      sumebond=0
      sumeang=0
      sumedih=0
      sumpang=0
      sumppair=0
      sumpbond=0
      sumpke=0
      sumpvir=0
      do i=1,Nst
        sumtemp = sumtemp  + temp(i)
        sumke   = sumke    + ke(i)
        sumpe   = sumpe    + pe(i)
        sumte   = sumte    + te(i)
        sumpress= sumpress + press(i)
        sumrg   = sumrg    + rg(i)
        sumrgsq = sumrgsq  + ( rg(i)*rg(i) )
        sumepair= sumepair + epair(i)
        sumebond= sumebond + ebond(i)
        sumeang = sumeang  + eang(i) 
        sumedih = sumedih  + edih(i) 
        sumpang = sumpang  + pang(i) 
        sumppair= sumppair + ppair(i)
        sumpbond= sumpbond + pbond(i)
        sumpke  = sumpke   + pke(i)
        sumpvir = sumpvir  + pvir(i)  
      End do
      
      avgtemp = sumtemp/dble(Nst)
      avgke   = sumke/dble(Nst)
      avgpe   = sumpe/dble(Nst)
      avgte   = sumte/dble(Nst)
      avgpress= sumpress/dble(Nst)
      avgrg   = sumrg/dble(Nst)
      avgrgsq = sumrgsq/dble(Nst)
      avgepair= sumepair/dble(Nst) 
      avgebond= sumebond/dble(Nst) 
      avgeang = sumeang/dble(Nst) 
      avgedih = sumedih/dble(Nst) 
      avgpang = sumpang/dble(Nst) 
      avgppair= sumppair/dble(Nst) 
      avgpbond= sumpbond/dble(Nst) 
      avgpke  = sumpke/dble(Nst)
      avgpvir = sumpvir/dble(Nst)

c      !Write the evolutions of pressure, total energy
      instte    = 0.0d0
      instrg    = 0.0d0
      instrgsq  = 0.0d0
      insttesq  = 0.0d0
      instpress = 0.0d0
      instepair = 0.0d0 
      instebond = 0.0d0 
      insteang  = 0.0d0
      instedih  = 0.0d0
      instpang  = 0.0d0
      instppair = 0.0d0 
      instpbond = 0.0d0 
      instpke   = 0.0d0     
      
      do i =1,Nst
         write(16,*) i , press(i)
         write(17,*) i , rg(i)
         write(18,*) i , rg(i)*rg(i)
         write(30,*) i , epair(i)
         write(31,*) i , ebond(i)
         write(32,*) i , eang(i)
         write(33,*) i , edih(i)
         write(34,*) i , pang(i)
         write(35,*) i , ppair(i)
         write(36,*) i , pbond(i)
         write(37,*) i , pke(i)
         write(38,*) i , pvir(i)

         instte = instte + te(i)
         instrg = instrg + rg(i)
         insttesq = insttesq + te(i)*te(i)
         instrgsq = instrgsq + rg(i)*rg(i)
         instpress = instpress + press(i)
         if ( mod(i,iprint) .eq. 0 ) then
            avginstte = instte/dble(i)
            avginsttesq = insttesq/dble(i)
            sigmaE = (avginsttesq - (avginstte*avginstte))
            Cv = (sigmaE/(avgtemp**2))*((4.814*1000)**2)/(kB*NA)/nchain ! J/mol/K 
            avginstrg = instrg/dble(i)
            avginstrgsq = instrgsq/dble(i)
            write(12,*) i , instte/dble(i)
            write(13,*) i , instpress/dble(i)
            write(15,*) i , Cv
            write(19,*) i , instrg/dble(i)
            write(20,*) i , instrgsq/dble(i)
         endif
      end do


      write(11,*)'average temp, K',  avgtemp
      write(11,*)'average kinetic energy, kcal/mol', avgke
      write(11,*)'average potential energy, kcal/mol', avgpe
      write(11,*)'average total energy, kcal/mol', avgte
      write(11,*)'average Press, atm',  avgpress
      write(11,*)'average Rg, Ang',  avgrg
      write(11,*)'average Rg^2, Ang',  avgrgsq
      write(11,*)'Cv, J/(mol K)', Cv
      write(11,*)'average pair energy, kcal/mol',     avgepair      
      write(11,*)'average bond energy, kcal/mol',     avgebond      
      write(11,*)'average angle energy, kcal/mol',    avgeang        
      write(11,*)'average dihedral energy, kcal/mol', avgedih           
      write(11,*)'average pair pressure, atm',        avgppair    
      write(11,*)'average bond pressure, atm',        avgpbond   
      write(11,*)'average kinetic pressure, atm',     avgpke      
      write(11,*)'average virial pressure, atm',      avgpvir       

      Close(11)
      Close(12)
      Close(13)
      Close(15)
      Close(16)
      Close(30)  
      Close(31)
      Close(32)
      Close(33)
      Close(34)
      Close(35)
      Close(36)
      Close(37)
      Close(38)
      
      End
