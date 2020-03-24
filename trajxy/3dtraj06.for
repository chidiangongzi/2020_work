c ---------------------------------------------------------------------------c
c   �����ר��Ϊ�Ͼ���ѧ������ѧϵ�±���׼����ȫLagrangeȺ�������й켣ģʽ   c
c                (2006��5��)      (from xytraj)                              c
c ---------------------------------------------------------------------------c 
      program alptraj
	parameter (ix=120,jy=120,kz=25,ixx=ix*jy*kz,mm=1500) 
      dimension tr(ix,jy,kz),p0(ix,jy,kz),qv(ix,jy,kz)
	DIMENSION V(ix,jy,kz),U(ix,jy,kz),W(ix,jy,kz)
	dimension al(kz)           
C  -----------------------------------------------------------
      dimension qc(ix,jy,kz),qr(ix,jy,kz)
	dimension treo(kz),po(kz)
	common /a1/dlx,dly,dlz,dlt
      common /a3/treo,po,al
      common /a4/iii,nn
      common /a9/ntex,itex,jtex,ktex,mtexb
	write(*,*) '  input dlt,dlx,dly,dlz=? '
	write(*,*) 'such as dlt=5.,dlx=50000.,dly=50000.,dlz=50000.,[cm]'
	read(*,*) dlt,dlx,dly,dlz
	write(*,*) '  input from model-output'
  	open(8,file='modelin.dat',form='unformatted',status='unknown') 
 55   format(i5)
      nfile=1
          do ll=1,nfile
      read(8) nn,nm
      read(8) time
      read(8) u,v,w
      read(8) tr,qv,qc,qr,p0
	     enddo
!      write(*,*)'qcmax='
!      write(*,*) maxVal(qc)
c ---------------------------------------------------------------------------------
      call xytraj0(u,v,w,tr,qv,qc,qr,p0)
	stop
	end
c ----------------------------------------------------------------------------- 
      subroutine xytraj0(u,v,w,tr,qv,qc,qr,p0)
      real m
      parameter (ix=120,jy=120,kz=25,ixx=ix*jy*kz,mm=1500)  
      dimension x(ixx),y(ixx),z(ixx),ali(ixx),d(ixx)
      dimension vr(ixx),vh(ixx),con(ixx)
      dimension m(ixx)
      dimension kh(ixx),kbig(ixx)
      dimension jd(ixx),jpr(ixx)
      dimension treo(kz),po(kz),al(kz)
      dimension tr(ix,jy,kz),p0(ix,jy,kz)
      dimension aqc(ix,jy,kz),aqr(ix,jy,kz),sh(ix,jy,kz)
      dimension qv(ix,jy,kz),qb(ix,jy,kz),r(ix,jy,kz),scr(ix,jy,kz)
      DIMENSION V(ix,jy,kz),U(ix,jy,kz),W(ix,jy,kz)           
C  -----------------------------------------------------------
      dimension qc(ix,jy,kz),qr(ix,jy,kz)
c ------
      common /a1/dlx,dly,dlz,dlt
      common /a3/treo,po,al
      common /a4/iii,nn
      common /a9/ntex,itex,jtex,ktex,mtex
      open(66,file='3rptxy.dat',status='unknown')
      open(99,file='3trajxy.dat',form='unformatted',status='unknown')
      open(77,file='xyz0.dat',status='unknown')
	nuw=30
      do 5 i=1,ixx
      x(i)=0.0
      y(i)=0.0
      z(i)=0.0
      jd(i)=1
      jpr(i)=0
      vr(i)=0.0
      vh(i)=0.0
      ali(i)=0.0
      d(i)=0.0
      m(i)=0.0
      kh(i)=0
   5  continue
      do 6 k=1,kz
      treo(k)=tr(1,1,k)
      po(k)=p0(1,1,k)
    6 continue
      do 7 k=1,kz
      do 8 j=1,jy
      do 9 i=1,ix
      sh(i,j,k)=0.0
!      qc(i,j,k)=0.0
!      qv(i,j,k)=0.0
      r(i,j,k)=0.0
      aqc(i,j,k)=0.0
      aqr(i,j,k)=0.0
    9 continue
    8 continue
    7 continue
      cd=0.45
      hbtex=0
      itex=1
      jtex=1
      ktex=3
      mtex=1
      ntex=1
      icp=2
      icr=0
      alg=0.9
      if(alg.ge.0.9) icr=1
      a=30.0
      iwrq=0
c   the unit of hd is cm.
      nf=0
      als=0.6
      alf=79.3
      alm=0.5
      aln=0.3
      vlfa=0.6
c      vlfa=1.0
      clfa=0.8
      rlfa=1.0
      ulfa=1.0
      wlfa=1.0
      dlt=5.0
      pai=3.14159264
      alio=0.9
      ecr=1.0e-6
      epcr=1.0e-5
      dlx=50000.0
      dly=50000.0
      dlz=50000.0
      iii=0
      nn=0
c ------------  get temp. , presure feild:
      do k=1,kz
      al(k)=po(k)/(2840.0*treo(k))
      enddo
c ---------------------------
      ni=ix
      if(ntex.eq.1) ni=ix
c ----- readin macro & micro data of hail cloud feild -----------
c      write(*,*) ' input nfile=? (i5) which is records of wmax case =3'
c      read(5,55) nfile
      write(*,*) ' input mstep=? (i5) mstep<= mm=1500'
      read(5,55) mstep
	write(*,*) ' input hd=? (f8.2) ,such as 0.05 cm.'
      read(5,56) hd
	write(*,*) ' input dbk=? (f8.2) ,such as 0.4-->1.0 cm.'
      read(5,56) dbk
      write(*,*) ' input dsh=? (f8.2) d>dsh to show information'
      read(5,56) dsh
      write(*,*) ' input hcon=? (f8.2) concentrtion of embryo ref:1.0 '
      write(*,*) '  the unit of hcon is 1.0/m**3 , such as 1000.0/m**3 '
      read(5,56) hcon
      write(*,*) ' input add=? (f8.2) add*qc or add*r; ref:1.0 '
      read(5,56) add
      kbx=0
	kgwth=1
 55   format(i5)
 505  format(2i3)
 56   format(f8.2)
 57   format(f8.5)
      write(6,11)itex,wlfa,vlfa,hd,kgwth,icp,add,hcon
      write(66,11)itex,wlfa,vlfa,hd,kgwth,icp,add,hcon
   11 format(2x,' itex=',i3,' wlfa=',f4.1,' vlfa=',f4.1,
     1' hd=',f5.2,' kgrowth=',i4,' icp=',i4/,' add=',f6.1,' hcon=',f8.2)
      write(6,12) alg,a,icr,clfa,rlfa,nf
   12 format(2x,' alg=',f4.1,' a=',f4.1,' icr=',i4,' clfa=',f4.1,
     1' rlfa=',f4.1,' nf=',i3)
c -----------------------------------------------------------------------------
	kmick=0
 9000 FORMAT(40E15.8)
 9009 FORMAT(40i4)
  930 FORMAT(2I4)
  931 FORMAT(F12.1)
      write(*,*) '  after read '
c   -------  get hydrometeor feild: qc - cloud water, r- rain
      do 80 k=1,kz
      do 81 j=1,ni
      do 82 i=1,ni
      r(i,j,k)=qr(i,j,k)
      scr(i,j,k)=qc(i,j,k)+r(i,j,k)
   82 continue
   81 continue
   80 continue
c ------  get moisture feild: 
      do 90 k=1,kz
      do 91 j=1,ni
      do 92 i=1,ni
      qb(I,j,k)=qv(I,j,k)
      if(qc(i,j,k).gt.1.0e-6) qv(i,j,k)=qv(i,j,k)
      if(qc(i,j,k).le.1.0e-6) qv(i,j,k)=0.5*qv(i,j,k)
   92 continue
   91 continue
   90 continue
c -------------------------------
      do k=1,kz
      do j=1,jy
      do i=1,ix
      qc(i,j,k)=qc(i,j,k)*add
      r(i,j,k)=r(i,j,k)*add
      enddo
      enddo
      enddo
c  -------------------------------
      write(*,*) '  after 2222 '
c   ����ʾ������Ⱥ���λ��
      nk=0
      il=2
      ih=119
      do 110 k=2,24
      do 111 j=il,ih
      do 112 i=il,ih
      nk=nk+1
      x(nk)=i*dlx
      y(nk)=j*dly
      z(nk)=k*dlz
      xn=x(nk)/dlx
      yn=y(nk)/dly
      zn=z(nk)/dlz
      write(77,888) nk,xn,yn,zn
  112 continue
  111 continue
  110 continue
888   format(i6,3f8.3)
      write(*,*) '  after 110 '
      do 100 i=1,nk
      d(i)=hd
      ali(i)=alg
	con(i)=hcon
	kbig(i)=0
      m(i)=pai*ali(i)*d(i)**3/6.0
      call avag(x,y,z,p0,apo,i)
      call avag(x,y,z,tr,atr,i)
      call avag(x,y,z,r,ar,i)
	if(atr.lt.258.0) kbig(i)=1
      alo=apo/(2870.0*atr)
      vh(i)=sqrt(4.0*980.0*d(i)*ali(i)/(3.0*cd*alo))
      vh(i)=vh(i)*vlfa
      vr(i)=0.0
      if(ar.gt.1.0e-7)
     1 vr(i)=6287.96*(alo*ar/(pai*0.08))**0.2
  100 continue
      write(*,*) '  after 100 & do '
      nn=nk
      xm=(ix-1)*dlx
      ym=(jy-1)*dly
      zm=(kz-1)*dlz
      amc=0.0
      amr=0.0
      do k=1,kz
      do j=1,jy
      do i=1,ix
      amc=amax1(amc,qc(i,j,k))
      amr=amax1(amr,r(i,j,k))
      enddo
      enddo
      enddo
      write(*,*) ' amc=',amc,' amr=',amr
      write(*,*) '  before loop '
      do k=1,kz-1
      do j=1,jy-1
      do i=1,ix-1
      aqc(i,j,k)=(qc(i,j,k)+qc(i+1,j,k)+qc(i,j+1,k)+qc(i+1,j+1,k)+
     1    qc(i,j,k+1)+qc(i+1,j,k+1)+qc(i,j+1,k+1)+qc(i+1,j+1,k+1))*
     2  (0.5*dlx)*(0.5*dly)*(0.5*dlz)
      aqr(i,j,k)=(r(i,j,k)+r(i+1,j,k)+r(i,j+1,k)+r(i+1,j+1,k)+
     1    r(i,j,k+1)+r(i+1,j,k+1)+r(i,j+1,k+1)+r(i+1,j+1,k+1))*
     2  (0.5*dlx)*(0.5*dly)*(0.5*dlz)
      enddo
      enddo
      enddo
c ------------------------------------------------
      iii=0
      nstep=0
	nout=0
  900 continue
      nstep=nstep+1
      iii=iii+1
      ni=(iii/15)*15
      nj=(iii/nuw)*nuw
      do 1000 i=1,nk
      if(kh(i).eq.1) go to 300
      if(x(i).lt.dlx.or.x(i).gt.xm.or.y(i).lt.dly.or.y(i).gt.ym.
     1 or.z(i).lt.dlz.or.z(i).gt.zm) go to 201
      call avag(x,y,z,u,au,i)
      call avag(x,y,z,v,av,i)
      call avag(x,y,z,w,aw,i)
      x(i)=x(i)+au*dlt
      y(i)=y(i)+av*dlt
      z(i)=z(i)+(aw-vh(i))*dlt
c ----------------------------------
      if(kbx.eq.1) then
      xin=18.0*dlx
      xma=23.0*dlx
      yin=18.0*dlx
      yma=23.0*dlx
      zin=6.0*dlz
      zma=15.0*dlz
      if(x(i).gt.xin.and.x(i).lt.xma)then
      if(y(i).gt.yin.and.y(i).lt.yma)then
      if(z(i).gt.zin.and.z(i).lt.zma)then
      if(d(i).gt.0.08) d(i)=0.08
      endif
      endif
      endif
       endif
c ----------------------------------
      if(x(i).lt.dlx.or.x(i).gt.xm.or.y(i).lt.dly.or.y(i).gt.ym.
     1 or.z(i).lt.dlz.or.z(i).gt.zm) go to 201
      go to 200
  201 kh(i)=1
      if(x(i).lt.dlx) x(i)=dlx
      if(y(i).lt.dly) y(i)=dly
      if(z(i).lt.dlz) z(i)=dlz
      if(x(i).gt.xm) x(i)=xm
      if(y(i).gt.ym) y(i)=ym
      if(z(i).gt.zm) z(i)=zm
  200 if(kh(i).eq.1) go to 300
      call avag(x,y,z,r,ar,i)
      call avag(x,y,z,qc,ac,i)
      call avag(x,y,z,p0,apo,i)
      call avag(x,y,z,tr,atr,i)
      call avag(x,y,z,qv,aq,i)
      ps=9.31e-3*(atr**1.8)/apo
	if(atr.lt.258.0.and.kbig(i).eq.0) kbig(i)=1
	if(atr.gt.273.0.and.kbig(i).eq.1) kbig(i)=0
      alo=apo/(2870.0*atr)
      bnew=1.4962e-5*(1.0/(atr+120.0))*atr**1.5/alo
      bka=0.33372*bnew*alo
      pc2h=0.25*pai*d(i)*d(i)*vh(i)*ac*alo
      pr2h=0.25*pai*d(i)*d(i)*abs(vh(i)-vr(i))*ar*alo
      phwet=-2.0*pai*d(i)*(bka*(atr-273.0)+597.3*ps*alo
     1*(aq-3.799/apo))*(1.6+0.3*sqrt(d(i)*vh(i)/bnew))/alf
      if(atr.lt.273.16) go to 303
  707 continue
      dlm=pc2h+pr2h
	m(i)=m(i)+dlm*dlt
	if(ali(i).le.1.0) then
	ali(i)=1.0
	d(i)=(6.0*m(i)/(pai*ali(i)))**(1.0/3.0)
	endif
      go to 400
  304 m(i)=m(i)+dlm*dlt
      if(m(i).le.ecr) m(i)=ecr
      ord=d(i)
      d(i)=(6.0*m(i)/(pai*alio))**(1.0/3.0)
      dld=(d(i)-ord)/dlt
      ali(i)=alio
      go to 400
  303 continue
      if(icp.eq.2) sp=pc2h+pr2h
      if(icp.eq.1) sp=pc2h
      if(phwet.lt.sp) go to 301
      go to 500
  301 continue
      if(itex.eq.1) go to 309
      dlm=phwet
      if(kgwth.eq.0) dlm=0.0
      all=aln
      dld=2.0*dlm/(pai*all*d(i)*d(i))
      m(i)=m(i)+dlm*dlt
      d(i)=d(i)+dld*dlt
      ali(i)=6.0*m(i)/(pai*d(i)**3)
      go to 400
  309 continue
      if(ali(i).ge.alio) go to 310
      dlm=sp
      if(kgwth.eq.0) dlm=0.0
      m(i)=m(i)+dlm*dlt
      if(m(i).le.ecr) m(i)=ecr
      ali(i)=6.0*m(i)/(pai*d(i)**3)
      if(ali(i).lt.alio) go to 38
      ord=d(i)
      ali(i)=alio
      d(i)=(6.0*m(i)/(pai*alio))**(1.0/3.0)
      dld=(d(i)-ord)/dlt
      go to 400
   38 dld=0.0
      d(i)=d(i)
      go to 400
  310 continue
      dlm=phwet
      if(kgwth.eq.0) dlm=0.0
      m(i)=m(i)+dlm*dlt
      ali(i)=alio
      if(m(i).le.ecr) m(i)=ecr
      ord=d(i)
      d(i)=(6.0*m(I)/(pai*alio))**(1.0/3.0)
      dld=(d(i)-ord)/dlt
      go to 400
  500 continue
      dlm=sp
      if(kgwth.eq.0) dlm=0.0
      if(itex.eq.2) alm=aln
      if(itex.eq.1) alm=alg
      if(atr.lt.273.0) alm=0.11*(-a*vh(i)/(atr-273.16)/100.0)**0.76
      if(alm.ge.alio) alm=alio
      dld=2*dlm/(pai*alm*d(i)*d(i))
      d(i)=d(i)+dld*dlt
      m(i)=m(i)+dlm*dlt
      if(d(i).le.ecr) d(i)=ecr
      if(m(i).le.ecr) m(i)=ecr
      ali(i)=6.0*m(i)/(pai*d(i)**3)
      if(ali(i).lt.alio) go to 400
      ord=d(i)
      d(i)=(6.0*m(i)/(pai*alio))**(1.0/3.0)
      d(i)=d(i)+(d(i)-ord)/dlt
      ali(i)=alio
  400 continue
      vr(i)=0.0
      if(d(i).le.ecr) d(i)=ecr
      if (ali(i).le.aln) ali(i)=aln
      if(m(i).le.ecr) m(i)=ecr
      if(ar.ge.ecr) vr(i)=6287.96*(alo*ar/(pai*0.08))**0.2
      vh(i)=sqrt(4.0*980.0*d(i)*ali(i)/(3.0*cd*alo))
      vh(i)=vh(i)*vlfa
      if(atr.gt.273.0) go to 405
      if(jd(i).eq.1) go to 401
      go to 402
  401 if(phwet.ge.sp) go to 403
      jd(i)=0
      go to 404
  403 continue
      go to 405
  402 if(phwet.lt.sp) go to 404
      jd(i)=1
  404 continue
  405 continue
      nr=abs(nj-iii)
      nq=abs(ni-iii)
      if(nr.gt.0.1) go to 700
      if(jpr(i).eq.0.and.d(i).gt.dsh) then
      xi=x(i)/dlx
      yi=y(i)/dly
      zi=z(i)/dlz
      ii=ifix(x(i)/dlx)
      jj=ifix(y(i)/dly)
      kk=ifix(z(i)/dlz)
      tim=iii*dlt/60.0
      write(6,600) i,xi,yi,zi,ali(i),d(i),vh(i),m(i),tim,
     1 tr(ii,jj,kk),nstep,con(i),kbig(i)
      write(66,600) i,xi,yi,zi,ali(i),d(i),vh(i),m(i),tim,
     1 tr(ii,jj,kk),nstep,con(i),kbig(i)
  600 format(1x,' i=',i5,' x=',f4.1,' y=',f4.1,' z=',f4.1,' ali=',
     1 f4.2,' d=',f6.2,' vh=',f6.1,' m=',e8.2/,' t=',f6.1,' tr=',
     2 f5.1,' nstep=',i6,' con=',f6.1,' kbig=',i4)
      endif
      go to 700
  300 continue
c-----
      if(jpr(i).eq.1) go to 700
	if(jpr(i).eq.0) then
      if(z(i).le.dlz)            then
      write(*,*) ' z(i)=',z(i),' i=',i,' con=',con(i),' m=',m(i)
      write(66,*) ' z(i)=',z(i),' i=',i,' con=',con(i),' m=',m(i)
                                endif
      tim=iii*dlt/60.0
      x1=x(i)/dlx
      y1=y(i)/dly
      z1=z(i)/dlz
	nout=nout+1
      write(6,77) i,tim,d(i),x1,y1,z1,nstep,nout
      write(66,77) i,tim,d(i),x1,y1,z1,nstep,nout
   77 format(5x,'The I=',i6,' was out off the domain of conputation',
     1'  time(min)=',f6.2/,'    d=',e10.4,' xi=',f6.1,' yi=',f6.1,
     2 ' zi=',f6.1,' nstep=',i6,' nout=',i5)
      jpr(i)=1
   66 format(1x,'i=',i3,6e10.3)
                     endif
  700 continue
      if(kbig(i).eq.0) then
	if(d(i).gt.dbk) then
	d(i)=d(i)/1.25992
	m(i)=0.5*m(i)
	con(i)=2.0*con(i)
	vh(i)=sqrt(4.0*980.0*d(i)*ali(i)/(3.0*cd*alo))
      vh(i)=vh(i)*vlfa
	endif
	endif
 1000 continue
       nn=0
      do 950 i=1,nk
      if(kh(i).eq.1) nn=nn+1
      if(nn.gt.nk.or.nstep.gt.mm.or.nstep.gt.mstep.or.nout.ge.nk) 
	1  go to 5000
  950 continue
      mpr=(nstep/10)*10
      if(mpr.eq.nstep) then
      write(99) nstep,tim
      write(99) x,y,z,d,ali
      endif
c -----------------------------------------------------------------------------
      nsh=(nstep/nuw)*nuw
c ----------------------------
      go to 900
 5000 continue
      write(*,*) ' nk=',nk,' nn=',nn,' nstep=',nstep,' mstep=',mstep
      write(*,*) ' nk=',nk,' nout=',nout
	  return
      end
c ---------------------
      subroutine avag(x,y,z,a,avg,i)
      parameter (ix=120,jy=120,kz=25,ixx=ix*jy*kz,mm=1500)  
      common /a1/dlx,dly,dlz,dlt
      dimension x(ixx),y(ixx),z(ixx),a(ix,jy,kz)
      nx=ifix(x(i)/dlx)
      ny=ifix(y(i)/dly)
      nz=ifix(z(i)/dlz)
      ddx=x(i)-nx*dlx
      ddy=y(i)-ny*dly
      ddz=z(i)-nz*dlz
      nx1=nx+1
      ny1=ny+1
      nz1=nz+1
c
      ax11=a(nx,ny,nz)+(a(nx1,ny,nz)-a(nx,ny,nz))*ddx/dlx
      ax12=a(nx,ny1,nz)+(a(nx1,ny1,nz)-a(nx,ny1,nz))*ddx/dlx
      ay1=ax11+(ax12-ax11)*ddy/dly
c
      ax21=a(nx,ny,nz1)+(a(nx1,ny,nz1)-a(nx,ny,nz1))*ddx/dlx
      ax22=a(nx,ny1,nz1)+(a(nx1,ny1,nz1)-a(nx,ny1,nz1))*ddx/dlx
      ay2=ax21+(ax22-ax21)*ddy/dly
c
      avg=ay1+(ay2-ay1)*ddz/dlz     
c
      return
      end
c ---------------------

