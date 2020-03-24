    program alptraj
    parameter (ix=120,jy=120,kz=25,ixx=ix*jy*kz,mm=1500) 
!定义空间三维变量温度、气压、水汽
    dimension tr(ix,jy,kz),p0(ix,jy,kz),qv(ix,jy,kz)
!定义空间三维变量风场U、V、W
    dimension v(ix,jy,kz),u(ix,jy,kz),w(ix,jy,kz)
! 定义一维变量垂直方向空气密度
    dimension al(kz)
C  -----------------------------------------------------------
! 定义云水、雨水比含水量
    dimension qc(ix,jy,kz),qr(ix,jy,kz)
! 定义一维垂直方向温度，气压
    dimension treo(kz),po(kz)
! 定义全局变量水平、垂直网格间距、时间间隔
    common /a1/dlx,dly,dlz,dlt
! 定义一维垂直方向温度，气压，密度 
    common /a3/treo,po,al
! 定义时间步数、示踪粒子数目
    common /a4/iii,nn
! 定义控制变量ntex控制nx取值，itex控制湿增长分支、jtex、ktex、mtexb未用到
    common /a9/ntex,itex,jtex,ktex,mtexb
! 输入时间间隔，网格距离
    write(*,*) '  input dlt,dlx,dly,dlz=? '
    write(*,*) 'such as dlt=5.,dlx=50000.,dly=50000.,dlz=50000.,[cm]'
    read(*,*) dlt,dlx,dly,dlz
    write(*,*) '  input from model-output'
! 输入数据
    open(8,file='modelin.dat',form='unformatted',status='unknown') 
! 数据格式整型数
 55   format(i5)
! 读取外部数据风场、温度、气压、水汽、云水雨水
! 只读取一组数据
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
! 将外部读取的气象场数据导入函数中
      call xytraj0(u,v,w,tr,qv,qc,qr,p0)
	stop
	end
c ----------------------------------------------------------------------------- 
! 子函数定义
      subroutine xytraj0(u,v,w,tr,qv,qc,qr,p0)
! 定义冰雹质量为实数
      real m
! 外网格数为120*120*25，总网格数为ixx，共1500个时刻
      parameter (ix=120,jy=120,kz=25,ixx=ix*jy*kz,mm=1500)  
! 将所有冰雹轨迹点的坐标、密度、直径都存储在矩阵中
      dimension x(ixx),y(ixx),z(ixx),ali(ixx),d(ixx)
! 雨滴下落末速度，数浓度
      dimension vr(ixx),vh(ixx),con(ixx)
! 冰雹粒子质量
      dimension m(ixx)
! kh用于标记即将跑出区域的粒子、kbig用于控制冰雹融化
      dimension kh(ixx),kbig(ixx)
! 定义干湿增长判断、jpr
      dimension jd(ixx),jpr(ixx)
! 定义一维垂直方向温度，气压，密度 
      dimension treo(kz),po(kz),al(kz)
! 三维数据温度场、气压场
      dimension tr(ix,jy,kz),p0(ix,jy,kz)
! 三维数据冰雹所在位置的云水比湿、雨水比湿、sh未用到变量
      dimension aqc(ix,jy,kz),aqr(ix,jy,kz),sh(ix,jy,kz)
! 三维数据网格点上的水汽、订正前的水汽场、r、水凝物场=云水+雨水
      dimension qv(ix,jy,kz),qb(ix,jy,kz),r(ix,jy,kz),scr(ix,jy,kz)
! 三维数据网格点上的风场U、V、W
      DIMENSION V(ix,jy,kz),U(ix,jy,kz),W(ix,jy,kz)           
C  -----------------------------------------------------------
! 三维数据网格点上的云水、雨水比含水量
      dimension qc(ix,jy,kz),qr(ix,jy,kz)
c ------
! 定义全局变量水平、垂直网格间距、时间间隔
      common /a1/dlx,dly,dlz,dlt
! 定义一维垂直方向温度，气压，密度 
      common /a3/treo,po,al
! 定义时间步数、示踪粒子数目
      common /a4/iii,nn
! 定义控制变量ntex控制nx取值，itex控制湿增长分支、jtex、ktex、mtexb未用到
      common /a9/ntex,itex,jtex,ktex,mtex
! 输出文件
      open(66,file='3rptxy.dat',status='unknown')
      open(99,file='3trajxy.dat',form='unformatted',status='unknown')
      open(77,file='xyz0.dat',status='unknown')
! 30个时间间隔为周期
	nuw=30
c  -------------------------------
      do 5 i=1,ixx
! 网格点坐标变量一维矩阵赋初值
      x(i)=0.0
      y(i)=0.0
      z(i)=0.0
! 默认初始时刻冰雹普遍为干增长
      jd(i)=1
! 
      jpr(i)=0
! 雨滴下落末速度初值赋为0
      vr(i)=0.0
! 冰雹下落末速度初值赋为0
      vh(i)=0.0
! 示踪雹胚初始密度初值赋为0
      ali(i)=0.0
! 示踪雹胚半径初值赋为0
      d(i)=0.0
! 示踪雹胚质量初值赋为0
      m(i)=0.0
! 默认初始时刻冰雹均未跑出边缘外
      kh(i)=0
   5  continue
c  -------------------------------
! 得到第一行第一列垂直方向所有网格点温度、气压数据
c  -------------------------------
      do 6 k=1,kz
      treo(k)=tr(1,1,k)
      po(k)=p0(1,1,k)
    6 continue
c  -------------------------------
! 网格点数据sh(未用到) 雨水比含水量赋予初值
! 云水含水量、雨水含水量赋初值
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
c  -------------------------------
! 定义阻力系数(coefficient of drag)
      cd=0.45
! 没有用到的变量
      hbtex=0
! itex控制湿增长的两个分支(A、B)
      itex=1
! 没有用到的变量
      jtex=1
! 没有用到的变量
      ktex=3
! 没有用到的变量
      mtex=1
! 
      ntex=1
! 碰并增长选择项目(是否考虑冰雹并冻结雨水率)
      icp=2
! 3rptxy.dat数据的头信息 
      icr=0
! 定义粒子群的初始密度为0.9
      alg=0.9
! 若密度是否大于等于0.9，icr=1
      if(alg.ge.0.9) icr=1
      a=30.0
! 没有用到的变量
      iwrq=0
! 冰雹直径单位是cm
c   the unit of hd is cm.
! 3rptxy.dat数据的头信息 
      nf=0
! 没有用到的变量
      als=0.6
! Lf冻结潜热
      alf=79.3
! 
      alm=0.5
! 干增长密度为0.3
      aln=0.3
! 下落末速度订正系数(0.6、1.0可选)
      vlfa=0.6
c      vlfa=1.0
! 3rptxy.dat数据的头信息 
      clfa=0.8
! 3rptxy.dat数据的头信息 
      rlfa=1.0
! 没有用到的变量
      ulfa=1.0
! 3rptxy.dat数据的头信息 
      wlfa=1.0
! 时间间隔为5s
      dlt=5.0
      pai=3.14159264
! 密度为0.9
      alio=0.9
! 去野点阈值(控制直径、质量不小于ecr)
      ecr=1.0e-6
! 没有用到的变量
      epcr=1.0e-5
! 网格间距
      dlx=50000.0
      dly=50000.0
      dlz=50000.0
! 时间步数初值赋为0，每随时间积分一次加一
      iii=0
! 普通的计数变量，应用于多个循环中
      nn=0
c ------------  get temp. , presure feild:
! 用状态方程通过气压场温度场推算空气密度场
      do k=1,kz
      al(k)=po(k)/(2840.0*treo(k))
      enddo
c ---------------------------
! ni为水平方向网格点数，如果ntex=1，ni=120
      ni=ix
      if(ntex.eq.1) ni=ix
c ----- readin macro & micro data of hail cloud feild -----------
! 读取冰雹云宏微观场数据
! 最大垂直风速场记载
c      write(*,*) ' input nfile=? (i5) which is records of wmax case =3'
c      read(5,55) nfile
! 输入时间步数
      write(*,*) ' input mstep=? (i5) mstep<= mm=1500'
      read(5,55) mstep
! 输入初始雹胚直径
	write(*,*) ' input hd=? (f8.2) ,such as 0.05 cm.'
      read(5,56) hd
! 输入直径bk
	write(*,*) ' input dbk=? (f8.2) ,such as 0.4-->1.0 cm.'
      read(5,56) dbk
      write(*,*) ' input dsh=? (f8.2) d>dsh to show information'
! 输入显示出终端的直径阈值
      read(5,56) dsh
! 输入雹胚数浓度，单位每立方米
      write(*,*) ' input hcon=? (f8.2) concentrtion of embryo ref:1.0 '
      write(*,*) '  the unit of hcon is 1.0/m**3 , such as 1000.0/m**3 '
      read(5,56) hcon
! 输入云水、雨水订正系数
      write(*,*) ' input add=? (f8.2) add*qc or add*r; ref:1.0 '
      read(5,56) add
! 控制是否使某区域内雹胚增长直径存在上限()
      kbx=0
! 冰雹质量增长率控制项
      kgwth=1
 55   format(i5)
 505  format(2i3)
 56   format(f8.2)
 57   format(f8.5)
! 将参数设置存储在3rptxy.dat第一行
      write(6,11)itex,wlfa,vlfa,hd,kgwth,icp,add,hcon
      write(66,11)itex,wlfa,vlfa,hd,kgwth,icp,add,hcon
   11 format(2x,' itex=',i3,' wlfa=',f4.1,' vlfa=',f4.1,
     1' hd=',f5.2,' kgrowth=',i4,' icp=',i4/,' add=',f6.1,' hcon=',f8.2)
      write(6,12) alg,a,icr,clfa,rlfa,nf
   12 format(2x,' alg=',f4.1,' a=',f4.1,' icr=',i4,' clfa=',f4.1,
     1' rlfa=',f4.1,' nf=',i3)
c -----------------------------------------------------------------------------
! 没有用到的变量
      kmick=0
 9000 FORMAT(40E15.8)
 9009 FORMAT(40i4)
  930 FORMAT(2I4)
  931 FORMAT(F12.1)
! 读取数据后
      write(*,*) '  after read '
c   -------  get hydrometeor feild: qc - cloud water, r- rain
! 得到水凝物场=云水+雨水
      do 80 k=1,kz
      do 81 j=1,ni
      do 82 i=1,ni
      r(i,j,k)=qr(i,j,k)
      scr(i,j,k)=qc(i,j,k)+r(i,j,k)
   82 continue
   81 continue
   80 continue
! 得到订正后的水汽场
! 订正方式，将小于10的-6次方的水汽场缩至原来一半
c ------  get moisture feild: 
      do 90 k=1,kz
      do 91 j=1,ni
      do 92 i=1,ni
      qb(i,j,k)=qv(i,j,k)
      if(qc(i,j,k).gt.1.0e-6) qv(i,j,k)=qv(i,j,k)
      if(qc(i,j,k).le.1.0e-6) qv(i,j,k)=0.5*qv(i,j,k)
   92 continue
   91 continue
   90 continue
c -------------------------------
! 云水、雨水乘以订正系数add
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
! 定义示踪冰雹轨迹点坐标初值位于网格点上
      nk=0
      il=2
      ih=119
      do 110 k=2,kz-1
      do 111 j=il,ih
      do 112 i=il,ih
      nk=nk+1
      x(nk)=i*dlx
      y(nk)=j*dly
      z(nk)=k*dlz
      xn=x(nk)/dlx
      yn=y(nk)/dly
      zn=z(nk)/dlz
! 将冰雹网格尺度轨迹信息存储于xyz0.dat中
      write(77,888) nk,xn,yn,zn
  112 continue
  111 continue
  110 continue
888   format(i6,3f8.3)
      write(*,*) '  after 110 '
c  -------------------------------
      do 100 i=1,nk
! 示踪雹胚初始直径0.05cm
      d(i)=hd
! 示踪雹胚初始密度0.9
      ali(i)=alg
! 示踪雹胚初始数浓度1000.0/m**3
     	con(i)=hcon
    	kbig(i)=0
! 示踪雹胚质量大小(冰密度0.9)
      m(i)=pai*ali(i)*d(i)**3/6.0
! 将网格点上的气压，温度，雨水比含水量插值到示踪粒子所在轨迹处
      call avag(x,y,z,p0,apo,i)
      call avag(x,y,z,tr,atr,i)
      call avag(x,y,z,r,ar,i)
! 如果示踪粒子处温度小于零下15度
	if(atr.lt.258.0) kbig(i)=1
! 由示踪粒子处的气压温度计算空气密度
      alo=apo/(2870.0*atr)
! 计算冰雹粒子的下落末速度
      vh(i)=sqrt(4.0*980.0*d(i)*ali(i)/(3.0*cd*alo))
! 下落末速度订正系数为0.6
      vh(i)=vh(i)*vlfa
! 雨滴下落末速度初始赋值为0
      vr(i)=0.0
! 如果雨水比含水量大于10的-7次方
      if(ar.gt.1.0e-7)
! 计算雨滴的下落末速度
     1 vr(i)=6287.96*(alo*ar/(pai*0.08))**0.2
  100 continue
c  -------------------------------
      write(*,*) '  after 100 & do '
! nn为示踪粒子数目
      nn=nk
! 冰雹轨迹的最大刚性边界
      xm=(ix-1)*dlx
      ym=(jy-1)*dly
      zm=(kz-1)*dlz
! 定义网格点上的最大云水比含水量，最大雨水比含水量并赋初值0.0
      amc=0.0
      amr=0.0
! 遍历所有网格点求云水、雨水比含水量最大值 
      do k=1,kz
      do j=1,jy
      do i=1,ix
      amc=amax1(amc,qc(i,j,k))
      amr=amax1(amr,r(i,j,k))
      enddo
      enddo
      enddo
! 将最大值输出至终端
      write(*,*) ' amc=',amc,' amr=',amr
      write(*,*) '  before loop '
      do k=1,kz-1
      do j=1,jy-1
      do i=1,ix-1
! 计算单位网格立方内的云水、雨水含水质量
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
c ------------------------------------------------
  900 continue
      nstep=nstep+1
      iii=iii+1
! 15的整数倍，时间间隔为75s
      ni=(iii/15)*15
! 30的整数倍，时间间隔为150s
      nj=(iii/nuw)*nuw
c ------------------------------------------------
      do 1000 i=1,nk
! 跑到边缘的粒子不参与循环计算，跳步至300
      if(kh(i).eq.1) go to 300
      if(x(i).lt.dlx.or.x(i).gt.xm.or.y(i).lt.dly.or.y(i).gt.ym.
     1 or.z(i).lt.dlz.or.z(i).gt.zm) go to 201
! 将网格点上的风场U、V、W插值到示踪粒子所在轨迹处
      call avag(x,y,z,u,au,i)
      call avag(x,y,z,v,av,i)
      call avag(x,y,z,w,aw,i)
! 通过风场、下落末速度更新示踪粒子的位置信息
      x(i)=x(i)+au*dlt
      y(i)=y(i)+av*dlt
      z(i)=z(i)+(aw-vh(i))*dlt
c ----------------------------------
! kbx=0 设置示踪粒子初始区域内某一子区域内雹胚增长上限为0.08cm
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
! 当粒子即将跑出区域外前往201，否则前往200
      if(x(i).lt.dlx.or.x(i).gt.xm.or.y(i).lt.dly.or.y(i).gt.ym.
     1 or.z(i).lt.dlz.or.z(i).gt.zm) go to 201
      go to 200
c ----------------------------------
! 将即将跑出区域的粒子拉回来，并通过kh(i)=1做标记
  201 kh(i)=1
      if(x(i).lt.dlx) x(i)=dlx
      if(y(i).lt.dly) y(i)=dly
      if(z(i).lt.dlz) z(i)=dlz
      if(x(i).gt.xm) x(i)=xm
      if(y(i).gt.ym) y(i)=ym
      if(z(i).gt.zm) z(i)=zm
c ----------------------------------
! 跑到边缘的粒子不再参与循环计算，跳过接下来步骤，前往300
  200 if(kh(i).eq.1) go to 300
! 对于没有跑出边缘的粒子进行如下步骤
! 将网格点上的雨水、云水比含水量，气压、温度、水汽场插值到示踪粒子所在轨迹处
      call avag(x,y,z,r,ar,i)
      call avag(x,y,z,qc,ac,i)
      call avag(x,y,z,p0,apo,i)
      call avag(x,y,z,tr,atr,i)
      call avag(x,y,z,qv,aq,i)
! 计算水汽扩散系数
      ps=9.31e-3*(atr**1.8)/apo
! 如果上一时刻温度并大于等于零下15度，当前示踪粒子温度小于零下15度，kbig(i)=1
	if(atr.lt.258.0.and.kbig(i).eq.0) kbig(i)=1
! 如果上一时刻温度小于等于零度，当前示踪粒子温度大于零度，kbig(i)=0，考虑融化
	if(atr.gt.273.0.and.kbig(i).eq.1) kbig(i)=0
! 由示踪粒子处的气压温度计算该时刻空气密度
      alo=apo/(2870.0*atr)
! 计算？？？
      bnew=1.4962e-5*(1.0/(atr+120.0))*atr**1.5/alo
! 计算热传导系数
      bka=0.33372*bnew*alo
! 冰雹并冻云水率
      pc2h=0.25*pai*d(i)*d(i)*vh(i)*ac*alo
! 冰雹并冻雨水率
      pr2h=0.25*pai*d(i)*d(i)*abs(vh(i)-vr(i))*ar*alo
! 计算热收支率
      phwet=-2.0*pai*d(i)*(bka*(atr-273.0)+597.3*ps*alo
     1*(aq-3.799/apo))*(1.6+0.3*sqrt(d(i)*vh(i)/bnew))/alf
! 如果温度温度大于0.16度前往303
       if(atr.lt.273.16) go to 303
c ----------------------------------
! 如果温度温度小于等于0.16度
  707 continue
! 单位时间冰雹的质量增长率=冰雹并冻云水率+冰雹并冻雨水率
      dlm=pc2h+pr2h
! 计算该时刻冰雹质量
      m(i)=m(i)+dlm*dlt
! 冰雹密度小于1(0.9)，转换为水的密度1
      if(ali(i).le.1.0) then
      ali(i)=1.0
! 冰雹全部融化后的直径
      d(i)=(6.0*m(i)/(pai*ali(i)))**(1.0/3.0)
      endif
      go to 400
c ----------------------------------
! 没有指向该段 ,应该属于309湿增长草稿
  304 m(i)=m(i)+dlm*dlt
! 控制质量不小于1.0e-6
      if(m(i).le.ecr) m(i)=ecr
      ord=d(i)
      d(i)=(6.0*m(i)/(pai*alio))**(1.0/3.0)
      dld=(d(i)-ord)/dlt
      ali(i)=alio
      go to 400
c ----------------------------------
! 当温度温度大于0.16度时
  303 continue
! icp=2，既考虑冰雹并冻云水率，也考虑冰雹并冻雨率
! icp=1，只考虑冰雹并冻云水率
! icp=2
      if(icp.eq.2) sp=pc2h+pr2h
      if(icp.eq.1) sp=pc2h
! 如果热收支率，小于冰雹碰并增长率，为湿度增长，否则为干增长 
      if(phwet.lt.sp) go to 301
      go to 500
c ----------------------------------
! 湿增长 (itex=1) itex控制湿增长的两个分支(A、B)
  301 continue
! 前往309，湿增长分支A
      if(itex.eq.1) go to 309
c ----------------------------------
! 如果itex不等于1，实际并未执行湿增长分支B
! 单位时间冰雹的质量增长率=热收支率
      dlm=phwet
! 如果kgwth=0，冰雹质量率为0.0，不增长(kgwth=1)
      if(kgwth.eq.0) dlm=0.0
! 如果kgwth=1，冰雹质量增长率为热收支率(kgwth=1)
! all=0.3
      all=aln
! 冰雹直径增长率
      dld=2.0*dlm/(pai*all*d(i)*d(i))
! 计算该时刻冰雹质量
      m(i)=m(i)+dlm*dlt
! 计算该时刻直径
      d(i)=d(i)+dld*dlt
! 通过该直径和质量推算出冰雹的密度
      ali(i)=6.0*m(i)/(pai*d(i)**3)
! 前往400
      go to 400
c ----------------------------------
! 湿增长
  309 continue
! 如果冰雹密度大于等于0.9，前往310，湿增长分支A下的分支1
      if(ali(i).ge.alio) go to 310
! 湿增长分支A下的分支2
! 质量增长仅仅为碰并增长
      dlm=sp
! 如果kgwth=0，冰雹质量率为0.0，不增长(kgwth=1)
      if(kgwth.eq.0) dlm=0.0
! 如果kgwth=1，冰雹质量增长率为热收支率(kgwth=1)
! 计算该时刻冰雹质量
      m(i)=m(i)+dlm*dlt
! 控制质量不小于1.0e-6
      if(m(i).le.ecr) m(i)=ecr
! 通过该直径和质量推算出冰雹的密度
      ali(i)=6.0*m(i)/(pai*d(i)**3)
! 如果冰雹密度比0.9还小，存在气泡
      if(ali(i).lt.alio) go to 38
! 如果冰雹密度大于等于0.9
! ord表示上一时刻冰雹直径
      ord=d(i)
! 密度从新调整为0.9
      ali(i)=alio
! 通过该密度和质量推算出该时刻冰雹的直径
      d(i)=(6.0*m(i)/(pai*alio))**(1.0/3.0)
! 由该时刻冰雹直径和上一时刻冰雹直径推算出冰雹直径增长率
      dld=(d(i)-ord)/dlt
! 前往400
      go to 400
c ----------------------------------
! 如果冰雹密度比0.9还小，存在气泡
! 冰雹直径不增加
   38 dld=0.0
      d(i)=d(i)
      go to 400
c ----------------------------------
! 湿增长分支A下的分支1，冰雹密度大于等于0.9的情况
  310 continue
! 单位时间冰雹的质量增长率=热收支率
      dlm=phwet
! 如果kgwth=0，冰雹质量率为0.0，不增长(kgwth=1)
      if(kgwth.eq.0) dlm=0.0
! 如果kgwth=1，冰雹质量增长率为热收支率(kgwth=1)
! 计算该时刻冰雹质量
      m(i)=m(i)+dlm*dlt
! 密度从新调整为0.9
      ali(i)=alio
! 控制质量不小于1.0e-6
      if(m(i).le.ecr) m(i)=ecr
! ord表示上一时刻冰雹直径
      ord=d(i)
! 通过该密度和质量推算出该时刻冰雹的直径
      d(i)=(6.0*m(i)/(pai*alio))**(1.0/3.0)
! 由该时刻冰雹直径和上一时刻冰雹直径推算出冰雹直径增长率
      dld=(d(i)-ord)/dlt
      go to 400
c ----------------------------------
! 干增长
  500 continue
! 质量增长仅仅为碰并增长
      dlm=sp
! 如果kgwth=0，冰雹质量率为0.0，不增长(kgwth=1)
      if(kgwth.eq.0) dlm=0.0
! itex=2，alm=0.3(alm最初值为0.5)
      if(itex.eq.2) alm=aln
! itex=1，alm=0.9(alm最初值为0.5)
      if(itex.eq.1) alm=alg
! 如果温度小于零度，alm
      if(atr.lt.273.0) alm=0.11*(-a*vh(i)/(atr-273.16)/100.0)**0.76
! 确保alm小于等于0.9
      if(alm.ge.alio) alm=alio
! 冰雹直径增长率
      dld=2*dlm/(pai*alm*d(i)*d(i))
! 计算该时刻直径
      d(i)=d(i)+dld*dlt
! 计算该时刻冰雹质量
      m(i)=m(i)+dlm*dlt
! 控制直径不小于1.0e-6
      if(d(i).le.ecr) d(i)=ecr
! 控制质量不小于1.0e-6
      if(m(i).le.ecr) m(i)=ecr
! 通过该直径和质量推算出该时刻冰雹的密度
      ali(i)=6.0*m(i)/(pai*d(i)**3)
! 如果密度小于0.9，前往400
      if(ali(i).lt.alio) go to 400
! 如果密度大于等于0.9
! ord表示上一时刻冰雹直径
      ord=d(i)
! 通过该密度和质量推算出该时刻冰雹的直径
      d(i)=(6.0*m(i)/(pai*alio))**(1.0/3.0)
! 由该时刻冰雹直径和上一时刻冰雹直径推算出下一时刻冰雹直径
      d(i)=d(i)+(d(i)-ord)/dlt
! 密度从新调整为0.9
      ali(i)=alio
c ----------------------------------
  400 continue
! 雨滴下落末速度重新赋值为0
      vr(i)=0.0
! 控制直径不小于1.0e-6
      if(d(i).le.ecr) d(i)=ecr
! 控制冰雹密度大于等于0.3
      if (ali(i).le.aln) ali(i)=aln
! 控制直径不小于1.0e-6
      if(m(i).le.ecr) m(i)=ecr
! 如果雨滴下落末速度大于1.0e-6，重新计算雨滴的下落末速度
      if(ar.ge.ecr) vr(i)=6287.96*(alo*ar/(pai*0.08))**0.2
! 粒子下落末速度往往与其尺度有关
      vh(i)=sqrt(4.0*980.0*d(i)*ali(i)/(3.0*cd*alo))
! 下落末速度订正系数为0.6
      vh(i)=vh(i)*vlfa
! 如果温度温度大于零度，前往405
      if(atr.gt.273.0) go to 405
c ----------------------------------
! 如果温度温度小于等于零度
! 判断上一时刻冰雹为何种增长，干增长前往401，湿增长前往402
      if(jd(i).eq.1) go to 401
      go to 402
c ----------------------------------
! 上一时刻为干增长，若该时刻为湿增长jd(i)=0
  401 if(phwet.ge.sp) go to 403
      jd(i)=0
      go to 404
c ----------------------------------
! 干增长则保持不变
  403 continue
      go to 405
c ----------------------------------
! 上一时刻为湿增长，若该时刻为干增长jd(i)=1
  402 if(phwet.lt.sp) go to 404
      jd(i)=1
c ----------------------------------
! 湿增长则保持不变
  404 continue
c ----------------------------------
! 输出3rptxy.dat模块
  405 continue
! 0~14时刻交替更新
      nr=abs(nj-iii)
! 0~29时刻交替更新
      nq=abs(ni-iii)
! 非零项前往700，即每隔15个时刻(75s)
c ----------------------------------
      if(nr.gt.0.1) go to 700
! jpr(i)=0 冰雹直径大于显示出在终端的直径阈值
c ----------------------------------
      if(jpr(i).eq.0.and.d(i).gt.dsh) then
! 得到整数和小数的网格点坐标信息
      xi=x(i)/dlx
      yi=y(i)/dly
      zi=z(i)/dlz
      ii=ifix(x(i)/dlx)
      jj=ifix(y(i)/dly)
      kk=ifix(z(i)/dlz)
! 模式模拟时间(单位分钟)
      tim=iii*dlt/60.0
! 写入内容输出到终端上
      write(6,600) i,xi,yi,zi,ali(i),d(i),vh(i),m(i),tim,
     1 tr(ii,jj,kk),nstep,con(i),kbig(i)
! 写入示踪粒子序号，网格坐标，密度，半径，冰雹下落末速度，冰雹质量，模拟时间(分钟)，
! 网格点上的温度，时间步数，数浓度，kbig
      write(66,600) i,xi,yi,zi,ali(i),d(i),vh(i),m(i),tim,
     1 tr(ii,jj,kk),nstep,con(i),kbig(i)
  600 format(1x,' i=',i5,' x=',f4.1,' y=',f4.1,' z=',f4.1,' ali=',
     1 f4.2,' d=',f6.2,' vh=',f6.1,' m=',e8.2/,' t=',f6.1,' tr=',
     2 f5.1,' nstep=',i6,' con=',f6.1,' kbig=',i4)
      endif
c ----------------------------------
      go to 700
c ----------------------------------
! 跑到边缘的粒子不参与循环计算，跳步至此
  300 continue
! 如果jpr=1，前往700(jpr(i)=0)
c ----------------------------------
      if(jpr(i).eq.1) go to 700
! 如果jpr=0(jpr(i)=0)
c ----------------------------------
      if(jpr(i).eq.0) then
! 筛选高度降落至地表的冰雹
c ----------------------------------
      if(z(i).le.dlz) then
! 把高度、示踪粒子序号、数浓度、冰雹质量输出到终端
      write(*,*) ' z(i)=',z(i),' i=',i,' con=',con(i),' m=',m(i)
! 把高度、示踪粒子序号、数浓度、冰雹质量写入3rptxy.dat
      write(66,*) ' z(i)=',z(i),' i=',i,' con=',con(i),' m=',m(i)
      endif
c ----------------------------------
! 模式模拟时间(单位分钟)
      tim=iii*dlt/60.0
      x1=x(i)/dlx
      y1=y(i)/dly
      z1=z(i)/dlz
! 计数跑出区域外的粒子数
      nout=nout+1
      write(6,77) i,tim,d(i),x1,y1,z1,nstep,nout
      write(66,77) i,tim,d(i),x1,y1,z1,nstep,nout
   77 format(5x,'The I=',i6,' was out off the domain of conputation',
     1'  time(min)=',f6.2/,'    d=',e10.4,' xi=',f6.1,' yi=',f6.1,
     2 ' zi=',f6.1,' nstep=',i6,' nout=',i5)
      jpr(i)=1
   66 format(1x,'i=',i3,6e10.3)
      endif
c ----------------------------------
  700 continue
! 由于温度大于零度考虑冰雹融化
      if(kbig(i).eq.0) then
! 如果直径大于dbk
      if(d(i).gt.dbk) then
! 直径除以1.25992
      d(i)=d(i)/1.25992
! 质量缩小为原来一半
      m(i)=0.5*m(i)
! 数浓度增加为原来两倍
      con(i)=2.0*con(i)
! 冰雹下落末速度 
      vh(i)=sqrt(4.0*980.0*d(i)*ali(i)/(3.0*cd*alo))
! 下落末速度订正系数为0.6
      vh(i)=vh(i)*vlfa
      endif
      endif
c ----------------------------------
 1000 continue
! 定义初始化用于统计边缘粒子数目
       nn=0
c ----------------------------------
      do 950 i=1,nk
! 粒子跑到边缘(kh(i)=1)，nn计数加一
      if(kh(i).eq.1) nn=nn+1
! 跑出边缘的粒子大于等于所有示踪粒子，模拟步数大于最大时间步数(1500)
! 模拟步数大于自定义变量(mstep<=mm=1500)，跑出边缘的粒子大于等于所有示踪粒子
      if(nn.gt.nk.or.nstep.gt.mm.or.nstep.gt.mstep.or.nout.ge.nk) 
! 跳出900时间当循环，前往5000
	1  go to 5000
  950 continue
c ----------------------------------
! 3trajxy.dat文件写入
! mpr为10的整数倍
      mpr=(nstep/10)*10
! if循环控制每隔50秒写入一次
      if(mpr.eq.nstep) then
! 写入总时间步数，模式模拟时间(分钟)
      write(99) nstep,tim
! 写入示踪粒子网格坐标半径，密度
      write(99) x,y,z,d,ali
      endif
c -----------------------------------------------------------------------------
      nsh=(nstep/nuw)*nuw
c ----------------------------
      go to 900
c -----------------------------------------------------------------------------
 5000 continue
      write(*,*) ' nk=',nk,' nn=',nn,' nstep=',nstep,' mstep=',mstep
      write(*,*) ' nk=',nk,' nout=',nout
	  return
      end
c ---------------------
! 将网格点上的物理场插值到示踪粒子所在轨迹处
      subroutine avag(x,y,z,a,avg,i)
! 外网格数为120*120*25，总网格数为ixx，共1500个时刻
      parameter (ix=120,jy=120,kz=25,ixx=ix*jy*kz,mm=1500)  
      common /a1/dlx,dly,dlz,dlt
! 将所有冰雹轨迹点的坐标、以及需要插值的物理场a(温压湿风)
      dimension x(ixx),y(ixx),z(ixx),a(ix,jy,kz)
! ifix去除小数位取整
! 获得离轨迹点较近的网格点坐标nx、ny、nz
      nx=ifix(x(i)/dlx)
      ny=ifix(y(i)/dly)
      nz=ifix(z(i)/dlz)
! 轨迹点距网格点相差的距离:dx、dy、dz
      ddx=x(i)-nx*dlx
      ddy=y(i)-ny*dly
      ddz=z(i)-nz*dlz
! 轨迹点x坐标位于网格点nx和nx1之间
      nx1=nx+1
! 轨迹点y坐标位于网格点ny和ny1之间
      ny1=ny+1
! 轨迹点z坐标位于网格点nz和nz1之间
      nz1=nz+1
! 插值部分
! 当z=nz，y=ny时x方向插值
      ax11=a(nx,ny,nz)+(a(nx1,ny,nz)-a(nx,ny,nz))*ddx/dlx
! 当z=nz，y=ny1时x方向插值
      ax12=a(nx,ny1,nz)+(a(nx1,ny1,nz)-a(nx,ny1,nz))*ddx/dlx
! 当z=nz时水平插值(在x方向插值基础上，对y方向插值)
      ay1=ax11+(ax12-ax11)*ddy/dly
! 当z=nz1，y=ny时x方向插值
      ax21=a(nx,ny,nz1)+(a(nx1,ny,nz1)-a(nx,ny,nz1))*ddx/dlx
! 当z=nz1，y=ny1时x方向插值
      ax22=a(nx,ny1,nz1)+(a(nx1,ny1,nz1)-a(nx,ny1,nz1))*ddx/dlx
! 当z=nz1时水平插值(在x方向插值基础上，对y方向插值)
      ay2=ax21+(ax22-ax21)*ddy/dly
! 将水平方向插值结果沿z方向插值得到最终的空间插值
      avg=ay1+(ay2-ay1)*ddz/dlz     
      return
      end
