      subroutine phases (pot,vv,s,u,a,b,aa,qq,eq)
c这个程序中vv,s,u,a,b...由调用程序分配，都是一维双精度实数组。其外部分配的大小远
c大于我们这个子程序需要用到的大小，所以不用传进来数组长度。又因为传递的只是指针，所以
c下面声明real*8 vv(1)在编译的时候不会出错
c****  in this version, the on-shell K-matrix is written out.
c
c
c        phases computes the r-matrix (or 'k-matrix') and  
c        the phase-shifts of nucleon-nucleon scattering
c        for a given nn-potential pot.
c        the calculations are performed in momentum space.
c        this package contains all subroutines needed, except for
c        subroutine gset that is contained in the package bonn.
c        furthermore, a potential subroutine is required, e.g. bonn.
c        an 'external' statement in the calling program has to specify
c        the name of the potential subroutine to be applied that will
c        replace 'pot' in this program.
c
c
c        this version of the code is to be published in "computational 
c        nuclear physics", vol. II, koonin, langanke, maruhn, eds.
c        (springer, heidelberg).
c
c
c        author: r. machleidt
c                department of physics
c                university of idaho
c                moscow, idaho 83843
c                u. s. a.
c                e-mail: machleid@phys.uidaho.edu
c
c                formerly:
c                institut fuer theoretische kernphysik bonn
c                nussallee 14-16
c                d-5300 bonn, w. germany
cgset需要用bonn package中的（是直接连接上还是重写是个问题），都行，可以先把bonn也放下去
cpot 是potential的简称，下面用external声明pot是个函数，所以如果不用bonn势call n3lo就行
      use phasecal
      use potential_global
      implicit real*8 (a-h,o-z)
      external pot
      real*8 elab(41)
      common /alpha/ melab
      common /einject/ elab
      common /crdwrt/ kread,kwrite,kpunch,kda(9)
c
c        arguments of the potential subroutine pot being called in this
c        program
c
      common /cstate/ jhang,heform,sing,trip,coup,endep,label
c        xmev and ymev are the final and initial relative momenta,
c        respectively, in units of mev.
c        v is the potential in units of mev**(-2).
c        if heform=.true., v must contain the 6 helicity matrix
c        elements associated with one j in the following order:
c        0v, 1v, 12v, 34v, 55v, 66v (helicity formalism).
c        if heform=.false., v must contain the six lsj-state matrix
c        elements associated with one j in the following order:
c        0v, 1v, v++, v--, v+-, v-+ (lsj formalism).
c        j is the total angular momentum. there is essentially
c        no upper limit for j.
c
c        specifications for these arguments
      logical heform,sing,trip,coup,endep
c
c        further specifications
c
c        all matrices are stored columnwise in vector arrays.
c        矩阵都是一列一列存的
c        see main program calling phases, cph, for the dimensions
c        of the following arrays
      real*8 vv(1),s(1),u(1),a(1),b(1),aa(1),qq(1),eq(1)
      real*8 q(97)
      real*8 delta(5)
      data delta/5*0.d0/
      real*8 rb(2)
      real*8 r(6)
      data uf/197.3286d0/
      data pih/1.570796326794897d0/
      real*4 ops
      data ops/1.e-15/
      data nalt/-1/
c        the following dimension for at most 40 elabs
      data memax/41/
      real*8 q0(40),q0q(40),eq0(40)
c!!!!!!放入mn，mp两个质量
      real*8 mn,mp 
      data mn/939.56563/
      data mp/938.27231/
c!!!!!!放入初值结束     
      integer nj(20)
c
      character*4 name(3),nname(15)
      character*1 state(2),state3
      character*1 multi(4)
      data multi/'1',3*'3'/
      integer ldel(4)
      data ldel/0,0,-1,1/
      character*1 spd(50)
      data spd/'s','p','d','f','g','h','i','k','l','m','n',
     1'o','q','r','t','u','v','w','x','y','z',29*' '/
      character*1 chars,chart
      data chars/'s'/,chart/'t'/
      character*1 lab1,lab2
      character*4 blanks
      data blanks/'    '/
      logical indbrn
      logical indrma
      logical indqua
      logical indwrt
      logical indpts
      logical inderg
      data indbrn/.false./
      data indrma/.false./
      data indqua/.false./
      data indwrt/.false./
      data indpts/.false./
      data inderg/.false./
      logical indj
c
c
c
c
9998  format (41f10.4)
9999  format (f10.6)
10000 format (2a4,a2,20i3)
10001 format (1h ,2a4,a2,20i3)
10002 format (2a4,a2,6f10.4)
10003 format (1h ,2a4,a2,6f10.4)
10004 format (//' input-parameters for phases'/1h ,27(1h-))
10005 format (//' transformed gauss points and weights for c =',f8.2,
     1'  and n =',i3//' points')
10006 format (7x,4f15.4)
10007 format (/' weights')
10008 format (2a4,a2,15a4)
10009 format (1h ,2a4,a2,15a4)
10013 format (///' error in phases. matrix inversion error index =',
     1 i5///)
10014 format (1h ,2a4,a2,6f10.4)
10016 format (/' low energy parameters    a',a1,' =',
     1 f10.4,'    r',a1,' =',f10.4)
10050 format (//' elab(mev)',5x,4(2a1,i2,10x),' e',i2/)
10051 format (1h ,f8.2,5e14.6)
10052 format (//' elab(mev)',5x,'1s0',25x,'3p0'/)
10053 format (1h1//' p h a s e - s h i f t s (radians)'/1h ,33(1h-))
10054 format (1h1//' p h a s e - s h i f t s (degrees)'/1h ,33(1h-))
10055 format (4f14.6)
10110 format (1h ,a4,2x,2a1,'-',a1,i1,3d20.12)
10111 format (i3)
10112 format (' elab (mev) ',f10.4)
10113 format (4d20.12)
      lambda=0.0d0
c
c上述是为了下面输出文件格式写得一些标注行
cwrite可以对照输出文件看，kwrite是6，代表输出文件
c
c        read and write input parameters for this program
c        ------------------------------------------------
c
c
c     write (kwrite,10004)
c
c        read and write comment

c        read and write the range for the total angular momentum j
c        for which phase shifts are to be calculated. there are
c        essentially no limitations for j.
      jb=0
      je=phasename(phasenum,1)
c      write (kwrite,10001) name,jb,je
c jb，je是输出相移分波中j量子数的上下限制
c        the born approximation will be used for j.ge.jborn.
c  在Fortran77中使用.ge.缩写左右加点得到逻辑判断符
      jborn=9
c      write (kwrite,10001) name,jborn
      jb1=jb+1
      je1=je+1
      jee1=je1
      if (jee1.gt.20) jee1=20
c  相移一共有三中计算方式，当j>=jborn时使用born近似计算散射矩阵元
c        read and write the number of gauss points to be used for the
c        matrix inversion; this is j-dependent.
      nj=0
      nj(1)=32
c      write (kwrite,10001) name,(nj(j1),j1=1,jee1)
c   这个不太清楚应该不用改
c        ihef=0: lsj formalism is used (heform=.false.),
c        ihef.ne.0: helicity formalism is used (heform=.true.).
      ihef=0
c      write (kwrite,10001) name,ihef
c   ihef就是核力文件中的heform
c        set ising=itrip=icoup=1, if you want that for each j
c        all possible states are considered.
      ising=1
      itrip=1
      icoup=1      
c      write (kwrite,10001) name,ising,itrip,icoup
c   这个核力文件写的很清楚
c        iprop=1: non-relativistic propagator in scattering equation;
c        iprop=2: relativistic propagator.
      iprop=1
c      write (kwrite,10001) name,iprop
c
c        iphrel=1: non-relativistic phase-relation;
c        iphrel=2: relativistic phase-relation.
      iphrel=1
c      write (kwrite,10001) name,iphrel
c   这两部分是计算相移是否用相对论方程
c        c is the factor involved in the transformation of the
c        gauss points.
      c=1000.0d0
c      write (kwrite,10003) name,c
c   这个也不清楚，应该是计算问题
c        wn is the nucleon mass.
      wn=938.9183
c      write (kwrite,10003) name,wn
c   wn是核子质量
      
c        read and write the elab's for which phase shifts are to be
c        calculated. any number of elab's between 1 and 40 is permissible.
c        the last elab has to be zero and is understood as the end
c        of the list of elab's.
c   读入入射动能，可以读40个数据，d0是dimension0乘10的零次方
      read (kread,*) melab
      allocate(phaseshifts(melab,phasenum))
      phaseshifts=0.0d0
      do 1 i=1,melab
    1 read (kread,*) elab(i)

c
c        irma=0:    the r-matrix is not written.
c        irma.ne.0: the r-matrix is written in terms of lsj states
c                   to unit kpunch,
c                   iqua=0: only the half off-shell r-matrix is written,
c                   iqua.ne.0: the (quadratic) fully off-shell r-matrix 
c                              is calculated and written.
c   irma是零不输出r矩阵，irma不是零把r矩阵输出到irma文件中；iqua是r矩阵内部的的信息
      irma=0
      iqua=0
c      write (kwrite,10001) name,irma,iqua
c
c        ipoint=0: the transformed gauss-points and -weights are
c                  not printed;
c        ipoint.ne.0: ... are printed.
      ipoint=0
c      write (kwrite,10001) name,ipoint
c
c        ideg=0: phase-shifts are printed in radians;
c        ideg=1: phase-shifts are printed in degrees.
      ideg=1
c      write (kwrite,10001) name,ideg
c
c
c        prepare constants
c        -----------------
c
c
      sing=.false.
      trip=.false.
      coup=.false.
      do i=1,24
      read (kread,9999) conta(i)
      end do
      read(kread,*) lambda
      if (ising.ne.0) sing=.true.
      if (itrip.ne.0) trip=.true.
      if (icoup.ne.0) coup=.true.
      heform=.false.
      if (ihef.ne.0) heform=.true.
      if (irma.ne.0) indrma=.true.
      if (irma.ne.0.and.iqua.ne.0) indqua=.true.
      if (ipoint.ne.0) indpts=.true.
c
      if(jb.le.1.and.elab(1).le.2.d0.and.elab(2).le.2.d0.and.melab.ge.2)
     1inderg=.true.
      wnh=wn*0.5d0
      wnq=wn*wn
      wp=wn*pih
      rd=90.d0/pih
      iideg=ideg+1
 
c**** label=blanks
c
c
c
c        prepare energies and on-shell momenta
c !!!!!!!!!!对于np计算我们把相移中动量修正的程序放进去，这段是原来的程序
c     do 10 k=1,melab
c     q0q(k)=wnh*elab(k)
c     q0(k)=dsqrt(q0q(k))
c  10 eq0(k)=dsqrt(q0q(k)+wnq)
c!!!!!!!!!!这里是原程序结尾
      do 10 k=1,melab
      q0q(k)=(mp*mp*elab(k)*(elab(k)+2.0d0*mn))
     &/((mn+mp)*(mn+mp)+2.0d0*elab(k)*mp)
      q0(k)=dsqrt(q0q(k))
c这一段是我们主要修改的
   10 eq0(k)=dsqrt(q0q(k)+wnq)

c
c     
c
c
c        loop of total angular momentum j
c        --------------------------------
c        --------------------------------
c
c
c
c
    
      do 2000 j1=jb1,je1
c
c
      indj=.false.
      j2=j1+1
      j=j1-1
      aj=dfloat(j)
      aj1=dfloat(j+1)
      a2j1=dfloat(2*j+1)
      d2j1=1.d0/a2j1
      aaj=dsqrt(aj*aj1)
c
c
      if (j.ge.jborn) indbrn=.true.
c
c
      if (j1.le.20) go to 105
      if (nalt.ge.1) go to 300
      n=16
      go to 115
c
c        number of gausspoints for this j
  105 if (nj(j1).eq.0) nj(j1)=nj(j)
      n=nj(j1)
      if (n.eq.nalt) go to 300
c
c
c        get gauss points and weights
c
  115 call gset (0.d0,1.d0,n,u,s)
c
      nalt=n
      n1=n+1
      n2=2*n1
      n3=3*n1
      n4=4*n1
      nx=n1
      if (indqua) nx=nx*nx
      nx2=2*nx
      nx2mn=nx2-n1
      nx2pn=nx2+n1
      nx4=4*nx
      nx4mn=nx4-n1
c
c        transform gauss points and weights
c
      do 201 i=1,n
      xx=pih*u(i)
c
c        transformed gauss point
      q(i)=dtan(xx)*c
      qq(i)=q(i)*q(i)
      eq(i)=dsqrt(qq(i)+wnq)
c
c        transformed gauss weight
      dc=1.d0/dcos(xx)
  201 s(i)=pih*c*dc*dc*s(i)
c
      if (.not.indqua) go to 205
c
c      write (kpunch,10111) n
c      write (kpunch,10113) (q(i),i=1,n)
c
  205 if (.not.indpts) go to 300
c
c        write gauss points and weights
c
c     write (kwrite,10005) c,n
c     write (kwrite,10006) (q(i),i=1,n)
c
c     write (kwrite,10007)
c     write (kwrite,10006) (s(i),i=1,n)
c
c
c
c
c
c
c        loop of elabs
c        -------------
c        -------------
c
c
c
c
  300 do 1000 k=1,melab
cfortran77 用行号标记end do
      q(n1)=q0(k)
c
c
      if (indrma)  continue!write (kpunch,10112) elab(k)
     
c
c
      if (indbrn.and..not.indqua) go to 500
c
c
c        check if right potential matrix does already exist
      if (indj) go to 500
      indj=.true.
c
c
c        compute potential matrix
c        ------------------------
c
c
      iii=0
      do 401 ix=1,n
c
c
      xmev=q(ix)
c
c
      do 401 iy=ix,n
c
c
      ymev=q(iy)
c
c
      call pot
c
c
      iaa=iii*6
      iii=iii+1
      do 401 iv=1,6
  401 aa(iv+iaa)=v(iv)
c
c
c        compute potential vector
c        ------------------------
c
c
  500 if (indbrn.and..not.indrma) go to 510
c
      ymev=q0(k)
c
c
      do 501 ix=1,n
c
c
      xmev=q(ix)
c
c
      call pot
c
c
      do 501 iv=1,6
      ivv=ix+(iv-1)*n
  501 vv(ivv)=v(iv)
c
c
c        compute potential element
c        -------------------------
c
c
  510 xmev=q0(k)
      ymev=q0(k)
c
c
      call pot
c
c
c
c        compute factor for the phase relation
c
      go to (601,602),iphrel
  601 wpq0=-wp*q0(k)
      go to 605
  602 wpq0=-pih*eq0(k)*q0(k)
  605 continue
c
c
      if (indbrn.and..not.indqua) go to 700
c
c
c        compute propagator
c        ------------------
c
c
      uq0=0.d0
      do 620 i=1,n
      sdq=s(i)/(qq(i)-q0q(k))
c
c        calculate propagator of lippmann-schwinger equation
c
      go to (621,622),iprop
  621 u(i)=sdq*qq(i)*wn
      go to 620
  622 u(i)=sdq*qq(i)*(eq(i)+eq0(k))*0.5d0
c
  620 uq0=uq0+sdq
c
c
      go to (631,632),iprop
  631 uq0=-uq0*q0q(k)*wn
      go to 700
  632 uq0=-uq0*q0q(k)*eq0(k)
c
c
c
c
c
c        build up matrix to be inverted
c        ------------------------------
c
c
  700 ni=0
      nii=0
      nv=n1
      mv=1
      if (indqua) mv=mv*n1
      ib=0
      eins=1.d0
c
c
      if (.not.sing) go to 720
      iv=1
      go to 770
c
c
  720 if (.not.trip.or.j.eq.0) go to 730
      iv=2
      go to 770
c
c
  730 if (.not.coup) go to 900
      iv=3
      if (j.eq.0) go to 770
      nv=n2
      mv=2
      if (indqua) mv=mv*n1
      go to 770
c
  740 if (j.eq.0) go to 800
      iv=4
      ib=n3
      ni=n1
      nii=n1
      go to 770
c
  750 iv=5
      ivi=6
      ib=n2
      ni=0
      nii=n1
      eins=0.d0
      go to 770
c
  760 iv=6
      ivi=5
      ib=n1
      ni=n1
      nii=0
c
c
c
c
  770 iii=0
      if (iv.le.4) ivi=iv
      igg=(iv-1)*n
      i1=(nii+n)*nv
      i2=(nii-1)*nv
c
c
      if (indbrn.and..not.indrma) go to 785
c
c
      do 780 i=1,n
      i3=i*nv
      i4=ni+i
c
c
      do 781 ii=i,n
      iaa=iii*6
      iii=iii+1
      i5=i2+i3+ni+ii
      i6=ivi+iaa
      i7=i2+i4+ii*nv
      i8=iv+iaa
      if (i.eq.ii) go to 782
c
c        matrix a
      a(i7)=aa(i8)*u(ii)
      a(i5)=aa(i6)*u(i)
      if (.not.indqua) go to 781
c
c        matrix b
      b(i7)=aa(i8)
      b(i5)=aa(i6)
      go to 781
c        diagonal element
  782 a(i7)=aa(i8)*u(i)+eins
      if (.not.indqua) go to 781
      b(i7)=aa(i8)
  781 continue
c
c        last column
      i9=i1+i4
      i10=i+igg
      a(i9)=vv(i10)*uq0
c        last row
      i11=i2+i3+ni+n1
      ivv=i+(ivi-1)*n
      a(i11)=vv(ivv)*u(i)
      if (.not.indqua) go to 783
      b(i9)=vv(i10)
      b(i11)=vv(ivv)
      go to 780
c
c        vector b
  783 b(ib+i)=vv(i+igg)
c
  780 continue
c
c
c        last element
      i12=i1+ni+n1
      a(i12)=v(iv)*uq0+eins
      if (.not.indqua) go to 785
      b(i12)=v(iv)
      go to 790
  785 b(ib+n1)=v(iv)
c
c
c
c
  790 go to (800,800,740,750,760,800),iv
c
c
c
c
c        invert matrix
c        -------------
c
c
c
c
  800 if (indbrn) go to 801
      call dgelg (b,a,nv,mv,ops,ier)
c
c
      if (ier.ne.0) continue!write(kwrite,10013) ier
      
c
c
c
c
c        compute phase shifts
c        --------------------
c
c
c
c
  801 if (iv.gt.2.and.j.ne.0) go to 820
c
c        uncoupled cases
c
      delta(iv)=datan(b(nx)*wpq0)
c
c        prepare for effective range
      if (inderg.and.j.eq.0.and.iv.eq.1.and.k.le.2) 
     1rb(k)=q0(k)/(b(nx)*wpq0)
c
c
      if (.not.indrma) go to 810
c
c
c
c        write r-matrix
c
c
      state(1)=multi(iv)
      ispd=j1+ldel(iv)
      if (j.eq.0.and.iv.eq.3) ispd=2
      state(2)=spd(ispd)
      state3=state(2)
      if (indqua) continue!write (kpunch,10110) label,state,state3,j 
      do 805 i=n1,n1
      if (indqua) go to 804
c
c        write half off-shell r-matrix
c      write (kpunch,10110) label,state,state3,j,q(i),q0(k),b(i)
      go to 805
c
c        write fully off-shell r-matrix (lower triangle)
  804 i1=(i-1)*n1
c      write (kpunch,10113) (b(i1+ii),ii=i,n1)
  805 continue
c
c
c
c
  810 go to (720,730,900),iv
c
c
c        coupled cases
c
  820 if (heform) go to 822
c
c        calculate phase shifts from lsj-state r-matrix elements
c
      r0=b(nx2) 
      r1=b(nx2mn)-b(nx4)
      r2=b(nx2mn)+b(nx4)
      rr=-2.d0*r0/r1
c        epsilon
      delta(5)=datan(rr)/2.d0
      rr=r1*dsqrt(1.d0+rr*rr)
c        prepare for effective range
      if (inderg.and.j.eq.1.and.k.le.2) 
     1rb(k)=2.d0*q0(k)/(wpq0*(r2-rr))
c        delta minus
      delta(3)=datan((r2-rr)*wpq0*0.5d0)
c        delta plus
      delta(4)=datan((r2+rr)*wpq0*0.5d0)
      go to 824
c
c
c        calculate phase shifts from helicity-state r-matrix elements
c
  822 r0=b(nx2)
      r1=b(nx2mn)-b(nx4)
      r2=b(nx2mn)+b(nx4)
      rr=-2.d0*(aaj*r1+r0)/(r1-4.d0*aaj*r0)
c        epsilon
      delta(5)=datan(rr)/2.d0
      rr=(r1-4.d0*aaj*r0)*dsqrt(1.d0+rr*rr)*d2j1
c        prepare for effective range
      if (inderg.and.j.eq.1.and.k.le.2) 
     1rb(k)=2.d0*q0(k)/(wpq0*(r2-rr))
c        delta minus
      delta(3)=datan((r2-rr)*wpq0*0.5d0)
c        delta plus
      delta(4)=datan((r2+rr)*wpq0*0.5d0)
c
c        so far the delta(..) have been the blatt-biedenharn phase-
c        shifts, transform the delta(..) now into bar-phase-shifts
c        according to stapp et al., phys. rev. 105 (1957) 302.
  824 if (delta(5).eq.0.d0) go to 829
      pp=delta(3)+delta(4)
      pm=delta(3)-delta(4)
      d52=2.d0*delta(5)
      pm2=dsin(d52)*dsin(pm)
      delta(5)=0.5d0*dasin(pm2)
      pm1=dtan(2.d0*delta(5))/dtan(d52)
      pm1=dasin(pm1)
      delta(3)=0.5d0*(pp+pm1)
      delta(4)=0.5d0*(pp-pm1)
  829 continue
c
c
      if (.not.indrma) go to 900
c
c
c
      na=1
      if (indqua) na=n1
      do 835 i=1,na
      i1=(i-1)*n2
      do 835 ii=i,n1
      i3=i1+ii
      i4=nx2pn+i3
      i5=nx2+i3
      i6=n1+i3
      r(3)=b(i3)
      r(4)=b(i4)
      r(5)=b(i5)
      r(6)=b(i6)
      if (.not.heform) go to 837
c
c        in case of heform=.true., transform into lsj-form
      r34=(r(3)-r(4))*aaj
      r56=(r(5)+r(6))*aaj
      b(i6)=(aj1*r(3)+aj*r(4)-r56)*d2j1
      b(i3)=(aj*r(3)+aj1*r(4)+r56)*d2j1
      b(i5)=(r34+aj1*r(5)-aj*r(6))*d2j1
      b(i4)=(r34-aj*r(5)+aj1*r(6))*d2j1
      go to 835
c
c        in case of heform=.false., reorganize
  837 b(i6)=r(3)
      b(i3)=r(4)
      b(i5)=r(5)
      b(i4)=r(6)
  835 continue
c
c        write r-matrix
      ivx=0
      do 840 iv1=3,4
      do 840 iv2=3,4
      ivx=ivx+1
      state(1)=multi(iv1)
      state(2)=spd(j1+ldel(iv1))
      state3=spd(j1+ldel(iv2))
      go to (831,832,833,834),ivx
  831 ny=0
      go to 836
  832 ny=nx2pn
      go to 836
  833 ny=nx2
      go to 836
  834 ny=n1
c
c
  836 if (indqua) continue!write (kpunch,10110) label,state,state3,j
      do 839 i=n1,n1
      if (indqua) go to 838
c
c        write half off-shell r-matrix
      ii1=ny+i
c      write (kpunch,10110) label,state,state3,j,q(i),q0(k),b(ii1)
      go to 839
c
c        write fully off-shell r-matrix (lower triangle)
  838 i1=(i-1)*n2+ny
c     write (kpunch,10113) (b(i1+ii),ii=i,n1)
  839 continue
  840 continue
c
c
c
c
  900 continue
      
c
c
c        write phase-shifts
c        ------------------
c
      if (indwrt) go to 921
      indwrt=.true.
      go to (931,932),iideg
  931 continue!write (kwrite,10053)
      go to 933
  932 continue!write (kwrite,10054)
  933 continue
c
  921 if (k.ne.1) go to 923
      if (j.ne.0) go to 922
c      write (kwrite,10052)
      go to 923
  922 continue !write (kwrite,10050) multi(1),spd(j1),j,
c    1                     multi(2),spd(j1),j,
c     2                     multi(2),spd(j),j,
c     3                     multi(2),spd(j2),j,
c     4                     j
  923 if (iideg.eq.1) go to 926
      do 925 iv=1,5
  925 delta(iv)=delta(iv)*rd
      do i=1,phasenum
            if (phasename(i,1).eq.j)then
                  phaseshifts(k,i)=delta(phasename(i,2))
            end if
      end do 
      go to 926
c      write (kwrite,10055) elab(k),delta
  926 continue!write (kwrite,10051) elab(k),delta
       

c
c
c
c
 1000 continue
c        this has been the end of the elab loop
c
c
c
c
c        calculate and write low energy parameters
c        -----------------------------------------
c
c
      if (.not.inderg) go to 2000
      if (j.gt.1) go to 2000
      rb2=4.d0/wn*(rb(1)-rb(2))/(elab(1)-elab(2))*uf
      rb1=wn/4.*elab(1)*rb2/uf-rb(1)
      rb1=1./rb1*uf
      if (j.ne.0) go to 1090
      lab1=chars
      lab2=chars
      go to 1091
 1090 lab1=chart
      lab2=chart
 1091 continue!write (kwrite,10016) lab1,rb1,lab2,rb2
c
c
c
c
 2000 continue
c        this has been the end of the j loop
c
c
      end
c*************************************************************
c
c name:      dgelg
c
c from:      programmbibliothek rhrz bonn, germany;   02/02/81
c            (free soft-ware)
c language:  fortran iv (fortran-77 compatible)
c
c purpose:
c
c to solve a general system of simultaneous linear equations.
c
c usage:   call dgelg(r,a,m,n,eps,ier)
c
c parameters:
c
c r:       double precision m by n right hand side matrix
c          (destroyed). on return r contains the solutions
c          of the equations.
c
c a:       double precision m by m coefficient matrix
c          (destroyed).
c
c m:       the number of equations in the system.
c
c n:       the number of right hand side vectors.
c
c eps:     single precision input constant which is used as
c          relative tolerance for test on loss of
c          significance.
c
c ier:     resulting error parameter coded as follows
c           ier=0  - no error,
c           ier=-1 - no result because of m less than 1 or
c                   pivot element at any elimination step
c                   equal to 0,
c           ier=k  - warning due to possible loss of signifi-
c                   cance indicated at elimination step k+1,
c                   where pivot element was less than or
c                   equal to the internal tolerance eps times
c                   absolutely greatest element of matrix a.
c
c remarks: (1) input matrices r and a are assumed to be stored
c              columnwise in m*n resp. m*m successive storage
c              locations. on return solution matrix r is stored
c              columnwise too.
c          (2) the procedure gives results if the number of equations m
c              is greater than 0 and pivot elements at all elimination
c              steps are different from 0. however warning ier=k - if
c              given indicates possible loss of significance. in case
c              of a well scaled matrix a and appropriate tolerance eps,
c              ier=k may be interpreted that matrix a has the rank k.
c              no warning is given in case m=1.
c
c method:
c
c solution is done by means of gauss-elimination with
c complete pivoting.
c
c
c author:         ibm, ssp iii
c
c**********************************************************************
      subroutine gset(ax,bx,n,z,w)
c
c
c        this code has been obtained from the CERN computer library
c        in the year of the lord 1972.
c
c
      implicit real*8 (a-h,o-z)
c
c     n-point gauss zeros and weights for the interval (ax,bx) are
c           stored in  arrays z and w respectively.
c
      dimension     a(273),x(273),ktab(96)
      dimension z(2),w(2)
c
c-----table of initial subscripts for n=2(1)16(4)96
      data ktab(2)/1/
      data ktab(3)/2/
      data ktab(4)/4/
      data ktab(5)/6/
      data ktab(6)/9/
      data ktab(7)/12/
      data ktab(8)/16/
      data ktab(9)/20/
      data ktab(10)/25/
      data ktab(11)/30/
      data ktab(12)/36/
      data ktab(13)/42/
      data ktab(14)/49/
      data ktab(15)/56/
      data ktab(16)/64/
      data ktab(20)/72/
      data ktab(24)/82/
      data ktab(28)/82/
      data ktab(32)/94/
      data ktab(36)/94/
      data ktab(40)/110/
      data ktab(44)/110/
      data ktab(48)/130/
      data ktab(52)/130/
      data ktab(56)/130/
      data ktab(60)/130/
      data ktab(64)/154/
      data ktab(68)/154/
      data ktab(72)/154/
      data ktab(76)/154/
      data ktab(80)/186/
      data ktab(84)/186/
      data ktab(88)/186/
      data ktab(92)/186/
      data ktab(96)/226/
c
c-----table of abscissae (x) and weights (a) for interval (-1,+1).
c
c**** n=2
      data x(1)/0.577350269189626  d0/, a(1)/1.000000000000000  d0/
c**** n=3
      data x(2)/0.774596669241483  d0/, a(2)/0.555555555555556  d0/
      data x(3)/0.000000000000000  d0/, a(3)/0.888888888888889  d0/
c**** n=4
      data x(4)/0.861136311594053  d0/, a(4)/0.347854845137454  d0/
      data x(5)/0.339981043584856  d0/, a(5)/0.652145154862546  d0/
c**** n=5
      data x(6)/0.906179845938664  d0/, a(6)/0.236926885056189  d0/
      data x(7)/0.538469310105683  d0/, a(7)/0.478628670499366  d0/
      data x(8)/0.000000000000000  d0/, a(8)/0.568888888888889  d0/
c**** n=6
      data x(9)/0.932469514203152  d0/, a(9)/0.171324492379170  d0/
      data x(10)/0.661209386466265 d0/, a(10)/0.360761573048139 d0/
      data x(11)/0.238619186083197 d0/, a(11)/0.467913934572691 d0/
c**** n=7
      data x(12)/0.949107912342759 d0/, a(12)/0.129484966168870 d0/
      data x(13)/0.741531185599394 d0/, a(13)/0.279705391489277 d0/
      data x(14)/0.405845151377397 d0/, a(14)/0.381830050505119 d0/
      data x(15)/0.000000000000000 d0/, a(15)/0.417959183673469 d0/
c**** n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c**** n=9
      data x(20)/0.968160239507626 d0/, a(20)/0.081274388361574 d0/
      data x(21)/0.836031107326636 d0/, a(21)/0.180648160694857 d0/
      data x(22)/0.613371432700590 d0/, a(22)/0.260610696402935 d0/
      data x(23)/0.324253423403809 d0/, a(23)/0.312347077040003 d0/
      data x(24)/0.000000000000000 d0/, a(24)/0.330239355001260 d0/
c**** n=10
      data x(25)/0.973906528517172 d0/, a(25)/0.066671344308688 d0/
      data x(26)/0.865063366688985 d0/, a(26)/0.149451349150581 d0/
      data x(27)/0.679409568299024 d0/, a(27)/0.219086362515982 d0/
      data x(28)/0.433395394129247 d0/, a(28)/0.269266719309996 d0/
      data x(29)/0.148874338981631 d0/, a(29)/0.295524224714753 d0/
c**** n=11
      data x(30)/0.978228658146057 d0/, a(30)/0.055668567116174 d0/
      data x(31)/0.887062599768095 d0/, a(31)/0.125580369464905 d0/
      data x(32)/0.730152005574049 d0/, a(32)/0.186290210927734 d0/
      data x(33)/0.519096129206812 d0/, a(33)/0.233193764591990 d0/
      data x(34)/0.269543155952345 d0/, a(34)/0.262804544510247 d0/
      data x(35)/0.000000000000000 d0/, a(35)/0.272925086777901 d0/
c**** n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c**** n=13
      data x(42)/0.984183054718588 d0/, a(42)/0.040484004765316 d0/
      data x(43)/0.917598399222978 d0/, a(43)/0.092121499837728 d0/
      data x(44)/0.801578090733310 d0/, a(44)/0.138873510219787 d0/
      data x(45)/0.642349339440340 d0/, a(45)/0.178145980761946 d0/
      data x(46)/0.448492751036447 d0/, a(46)/0.207816047536889 d0/
      data x(47)/0.230458315955135 d0/, a(47)/0.226283180262897 d0/
      data x(48)/0.000000000000000 d0/, a(48)/0.232551553230874 d0/
c**** n=14
      data x(49)/0.986283808696812 d0/, a(49)/0.035119460331752 d0/
      data x(50)/0.928434883663574 d0/, a(50)/0.080158087159760 d0/
      data x(51)/0.827201315069765 d0/, a(51)/0.121518570687903 d0/
      data x(52)/0.687292904811685 d0/, a(52)/0.157203167158194 d0/
      data x(53)/0.515248636358154 d0/, a(53)/0.185538397477938 d0/
      data x(54)/0.319112368927890 d0/, a(54)/0.205198463721296 d0/
      data x(55)/0.108054948707344 d0/, a(55)/0.215263853463158 d0/
c**** n=15
      data x(56)/0.987992518020485 d0/, a(56)/0.030753241996117 d0/
      data x(57)/0.937273392400706 d0/, a(57)/0.070366047488108 d0/
      data x(58)/0.848206583410427 d0/, a(58)/0.107159220467172 d0/
      data x(59)/0.724417731360170 d0/, a(59)/0.139570677926154 d0/
      data x(60)/0.570972172608539 d0/, a(60)/0.166269205816994 d0/
      data x(61)/0.394151347077563 d0/, a(61)/0.186161000015562 d0/
      data x(62)/0.201194093997435 d0/, a(62)/0.198431485327111 d0/
      data x(63)/0.000000000000000 d0/, a(63)/0.202578241925561 d0/
c**** n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/
c**** n=20
      data x(72)/0.993128599185094 d0/, a(72)/0.017614007139152 d0/
      data x(73)/0.963971927277913 d0/, a(73)/0.040601429800386 d0/
      data x(74)/0.912234428251325 d0/, a(74)/0.062672048334109 d0/
      data x(75)/0.839116971822218 d0/, a(75)/0.083276741576704 d0/
      data x(76)/0.746331906460150 d0/, a(76)/0.101930119817240 d0/
      data x(77)/0.636053680726515 d0/, a(77)/0.118194531961518 d0/
      data x(78)/0.510867001950827 d0/, a(78)/0.131688638449176 d0/
      data x(79)/0.373706088715419 d0/, a(79)/0.142096109318382 d0/
      data x(80)/0.227785851141645 d0/, a(80)/0.149172986472603 d0/
      data x(81)/0.076526521133497 d0/, a(81)/0.152753387130725 d0/
c**** n=24
      data x(82)/0.995187219997021 d0/, a(82)/0.012341229799987 d0/
      data x(83)/0.974728555971309 d0/, a(83)/0.028531388628933 d0/
      data x(84)/0.938274552002732 d0/, a(84)/0.044277438817419 d0/
      data x(85)/0.886415527004401 d0/, a(85)/0.059298584915436 d0/
      data x(86)/0.820001985973902 d0/, a(86)/0.073346481411080 d0/
      data x(87)/0.740124191578554 d0/, a(87)/0.086190161531953 d0/
      data x(88)/0.648093651936975 d0/, a(88)/0.097618652104113 d0/
      data x(89)/0.545421471388839 d0/, a(89)/0.107444270115965 d0/
      data x(90)/0.433793507626045 d0/, a(90)/0.115505668053725 d0/
      data x(91)/0.315042679696163 d0/, a(91)/0.121670472927803 d0/
      data x(92)/0.191118867473616 d0/, a(92)/0.125837456346828 d0/
      data x(93)/0.064056892862605 d0/, a(93)/0.127938195346752 d0/
c**** n=32
      data x(94)/0.997263861849481 d0/, a(94)/0.007018610009470 d0/
      data x(95)/0.985611511545268 d0/, a(95)/0.016274394730905 d0/
      data x(96)/0.964762255587506 d0/, a(96)/0.025392065309262 d0/
      data x(97)/0.934906075937739 d0/, a(97)/0.034273862913021 d0/
      data x(98)/0.896321155766052 d0/, a(98)/0.042835898022226 d0/
      data x(99)/0.849367613732569 d0/, a(99)/0.050998059262376 d0/
      data x(100)/0.794483795967942d0/, a(100)/0.058684093478535d0/
      data x(101)/0.732182118740289d0/, a(101)/0.065822222776361d0/
      data x(102)/0.663044266930215d0/, a(102)/0.072345794108848d0/
      data x(103)/0.587715757240762d0/, a(103)/0.078193895787070d0/
      data x(104)/0.506899908932229d0/, a(104)/0.083311924226946d0/
      data x(105)/0.421351276130635d0/, a(105)/0.087652093004403d0/
      data x(106)/0.331868602282127d0/, a(106)/0.091173878695763d0/
      data x(107)/0.239287362252137d0/, a(107)/0.093844399080804d0/
      data x(108)/0.144471961582796d0/, a(108)/0.095638720079274d0/
      data x(109)/0.048307665687738d0/, a(109)/0.096540088514727d0/
c**** n=40
      data x(110)/0.998237709710559d0/, a(110)/0.004521277098533d0/
      data x(111)/0.990726238699457d0/, a(111)/0.010498284531152d0/
      data x(112)/0.977259949983774d0/, a(112)/0.016421058381907d0/
      data x(113)/0.957916819213791d0/, a(113)/0.022245849194166d0/
      data x(114)/0.932812808278676d0/, a(114)/0.027937006980023d0/
      data x(115)/0.902098806968874d0/, a(115)/0.033460195282547d0/
      data x(116)/0.865959503212259d0/, a(116)/0.038782167974472d0/
      data x(117)/0.824612230833311d0/, a(117)/0.043870908185673d0/
      data x(118)/0.778305651426519d0/, a(118)/0.048695807635072d0/
      data x(119)/0.727318255189927d0/, a(119)/0.053227846983936d0/
      data x(120)/0.671956684614179d0/, a(120)/0.057439769099391d0/
      data x(121)/0.612553889667980d0/, a(121)/0.061306242492928d0/
      data x(122)/0.549467125095128d0/, a(122)/0.064804013456601d0/
      data x(123)/0.483075801686178d0/, a(123)/0.067912045815233d0/
      data x(124)/0.413779204371605d0/, a(124)/0.070611647391286d0/
      data x(125)/0.341994090825758d0/, a(125)/0.072886582395804d0/
      data x(126)/0.268152185007253d0/, a(126)/0.074723169057968d0/
      data x(127)/0.192697580701371d0/, a(127)/0.076110361900626d0/
      data x(128)/0.116084070675255d0/, a(128)/0.077039818164247d0/
      data x(129)/0.038772417506050d0/, a(129)/0.077505947978424d0/
c**** n=48
      data x(130)/0.998771007252426d0/, a(130)/0.003153346052305d0/
      data x(131)/0.993530172266350d0/, a(131)/0.007327553901276d0/
      data x(132)/0.984124583722826d0/, a(132)/0.011477234579234d0/
      data x(133)/0.970591592546247d0/, a(133)/0.015579315722943d0/
      data x(134)/0.952987703160430d0/, a(134)/0.019616160457355d0/
      data x(135)/0.931386690706554d0/, a(135)/0.023570760839324d0/
      data x(136)/0.905879136715569d0/, a(136)/0.027426509708356d0/
      data x(137)/0.876572020274247d0/, a(137)/0.031167227832798d0/
      data x(138)/0.843588261624393d0/, a(138)/0.034777222564770d0/
      data x(139)/0.807066204029442d0/, a(139)/0.038241351065830d0/
      data x(140)/0.767159032515740d0/, a(140)/0.041545082943464d0/
      data x(141)/0.724034130923814d0/, a(141)/0.044674560856694d0/
      data x(142)/0.677872379632663d0/, a(142)/0.047616658492490d0/
      data x(143)/0.628867396776513d0/, a(143)/0.050359035553854d0/
      data x(144)/0.577224726083972d0/, a(144)/0.052890189485193d0/
      data x(145)/0.523160974722233d0/, a(145)/0.055199503699984d0/
      data x(146)/0.466902904750958d0/, a(146)/0.057277292100403d0/
      data x(147)/0.408686481990716d0/, a(147)/0.059114839698395d0/
      data x(148)/0.348755886292160d0/, a(148)/0.060704439165893d0/
      data x(149)/0.287362487355455d0/, a(149)/0.062039423159892d0/
      data x(150)/0.224763790394689d0/, a(150)/0.063114192286254d0/
      data x(151)/0.161222356068891d0/, a(151)/0.063924238584648d0/
      data x(152)/0.097004699209462d0/, a(152)/0.064466164435950d0/
      data x(153)/0.032380170962869d0/, a(153)/0.064737696812683d0/
c**** n=64
      data x(154)/0.999305041735772d0/, a(154)/0.001783280721696d0/
      data x(155)/0.996340116771955d0/, a(155)/0.004147033260562d0/
      data x(156)/0.991013371476744d0/, a(156)/0.006504457968978d0/
      data x(157)/0.983336253884625d0/, a(157)/0.008846759826363d0/
      data x(158)/0.973326827789910d0/, a(158)/0.011168139460131d0/
      data x(159)/0.961008799652053d0/, a(159)/0.013463047896718d0/
      data x(160)/0.946411374858402d0/, a(160)/0.015726030476024d0/
      data x(161)/0.929569172131939d0/, a(161)/0.017951715775697d0/
      data x(162)/0.910522137078502d0/, a(162)/0.020134823153530d0/
      data x(163)/0.889315445995114d0/, a(163)/0.022270173808383d0/
      data x(164)/0.865999398154092d0/, a(164)/0.024352702568710d0/
      data x(165)/0.840629296252580d0/, a(165)/0.026377469715054d0/
      data x(166)/0.813265315122797d0/, a(166)/0.028339672614259d0/
      data x(167)/0.783972358943341d0/, a(167)/0.030234657072402d0/
      data x(168)/0.752819907260531d0/, a(168)/0.032057928354851d0/
      data x(169)/0.719881850171610d0/, a(169)/0.033805161837141d0/
      data x(170)/0.685236313054233d0/, a(170)/0.035472213256882d0/
      data x(171)/0.648965471254657d0/, a(171)/0.037055128540240d0/
      data x(172)/0.611155355172393d0/, a(172)/0.038550153178615d0/
      data x(173)/0.571895646202634d0/, a(173)/0.039953741132720d0/
      data x(174)/0.531279464019894d0/, a(174)/0.041262563242623d0/
      data x(175)/0.489403145707052d0/, a(175)/0.042473515123653d0/
      data x(176)/0.446366017253464d0/, a(176)/0.043583724529323d0/
      data x(177)/0.402270157963991d0/, a(177)/0.044590558163756d0/
      data x(178)/0.357220158337668d0/, a(178)/0.045491627927418d0/
      data x(179)/0.311322871990210d0/, a(179)/0.046284796581314d0/
      data x(180)/0.264687162208767d0/, a(180)/0.046968182816210d0/
      data x(181)/0.217423643740007d0/, a(181)/0.047540165714830d0/
      data x(182)/0.169644420423992d0/, a(182)/0.047999388596458d0/
      data x(183)/0.121462819296120d0/, a(183)/0.048344762234802d0/
      data x(184)/0.072993121787799d0/, a(184)/0.048575467441503d0/
      data x(185)/0.024350292663424d0/, a(185)/0.048690957009139d0/
c**** n=80
      data x(186)/0.999553822651630d0/, a(186)/0.001144950003186d0/
      data x(187)/0.997649864398237d0/, a(187)/0.002663533589512d0/
      data x(188)/0.994227540965688d0/, a(188)/0.004180313124694d0/
      data x(189)/0.989291302499755d0/, a(189)/0.005690922451403d0/
      data x(190)/0.982848572738629d0/, a(190)/0.007192904768117d0/
      data x(191)/0.974909140585727d0/, a(191)/0.008683945269260d0/
      data x(192)/0.965485089043799d0/, a(192)/0.010161766041103d0/
      data x(193)/0.954590766343634d0/, a(193)/0.011624114120797d0/
      data x(194)/0.942242761309872d0/, a(194)/0.013068761592401d0/
      data x(195)/0.928459877172445d0/, a(195)/0.014493508040509d0/
      data x(196)/0.913263102571757d0/, a(196)/0.015896183583725d0/
      data x(197)/0.896675579438770d0/, a(197)/0.017274652056269d0/
      data x(198)/0.878722567678213d0/, a(198)/0.018626814208299d0/
      data x(199)/0.859431406663111d0/, a(199)/0.019950610878141d0/
      data x(200)/0.838831473580255d0/, a(200)/0.021244026115782d0/
      data x(201)/0.816954138681463d0/, a(201)/0.022505090246332d0/
      data x(202)/0.793832717504605d0/, a(202)/0.023731882865930d0/
      data x(203)/0.769502420135041d0/, a(203)/0.024922535764115d0/
      data x(204)/0.744000297583597d0/, a(204)/0.026075235767565d0/
      data x(205)/0.717365185362099d0/, a(205)/0.027188227500486d0/
      data x(206)/0.689637644342027d0/, a(206)/0.028259816057276d0/
      data x(207)/0.660859898986119d0/, a(207)/0.029288369583267d0/
      data x(208)/0.631075773046871d0/, a(208)/0.030272321759557d0/
      data x(209)/0.600330622829751d0/, a(209)/0.031210174188114d0/
      data x(210)/0.568671268122709d0/, a(210)/0.032100498673487d0/
      data x(211)/0.536145920897131d0/, a(211)/0.032941939397645d0/
      data x(212)/0.502804111888784d0/, a(212)/0.033733214984611d0/
      data x(213)/0.468696615170544d0/, a(213)/0.034473120451753d0/
      data x(214)/0.433875370831756d0/, a(214)/0.035160529044747d0/
      data x(215)/0.398393405881969d0/, a(215)/0.035794393953416d0/
      data x(216)/0.362304753499487d0/, a(216)/0.036373749905835d0/
      data x(217)/0.325664370747701d0/, a(217)/0.036897714638276d0/
      data x(218)/0.288528054884511d0/, a(218)/0.037365490238730d0/
      data x(219)/0.250952358392272d0/, a(219)/0.037776364362001d0/
      data x(220)/0.212994502857666d0/, a(220)/0.038129711314477d0/
      data x(221)/0.174712291832646d0/, a(221)/0.038424993006959d0/
      data x(222)/0.136164022809143d0/, a(222)/0.038661759774076d0/
      data x(223)/0.097408398441584d0/, a(223)/0.038839651059051d0/
      data x(224)/0.058504437152420d0/, a(224)/0.038958395962769d0/
      data x(225)/0.019511383256793d0/, a(225)/0.039017813656306d0/
c**** n=96
      data x(226)/0.999689503883230d0/, a(226)/0.000796792065552d0/
      data x(227)/0.998364375863181d0/, a(227)/0.001853960788946d0/
      data x(228)/0.995981842987209d0/, a(228)/0.002910731817934d0/
      data x(229)/0.992543900323762d0/, a(229)/0.003964554338444d0/
      data x(230)/0.988054126329623d0/, a(230)/0.005014202742927d0/
      data x(231)/0.982517263563014d0/, a(231)/0.006058545504235d0/
      data x(232)/0.975939174585136d0/, a(232)/0.007096470791153d0/
      data x(233)/0.968326828463264d0/, a(233)/0.008126876925698d0/
      data x(234)/0.959688291448742d0/, a(234)/0.009148671230783d0/
      data x(235)/0.950032717784437d0/, a(235)/0.010160770535008d0/
      data x(236)/0.939370339752755d0/, a(236)/0.011162102099838d0/
      data x(237)/0.927712456722308d0/, a(237)/0.012151604671088d0/
      data x(238)/0.915071423120898d0/, a(238)/0.013128229566961d0/
      data x(239)/0.901460635315852d0/, a(239)/0.014090941772314d0/
      data x(240)/0.886894517402420d0/, a(240)/0.015038721026994d0/
      data x(241)/0.871388505909296d0/, a(241)/0.015970562902562d0/
      data x(242)/0.854959033434601d0/, a(242)/0.016885479864245d0/
      data x(243)/0.837623511228187d0/, a(243)/0.017782502316045d0/
      data x(244)/0.819400310737931d0/, a(244)/0.018660679627411d0/
      data x(245)/0.800308744139140d0/, a(245)/0.019519081140145d0/
      data x(246)/0.780369043867433d0/, a(246)/0.020356797154333d0/
      data x(247)/0.759602341176647d0/, a(247)/0.021172939892191d0/
      data x(248)/0.738030643744400d0/, a(248)/0.021966644438744d0/
      data x(249)/0.715676812348967d0/, a(249)/0.022737069658329d0/
      data x(250)/0.692564536642171d0/, a(250)/0.023483399085926d0/
      data x(251)/0.668718310043916d0/, a(251)/0.024204841792364d0/
      data x(252)/0.644163403784967d0/, a(252)/0.024900633222483d0/
      data x(253)/0.618925840125468d0/, a(253)/0.025570036005349d0/
      data x(254)/0.593032364777572d0/, a(254)/0.026212340735672d0/
      data x(255)/0.566510418561397d0/, a(255)/0.026826866725591d0/
      data x(256)/0.539388108324357d0/, a(256)/0.027412962726029d0/
      data x(257)/0.511694177154667d0/, a(257)/0.027970007616848d0/
      data x(258)/0.483457973920596d0/, a(258)/0.028497411065085d0/
      data x(259)/0.454709422167743d0/, a(259)/0.028994614150555d0/
      data x(260)/0.425478988407300d0/, a(260)/0.029461089958167d0/
      data x(261)/0.395797649828908d0/, a(261)/0.029896344136328d0/
      data x(262)/0.365696861472313d0/, a(262)/0.030299915420827d0/
      data x(263)/0.335208522892625d0/, a(263)/0.030671376123669d0/
      data x(264)/0.304364944354496d0/, a(264)/0.031010332586313d0/
      data x(265)/0.273198812591049d0/, a(265)/0.031316425596861d0/
      data x(266)/0.241743156163840d0/, a(266)/0.031589330770727d0/
      data x(267)/0.210031310460567d0/, a(267)/0.031828758894411d0/
      data x(268)/0.178096882367618d0/, a(268)/0.032034456231992d0/
      data x(269)/0.145973714654896d0/, a(269)/0.032206204794030d0/
      data x(270)/0.113695850110665d0/, a(270)/0.032343822568575d0/
      data x(271)/0.081297495464425d0/, a(271)/0.032447163714064d0/
      data x(272)/0.048812985136049d0/, a(272)/0.032516118713868d0/
      data x(273)/0.016276744849602d0/, a(273)/0.032550614492363d0/
c
c
c-----test n
      alpha=0.5d0*(ax+bx)
      beta=0.5d0*(bx-ax)
      if( n.lt.1 .or. n.gt.96 ) go to 100
      if(n.ne.1) go to 1
      z(1)=alpha
      w(1)=bx-ax
      return
c
    1 if (n.le.16) go to 3
      if (n.gt.24) go to 4
      n=4*(n/4)
      go to 3
    4 if (n.gt.48) go to 5
      n=8*(n/8)
      go to 3
    5 n=16*(n/16)
c
c----- set k equal to initial subscript and store results
    3 k=ktab(n)
      m=n/2
      do 2 j=1,m
      jtab=k-1+j
      wtemp=beta*a(jtab)
      delta=beta*x(jtab)
      z(j)=alpha-delta
      w(j)=wtemp
      jp=n+1-j
      z(jp)=alpha+delta
      w(jp)=wtemp
    2 continue
      if((n-m-m).eq.0) return
      z(m+1)=alpha
      jmid=k+m
      w(m+1)=beta*a(jmid)
      return
c
  100 zn=n
      write(6,200) zn
  200 format(/////' error in gset. n has the non-permissible value',
     1e11.3/' execution terminated.')
      stop
      end
c name:    dgelg
c        programmbibliothek rhrz bonn        02/02/81       dgelg
c                                            fortran iv     ibm 370/168
c
c purpose:
c
c to solve a general system of simultaneous linear equations.
c
c usage:   call dgelg(r,a,m,n,eps,ier)
c
c parameters:
c
c r:       double precision m by n right hand side matrix
c          (destroyed). on return r contains the solutions
c          of the equations.
c
c a:       double precision m by m coefficient matrix
c          (destroyed).
c
c m:       the number of equations in the system.
c
c n:       the number of right hand side vectors.
c
c eps:     single precision input constant which is used as
c          relative tolerance for test on loss of
c          significance.
c
c ier:     resulting error parameter coded as follows
c           ier=0  - no error,
c           ier=-1 - no result because of m less than 1 or
c                   pivot element at any elimination step
c                   equal to 0,
c           ier=k  - warning due to possible loss of signifi-
c                   cance indicated at elimination step k+1,
c                   where pivot element was less than or
c                   equal to the internal tolerance eps times
c                   absolutely greatest element of matrix a.
c
c remarks: (1) input matrices r and a are assumed to be stored
c              columnwise in m*n resp. m*m successive storage
c              locations. on return solution matrix r is stored
c              columnwise too.
c          (2) the procedure gives results if the number of equations m
c              is greater than 0 and pivot elements at all elimination
c              steps are different from 0. however warning ier=k - if
c              given indicates possible loss of significance. in case
c              of a well scaled matrix a and appropriate tolerance eps,
c              ier=k may be interpreted that matrix a has the rank k.
c              no warning is given in case m=1.
c
c method:
c
c solution is done by means of gauss-elimination with
c complete pivoting.
c
c programs required:
c          none
c
c access:
c
c load module:    sys3.fortlib(dgelg)
c source module:  sys3.symlib.fortran(dgelg)
c description:    sys3.infolib(dgelg)
c
c author:         ibm, ssp iii
c installation:   ibm 370/168, mvs-jes2, fortran iv (h ext. enh.)
c
c**********************************************************************
      subroutine dgelg(r,a,m,n,eps,ier)
c
c
      implicit real*8 (a-h,o-z)
      dimension a(1),r(1)
      real*4 eps
c
c
c
c
      if(m)23,23,1
c
c     search for greatest element in matrix a
    1 ier=0
      piv=0.d0
      mm=m*m
      nm=n*m
      do 3 l=1,mm
      tb=dabs(a(l))
      if(tb-piv)3,3,2
    2 piv=tb
      i=l
    3 continue
      tol=eps*piv
c     a(i) is pivot element. piv contains the absolute value of a(i).
c
c
c     start elimination loop
      lst=1
      do 17 k=1,m
c
c     test on singularity
      if(piv)23,23,4
    4 if(ier)7,5,7
    5 if(piv-tol)6,6,7
    6 ier=k-1
    7 pivi=1.d0/a(i)
      j=(i-1)/m
      i=i-j*m-k
      j=j+1-k
c     i+k is row-index, j+k column-index of pivot element
c
c     pivot row reduction and row interchange in right hand side r
      do 8 l=k,nm,m
      ll=l+i
      tb=pivi*r(ll)
      r(ll)=r(l)
    8 r(l)=tb
c
c     is elimination terminated
      if(k-m)9,18,18
c
c     column interchange in matrix a
    9 lend=lst+m-k
      if(j)12,12,10
   10 ii=j*m
      do 11 l=lst,lend
      tb=a(l)
      ll=l+ii
      a(l)=a(ll)
   11 a(ll)=tb
c
c     row interchange and pivot row reduction in matrix a
   12 do 13 l=lst,mm,m
      ll=l+i
      tb=pivi*a(ll)
      a(ll)=a(l)
   13 a(l)=tb
c
c     save column interchange information
      a(lst)=j
c
c     element reduction and next pivot search
      piv=0.d0
      lst=lst+1
      j=0
      do 16 ii=lst,lend
      pivi=-a(ii)
      ist=ii+m
      j=j+1
      do 15 l=ist,mm,m
      ll=l-j
      a(l)=a(l)+pivi*a(ll)
      tb=dabs(a(l))
      if(tb-piv)15,15,14
   14 piv=tb
      i=l
   15 continue
      do 16 l=k,nm,m
      ll=l+j
   16 r(ll)=r(ll)+pivi*r(l)
   17 lst=lst+m
c     end of elimination loop
c
c
c     back substitution and back interchange
   18 if(m-1)23,22,19
   19 ist=mm+m
      lst=m+1
      do 21 i=2,m
      ii=lst-i
      ist=ist-lst
      l=ist-m
      l=a(l)+.5d0
      do 21 j=ii,nm,m
      tb=r(j)
      ll=j
      do 20 k=ist,mm,m
      ll=ll+1
   20 tb=tb-a(k)*r(ll)
      k=j+l
      r(j)=r(k)
   21 r(k)=tb
   22 return
c
c
c     error return
   23 ier=-1
      return
      end