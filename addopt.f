c***********************************************************************
c    version one :March 2022
c    version two :June 2022 (change the representation of theparameters
c    in lsj representation)
c    version three: 2023.3.29. LO
c    version four: 2023.4.6 NLO 
c    version five:2023.4.7 N2LO
c    version six:2023.4.10 N2LO rebulid codes
c    version seven:2023.4.14 N3LO except two loop terms
c    version eight:2023.4.18 N3LO pion terms finished
c    version nine:2023.4.20 N3LO all
c
c***********************************************************************
c    
c    in this code you can add any terms you want in the chiral force
c    using module addterms,note that you have to use mass of nucleon as
c    the unit of energy
c    
c***********************************************************************  
c
c    author: shang haoyu 
c            physics school
c            Peking university
c    email:  shy@stu.pku.edu.cn
c
c
c***********************************************************************   
c    module kqxy contains the variable we use mostly in this code

c    function initialize evaluate the variables we use in this code

c    function normk, normq means the length of vector k and q
c    thier variable z means cos(theta)
      module const

c   variable:pi,ga,mpi,mpi0,mpipm,mass,fpi,tidelambda,c1,c2,c3,c4
c   subroutine :ini_const

         real*8 ::alpha=137.03599976d0
         real*8 ::pi=3.141592653589793d0
         real*8 ::ga=1.29d0
         real*8 ::mpi=138.0390d0
         real*8 ::mpi0=134.9766d0
         real*8 ::mpipm=139.5702d0
         real*8 ::mass=938.9182d0
         real*8 ::fpi=92.4d0
         real*8 ::tidelambda=650.0d0
         real*8 ::c1=-1.07d0
         real*8 ::c2=3.20d0
         real*8 ::c3=-5.32d0
         real*8 ::c4=3.56d0
         real*8 ::d15m14=1.90d0
         real*8 ::d1p2=1.04d0 
         real*8 ::d3=-0.48d0 
         real*8 ::d5=0.14d0
         contains
         subroutine ini_const
c      this subroutine should been used in the main programm            
            mpi=mpi/mass
            mpi0=mpi0/mass
            mpipm=mpipm/mass
            fpi=fpi/mass
            tidelambda=tidelambda/mass
            c1=c1*1.0d-3*mass
            c2=c2*1.0d-3*mass
            c3=c3*1.0d-3*mass
            c4=c4*1.0d-3*mass
            d15m14=d15m14*1.0d-6*mass*mass
            d1p2=d1p2*1.0d-6*mass*mass
            d3=d3*1.0d-6*mass*mass
            d5=d5*1.0d-6*mass*mass
         end subroutine
      end module
     
      module potential_global
c     this module contains the global variable of the subroutine potential

c     variables:xmev,ymev,conta(24),lambda,j,v(6)
c     v(6) means
c     in the following order:
c     0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled)
c     input
      real*8 xmev,ymev
      real*8 conta(24)
      real*8 lambda
      integer j
c     output
      real*8 v(6)
             
      end module
      module paravari

c variable:  pi,ga,mpi,mpi0,mpipm,mass,fpi,tidelambda,c1,c2,c3,c4(const)
c             x,y,dwn,wnq,wn3,dwnq,x2,y2,c(24)
c subroutine: ini_paravari      
       
        use const
        real*8,save :: x,y,dwn,wnq,wn3,dwnq,x2,y2,c(24)

        contains

          subroutine ini_paravari
            use potential_global,only:xmev,ymev,conta
c this subroutine should be used in potential subroutine
cas it contains the variables often been used             
            implicit real*8 (a-h,o-z)
            real*8 matrix1(9,9),matrix2(15,15)
            real*8 t(24)
            logical :: parlsj=.true.          
            data((matrix1(j,i),i=1,9),j=1,9)/
     1 0.0198943d0,       0.0d0,       0.0d0,       0.0d0,        0.0d0,
     * 0.0596831d0,       0.0d0,       0.0d0,       0.0d0,
     2       0.0d0, 0.0099471d0,-0.0049735d0,-0.0149207d0, -0.0149207d0,
     *       0.0d0, 0.0298415d0,       0.0d0, -0.024867d0,
     3       0.0d0, 0.0397887d0, 0.0198943d0, 0.0596831d0,  0.0596831d0,
     *       0.0d0, 0.1193662d0,       0.0d0, 0.0994718d0,
     4 -0.019894d0,       0.0d0,       0.0d0,       0.0d0,        0.0d0,
     * 0.0198943d0,       0.0d0,       0.0d0,       0.0d0,
     5       0.0d0,-0.0099471d0,-0.0049735d0, 0.0149207d0,        0.0d0,
     *       0.0d0, 0.0099471d0, 0.0140674d0,-0.0099471d0,
     6       0.0d0, -0.039788d0, 0.0198943d0,-0.0596831d0,        0.0d0,
     *       0.0d0, 0.0397887d0, 0.0562697d0, 0.0397887d0,
     7       0.0d0,       0.0d0,-0.0397887d0,       0.0d0, -0.0596831d0,
     *       0.0d0,       0.0d0,       0.0d0, 0.0994718d0,
     8       0.0d0,       0.0d0, 0.0099471d0,       0.0d0, -0.0149207d0,
     *       0.0d0,       0.0d0,-0.0422023d0, 0.0049735d0,
     9       0.0d0,       0.0d0,-0.0397887d0,       0.0d0,  0.0596831d0,
     *       0.0d0,       0.0d0,-0.1688093d0,-0.0198943d0/
c   here we evaluate the number we use
            data((matrix2(j,i),i=1,15),j=1,15)/
     1  0.002486796d0,	0.001243398d0,	-0.002486796d0,	-0.007460388d0,
     *  -0.007460388d0,	0.007460388d0,	0.003730194d0,	0.003730194d0,	
     *  0.0d0,	0.0d0,	0.00621699d0,	0.00621699d0,	-0.01243398d0,	
     *  0.0d0,	0.008703786d0,	
     2  0.039788736d0,	0.019894368d0,	0.039788736d0,	0.119366207d0,	
     *  0.119366207d0,	0.119366207d0,	0.059683104d0,	0.059683104d0,	
     *  0.0d0,	0.0d0,	0.099471839d0,	0.099471839d0,	0.198943679d0,	
     *  0.0d0,	0.139260575d0,	
     3  0.059683104d0,	-0.009947184d0,	0.0d0,	0.0d0,	0.0d0,	
     *  0.179049311d0,	-0.029841552d0,	-0.029841552d0,	0.0d0,
     *  0.0d0,	-0.04973592d0,	-0.04973592d0,	0.0d0,	0.0d0,	
     * -0.069630288d0,	
     4 -0.039788736d0,	0.019894368d0,	0.0d0,	0.0d0,	0.0d0,	
     * -0.119366207d0,	0.059683104d0,	-0.029841552d0,	0.0d0,	
     * 0.0d0,	-0.04973592d0,	-0.04973592d0,	0.0d0,	0.0d0,
     *	-0.069630288d0,	
     5 -0.002486796d0,	-0.001243398d0,	-0.002486796d0,
     * 0.007460388d0,	0.0d0,	0.002486796d0,	0.001243398d0,
     * 0.002486796d0,	0.003516861d0,	0.003516861d0,	-0.00621699d0,
     * 0.0d0,	-0.004973592d0,	-0.006091381d0,	0.003730194d0,	
     6 -0.039788736d0,	-0.019894368d0,	0.039788736d0,	
     * -0.119366207d0,	0.0d0,	0.039788736d0,	0.019894368d0,
     * 0.039788736d0,	0.05626977d0,	0.05626977d0,	-0.099471839d0,
     *	0.0d0,	0.079577472d0,	0.0974621d0,	0.059683104d0,	
     7 -0.059683104d0,	0.009947184d0,	0.0d0,	0.0d0,	0.0d0,	
     * 0.059683104d0,	-0.009947184d0,	-0.019894368d0,
     * 0.084404655d0,	-0.028134885d0,	0.04973592d0,	0.0d0,
     * 0.0d0,	0.0d0,	-0.029841552d0,	
     8 0.039788736d0,	-0.019894368d0,	0.0d0,	0.0d0,	
     * 0.0d0,	-0.039788736d0,	0.019894368d0,	0.009947184d0,	
     * 0.028134885d0,	-0.028134885d0,	0.04973592d0,	-0.04973592d0,
     * 0.0d0,	0.0d0,	-0.009947184d0,	
     9 0.0d0,	0.0d0,	-0.019894368d0,	0.0d0,	-0.029841552d0,	
     * 0.0d0,	0.0d0,	0.044762328d0,	0.0d0,	0.0d0,	0.0d0,
     * 0.02486796d0,	0.04973592d0,	0.0d0,	-0.069630288d0,	
     * 0.0d0,	0.0d0,	-0.079577472d0,	0.0d0,	-0.119366207d0,
     * 0.0d0,	0.0d0,	-0.179049311d0,	0.0d0,	0.0d0,	0.0d0,
     *	-0.099471839d0,	0.198943679d0,	0.0d0,	0.27852115d0,	
     1 0.0d0,	0.0d0,	0.004973592d0,	0.0d0,	-0.007460388d0,	
     * 0.0d0,	0.0d0,	-0.003730194d0,	-0.010550582d0,	
     * -0.010550582d0,	0.0d0,	0.00621699d0,	0.002486796d0,
     * 0.018274144d0,	-0.002486796d0,	
     2 0.0d0,	0.0d0,	0.019894368d0,	0.0d0,	-0.029841552d0,	
     * 0.0d0,	0.0d0,	0.014920776d0,	-0.126606982d0,	0.042202327d0,
     * 0.0d0,	-0.02486796d0,	0.009947184d0,	-0.170558675d0,
     * 0.009947184d0,	
     3 0.0d0,	0.0d0,	-0.019894368d0,	0.0d0,	0.029841552d0,
     * 0.0d0,	0.0d0,	0.014920776d0,	-0.126606982d0,	0.042202327d0,
     * 0.0d0,	-0.02486796d0,	-0.009947184d0,	0.170558675d0,
     * 0.009947184d0,	
     4 0.0d0,	0.0d0,	-0.079577472d0,	0.0d0,	0.119366207d0,	
     * 0.0d0, 0.0d0,	-0.059683104d0,	-0.168809309d0,	
     * -0.168809309d0,0.0d0,	0.099471839d0,	-0.039788736d0,
     * -0.2923863d0,	-0.039788736d0,	
     5 0.0d0,	0.0d0,	0.0d0,	0.0d0,	0.0d0,	0.0d0,	0.0d0,
     *	-0.059683104d0,	-0.084404655d0,	0.084404655d0,	0.0d0,	
     * 0.099471839d0,	0.0d0,	0.0d0,	-0.039788736d0/
   
                     
            dwn=1.0d0/mass
            wnq=mass*mass
            wn3=mass*mass*mass
            dwnq=dwn*dwn

            x=xmev*dwn
            y=ymev*dwn
            x2=x*x
            y2=y*y
c    the initial of contact parameters:if parlsj is .true.,using 
c    parameters in lsj;if it's .false.,using parameters in pphase           
            c=0.0d0
            c=conta
            if(parlsj)then
            t(1)=c(1)
            t(2)=c(2)
            t(10)=c(3)
            t(11)=c(4)
            t(3)=c(5)
            t(12)=c(6)
            t(4)=c(7)
            t(13)=c(8)
            t(5)=c(9)
            t(14)=c(10)
            t(6)=c(11)
            t(7)=c(12)
            t(15)=c(13)
            t(16)=c(14)
            t(17)=c(15)
            t(8)=c(16)
            t(18)=c(17)
            t(19)=c(18)
            t(20)=c(19)
            t(21)=c(20)
            t(9)=c(21)
            t(22)=c(22)
            t(23)=c(23)
            t(24)=c(24)
            c=0.0d0
            do i=1,9
            do j=1,9
               c(i)=c(i)+matrix1(i,j)*t(j)
            end do 
            end do
            do i=10,24
            do j=10,24
               c(i)=c(i)+matrix2(i-9,j-9)*t(j)
            end do
            end do
         else
            t(1)=c(2)
            t(2)=c(3)
            t(3)=c(4)
            c(2)=t(2)
            c(3)=t(3)
            c(4)=t(1)            
         end if
         c(1)=c(1)*0.01d0*wnq
         c(4)=c(4)*0.01d0*wnq
         c(2)=c(2)*wnq*wnq*1.0d-8
         c(3)=c(3)*wnq*wnq*1.0d-8
         c(5:9)=c(5:9)*wnq*wnq*1.0d-8
         c(10:24)=c(10:24)*wn3*wn3*1.0d-14
          end subroutine
         end module
         
         module genfunc

c variable: none
c function:normk,normq,qdotk,kcrossq2,wfunc,lfunc,afunc
            
            use paravari
            private pi,ga,mpi,mpi0,mpipm,mass,fpi,tidelambda,
     +       c1,c2,c3,c4,x,y,dwn,wnq,wn3,dwnq,x2,y2,c
c this module contains general functions will be used in
c potential subroutine
         contains
          real*8 function normk(z)
            implicit real*8 (a-h,o-z)
            real*8 z
            normk=dsqrt(x*x+y*y+2.0d0*x*y*z)/2.0d0
            return
          end function
        
          real*8 function normq(z)
            implicit real*8 (a-h,o-z)
            real*8 z    
            normq=dsqrt(x*x+y*y-2.0d0*x*y*z)
            return
          end function

          real*8 function qdotk()
            implicit real*8 (a-h,o-z)
            qdotk=(x*x-y*y)/2.0d0
            return
          end function

          real*8 function kcrossq2(z)
            implicit real*8 (a-h,o-z)
            real*8 z
            kcrossq2=normq(z)**2*normk(z)**2
     1      -(qdotk())**2
            return
          end function

c   function w(q) in (20) PRC 66,014002(2002)  
         
          real*8 function wfunc(z)
             implicit real*8 (a-h,o-z)
             real*8 z    
             wfunc=dsqrt(4.0d0*(mpi)**2+normq(z)**2)
             return
          end function
c   function L(q) in (19) PRC 66,014002(2002)         
          real*8 function lfunc(z)
            implicit real*8 (a-h,o-z)
            real*8 z  
            logical :: SFC=.true.
c   SFC means spectrum function cutoff            
            if (SFC)then
            lfunc=wfunc(z)/(2.0d0*normq(z))*dlog(((tidelambda)**2
     1       *(2.0d0*(mpi)**2+normq(z)**2)-2.0d0*(mpi)**2
     2      *normq(z)**2+tidelambda*dsqrt((tidelambda)**2
     3      -4.0d0*(mpi)**2)*normq(z)*wfunc(z))
     4      /(2.0d0*(mpi)**2*((tidelambda)**2+normq(z)**2)))
            else  
            lfunc=wfunc(z)/normq(z)*dlog((wfunc(z)+normq(z))
     1       /(2.0d0*mpi))
            end if 
            return
          end function
         
          real*8 function afunc(z)
          real*8 z
           logical :: SFC=.true.
           if (SFC)then 
            afunc=datan((normq(z)*(tidelambda-2.0d0*mpi))/
     1       (normq(z)**2+2.0d0*tidelambda*mpi))
     2      /(2.0d0*normq(z))
           else
c   tidelambda = infinity 
            afunc=datan(normq(z)/(2.0d0*mpi))/(2.0d0*normq(z))
           end if 
           return
           end function
        end module

c    in module add-terms we can add potential terms arbitrary
      
      module addterms
        use genfunc
        use paravari
        use velementm
c    we set some logical variable to control whether the terms 
c    are contained(private),true means they are contained
        logical :: vc=.true.
        logical :: vss=.true.
        logical :: vso=.true.
        logical :: vsigl=.true.
        logical :: vt=.true.
        logical :: vsigk=.true.

        type(velement) :: nlo_ct
        type(velement) ::nlo_tp !two pi


        private vc,vt,vss,vso,vsigl,vsigk
c    in the following terms w means terms with tau1*tau2
c    z is cos(theta) 

!!!! notice here we use mass of nucleon as the uint of energy

        contains
          real*8 function vcentral(z)
            real*8 z
            real*8 xr(6)
            real*8 a00,a10,a01,a11,b01,b11,b10
c   using type to determine the pade form we use             
            integer type
            xr(1)=0.01d0*wnq*c(1)
            xr(2)=c(2)*wnq*wnq*1.0d-8
            xr(3)=c(3)*wnq*wnq*1.0d-8
            xr(4)=c(10)*wn3*wn3*1.0d-14
            xr(5)=c(11)*wn3*wn3*1.0d-14
            xr(6)=c(12)*wn3*wn3*1.0d-14
            b01=-xr(5)/xr(3)
            b10=-xr(4)/xr(2)
            b11=(-b01*(xr(4)+xr(6))-b10*(xr(5)+xr(6)))
            a00=xr(1)
            a10=xr(2)+b10*xr(1)
            a01=xr(3)+b01*xr(1)
            a11=xr(6)+b01*xr(2)+b10*xr(3)+b11*xr(1)
            type=1
            if(vc) then
               select case(type)
c type 1 means original n3lo
               case(1)
               vcentral=0.01d0*wnq*c(1)+c(2)*wnq*wnq*1.0d-8*normq(z)**2
     1      +c(3)*wnq*wnq*1.0d-8*normk(z)**2
     2      +c(10)*wn3*wn3*1.0d-14*normq(z)**4+c(11)*wn3*wn3*1.0d-14
     3      *normk(z)**4+c(12)*wn3*wn3*1.0d-14*normq(z)**2*normk(z)**2
     4      +c(13)*kcrossq2(z)*wn3*wn3*1.0d-14
               case(2)
c type 2 using single variable form pade                   
              vcentral=0.01d0*wnq*c(1)+(c(2)*wnq*wnq*1.0d-8*normq(z)**2
     1      +c(3)*wnq*wnq*1.0d-8*normk(z)**2)/(1.0d0-
     2      (c(10)*wn3*wn3*1.0d-14*normq(z)**4+c(11)*wn3*wn3*1.0d-14
     3      *normk(z)**4+c(12)*wn3*wn3*1.0d-14*normq(z)**2*normk(z)**2
     4      )/(c(2)*wnq*wnq*1.0d-8*normq(z)**2+c(3)*wnq*wnq*1.0d-8
     5      *normk(z)**2)) +c(13)*kcrossq2(z)*wn3*wn3*1.0d-14
               case(3)
c type 3 using 2 variable form pade (1,1)/(1,1)
            vcentral=(a00+a10*normq(z)**2+a01*normk(z)**2+a11*
     +      normq(z)**2*normk(z)**2)/(1.0d0+b10*normq(z)**2+b01
     +      *normk(z)**2+a11*normq(z)**2*normk(z)**2)+c(13)
     +      *kcrossq2(z)*wn3*wn3*1.0d-14
            end select
            else
            vcentral=0.0d0
            end if
          return
          end function

          real*8 function nlowc(z)
          real*8 z
          nlowc=-lfunc(z)/(384.0d0*pi**2*(fpi)**4)
     +    *(4.0d0*(mpi)**2*(5.0d0*ga**4-4.0d0*ga**2-1.0d0)
     +    +normq(z)**2*(23.0d0*ga**4-10.0d0*ga**2-1.0d0)
     +    +48.0d0*ga**4*(mpi)**4/wfunc(z)**2)
          return
         end function

          real*8 function n2lovc(z)
          real*8 z
          n2lovc=3.0d0*ga**2/(16.0d0*pi*(fpi)**4)
     1     *(2.0d0*(mpi)**2*(c3-2.0d0*c1)
     2     +c3*normq(z)**2)*(2.0d0*(mpi)**2
     3     +normq(z)**2)*afunc(z)
          return
      end function

          real*8 function vspinspin(z)
          real*8 z,ct,c3,d5,q2
          integer :: type=1
          ct=0.01d0*wnq*c(4)
          c3=c(5)*wnq*wnq*1.0d-8
          d5=c(14)*wn3*wn3*1.0d-14
          q2=normq(z)**2
          if(vss) then
            select case(type)
            case(1)
c        the original n3lo               
          vspinspin=0.01d0*wnq*c(4)+c(5)*wnq*wnq*1.0d-8*normq(z)**2
     1      +c(6)*wnq*wnq*1.0d-8*normk(z)**2 +c(14)*wn3*wn3*1.0d-14
     2      *normq(z)**4+c(15)*wn3*wn3*1.0d-14*normk(z)**4+c(16)
     3      *wn3*wn3*1.0d-14*normq(z)**2*normk(z)**2+c(17)
     4      *kcrossq2(z)*wn3*wn3*1.0d-14
            case(2)
c        mix q terms using pade approximant
           vspinspin=(ct*c3+c3**2*q2-ct*d5*q2)/(c3-d5*q2)
     1    + c(6)*wnq*wnq*1.0d-8*normk(z)**2
     2    + c(15)*wn3*wn3*1.0d-14*normk(z)**4
     3    + c(16)*wn3*wn3*1.0d-14*normq(z)**2*normk(z)**2
     4    + c(17)*kcrossq2(z)*wn3*wn3*1.0d-14
            case(3)
            vspinspin=ct/(1-c3*q2/ct)+d5*q2*2
     1    + c(6)*wnq*wnq*1.0d-8*normk(z)**2
     2    + c(15)*wn3*wn3*1.0d-14*normk(z)**4
     3    + c(16)*wn3*wn3*1.0d-14*normq(z)**2*normk(z)**2
     4    + c(17)*kcrossq2(z)*wn3*wn3*1.0d-14   
            end select
          else
          vspinspin=0.0d0
          end if
          return
          end function

          real*8 function nlovss(z)
          real*8 z
          nlovss=-normq(z)**2*nlovt(z)
          return
         end function
           
          real*8 function n2lows(z)
          real*8 z
          n2lows=-normq(z)**2*n2lowt(z)
          return
          end function

          real*8 function vspinobit(z)
          real*8 z,c5,d9,d10,k2
          integer ::type=1
          c5=c(7)*wnq*wnq*1.0d-8
          d9=c(18)*wn3*wn3*1.0d-14
          d10=c(19)*wn3*wn3*1.0d-14
          k2=normk(z)**2
          if(vso) then
            select case(type)
            case(1)
          vspinobit=c(7)*wnq*wnq*1.0d-8+c(18)*wn3*wn3*1.0d-14
     1     *normq(z)**2+c(19)*wn3*wn3*1.0d-14*normk(z)**2
            case(2)
          vspinobit=c5/(1.0d0-d10*k2/c5)+d9*normq(z)**2
          end select
          else
          vspinobit=0.0d0
          end if
          return
          end function

          real*8 function vsigmaL(z)
          real*8 z
          if(vsigl) then
          vsigmaL=c(24)*wn3*wn3*1.0d-14
          else
          vsigmaL=0.0d0
          end if
          return
          end function

          real*8 function vtensor(z)
            real*8 z
            integer type
            type=1
            if(vt) then
                select case(type)
c  the original case (n3lo)
                case(1)
            vtensor=c(8)*wnq*wnq*1.0d-8+c(20)*wn3*wn3*1.0d-14
     +     *normq(z)**2+c(21)*wn3*wn3*1.0d-14*normk(z)**2

c the nlo       
c  yukawwa form a/(1+b*q^2) + c*k^2
                case(2)
                    vtensor=c(8)*wnq*wnq*1.0d-8/(1.0d0-c(20)*wn3*wn3
     +              *1.0d-14*normq(z)**2/c(8)*wnq*wnq*1.0d-8)+c(21)
     +              *wn3*wn3*1.0d-14*normk(z)**2
c  form (a+b*k^2)/(c+d*q^2)
                case(3)
                    vtensor=(c(8)*wnq*wnq*1.0d-8+c(21)*wn3*wn3
     +              *1.0d-14*normk(z)**2)/(1.0d0-c(20)*wn3*wn3
     +              *1.0d-14*normq(z)**2/c(8)*wnq*wnq*1.0d-8)
c  form 
                end select    
            else
            vtensor=0.0d0
            end if
          return
          end function
c   nlo vt
         
           
          real*8 function nlovt(z)
          real*8 z
          nlovt=-3.0d0*ga**4/(64.0d0*pi**2*(fpi)**4)*lfunc(z)
          return
         end function
         


          real*8 function nlo_vc_ct(z)
          real*8 z
          nlo_vc_ct=c(2)*normq(z)**2+c(3)*normk(z)**2
          return
          end function

          real*8 function nlo_vs_ct(z)
          real*8 z
          nlo_vs_ct=c(5)*normq(z)**2+c(6)*normk(z)**2 
          return
         end function

          real*8 function nlo_vls_ct(z)
          real*8 z
          nlo_vls_ct=c(7)
          return
          end function
         
          real*8 function nlo_vt_ct(z)
          real*8 z
          nlo_vt_ct=c(8)
          return
         end function

          real*8 function nlo_vsk_ct(z)
          real*8 z
          nlo_vsk_ct=c(9)
          return
          end function

          real*8 function n2lowt(z)
          real*8 z
          n2lowt=-ga**2/(32.0d0*pi*(fpi)**4)
     1   *c4*wfunc(z)**2*afunc(z)
          return
          end function
          

          real*8 function vsigmak(z)
          real*8 z
          if(vsigk) then
          vsigmak=c(9)*wnq*wnq*1.0d-8+c(22)*wn3*wn3*1.0d-14
     1     *normq(z)**2+c(23)*wn3*wn3*1.0d-14*normk(z)**2
          else
          vsigmak=0.0d0
          end if
          return
          end function
      end module 
      
      module velementm
         private
         public velement
         type velement
         real*8 ::vc(6)
         real*8 ::wc(6)
         real*8 ::vss(6)
         real*8 ::wss(6)
         real*8 ::vt(6)
         real*8 ::wt(6)
         real*8 ::vls(6)
         real*8 ::wls(6)
         real*8 ::vsk(6)
         real*8 ::vslsl(6)
         real*8 ::sum(6)
         contains 
         procedure :: init => ini_velement
         procedure :: add => add_velement
         end type
         contains
         subroutine ini_velement(this)
            class(velement) ::this
            this%vc=0.0d0
            this%vls=0.0d0
            this%vss=0.0d0
            this%vt=0.0d0
            this%wc=0.0d0
            this%wls=0.0d0
            this%wss=0.0d0
            this%wt=0.0d0
            this%vsk=0.0d0
            this%vslsl=0.0d0
            this%sum=0.0d0
            end subroutine
         subroutine add_velement(this)
            class(velement) ::this
            this%sum=this%vc+this%vls+this%vss+this%vt
     +    +this%wc+this%wls+this%wss+this%wt+this%vsk+this%vslsl
         end subroutine

      end module

      module lopot 
         use genfunc
         use paravari,only:c
         use velementm
         use const
         
         implicit none 
         private

c    the interface
c    variable
        public lo_ct,lo_onepi 
c    function
        public  onepii0,onepii1
        public lo_vc_ct,lo_vs_ct      

        type(velement) ::lo_ct
        type(velement) ::lo_onepi
        contains
c    one-pi exchange term         
         real*8 function onepi(z,masspi)
          real*8 z,masspi
          onepi=-ga**2/(4.0d0*(fpi)**2
     +     *(normq(z)**2+(masspi)**2))
          return
          end function

          real*8 function onepi1(z)
          real*8 z
          onepi1=onepi(z,mpi)
          return
         end function
c     charge dependent one-pi terms
c         I=0
          real*8 function onepii0(z)
          real*8 z
          onepii0=-onepi(z,mpi0)-2.0d0*onepi(z,mpipm)
          return
          end function
c         I=1
          real*8 function onepii1(z)
          real*8 z
          onepii1=-onepi(z,mpi0)+2.0d0*onepi(z,mpipm)
          return
          end function

c       contact terms
          real*8 function lo_vc_ct(z)
          real*8 z
          lo_vc_ct=c(1)
          return
          end function

          real*8 function lo_vs_ct(z)
          real*8 z
          lo_vs_ct=c(4)
          return
          end function          
      
      end module      
      module n3lopot
c to improve the speed of calculation,you should use 
c n3lo_ini before using subroutine n3lo500new(in the main programm)                 
         use velementm
         use const
         use genfunc
         use paravari,only:c

         implicit none 
         private

c     the interface
c     variable
         public n3lo_fd,n3lo_rc,n3lo_tl,n3lo_cM,n3lo_ct
c     function
         public n3lo_vc_fd,n3lo_wt_fd,n3lo_ws_fd 
         public n3lo_vc_rc,n3lo_wc_rc,n3lo_vt_rc,n3lo_wt_rc,
     +    n3lo_vs_rc,n3lo_ws_rc,n3lo_vls_rc,n3lo_wls_rc
         public n3lo_vc_tl,n3lo_ws_tl,n3lo_wt_tl,n3lo_vs_tl,n3lo_vt_tl,
     +    n3lo_wc_tl
         public n3lo_vc_cM,n3lo_wc_cM,n3lo_wt_cM,n3lo_ws_cM,
     +    n3lo_vls_cM,n3lo_wls_cM
         public n3lo_vc_ct,n3lo_vt_ct,n3lo_vls_ct,n3lo_vs_ct,
     +    n3lo_vsk_ct,n3lo_vslsl_ct
         public  pi_gamma  
c     subroutine         
         public n3lo_ini

         type(velement) ::n3lo_fd
         type(velement) ::n3lo_rc
         type(velement) ::n3lo_tl
         type(velement) ::n3lo_cM
         type(velement) ::n3lo_ct
         real*8 :: imvs(96),imwc(96),imvc(96),
     +   imws(96),imwt(96),imvt(96)    

         contains
         subroutine n3lo_ini
            implicit none
            real*8 wt(96),ct(96)
            real*8 xlb
            integer i
            xlb=2.0d0*mpi
            call fset(ct,wt,xlb,tidelambda,96)
            do i=1,96
            imvs(i)=n3lo_imvs(ct(i))
            imwc(i)=n3lo_imwc(ct(i))
            imvc(i)=n3lo_imvc(ct(i))
            imws(i)=n3lo_imws(ct(i))
            imwt(i)=n3lo_imwt(ct(i))
            imvt(i)=n3lo_imvt(ct(i))            
            end do
         end subroutine
c     football diagram(the 'fd')
         real*8 function n3lo_vc_fd(z)
         real*8 z
      n3lo_vc_fd=3.0d0/(16.0d0*pi**2*fpi**4)*((c2/6.0d0*wfunc(z)**2
     +   +c3*(2.0d0*mpi**2+normq(z)**2)-4.0d0*c1*mpi**2)**2+c2**2
     +   /45.0d0*wfunc(z)**4)*lfunc(z)
        return
        end function

        real*8 function n3lo_wt_fd(z)
        real*8 z
        n3lo_wt_fd=c4**2/(96.0d0*pi**2*fpi**4)*wfunc(z)**2*lfunc(z)
        return
      end function

        real*8 function n3lo_ws_fd(z)
        real*8 z
        n3lo_ws_fd=-normq(z)**2*n3lo_wt_fd(z)
        return
        end function
c     relativistic corrections(the 'rc')
        
        real*8 function n3lo_vc_rc(z)
        real*8 z
        n3lo_vc_rc=3.0d0*ga**4/(128.0d0*pi*fpi**4)
     +  *(mpi**5/(2.0d0*wfunc(z)**2)+(2.0d0*mpi**2
     +  +normq(z)**2)*(normq(z)**2-mpi**2)*afunc(z))
        return
      end function
      
        real*8 function n3lo_wc_rc(z)
        real*8 z
        n3lo_wc_rc=ga**2/(64.0d0*pi*fpi**4)*(3*ga**2
     +  *mpi**5/(2.0d0*wfunc(z)**2)+(ga**2*(3.0d0*mpi**2
     + +2.0d0*normq(z)**2)-2.0d0*mpi**2-normq(z)**2)
     +  *(2.0d0*mpi**2+normq(z)**2)*afunc(z))
        return
       end function
      
        real*8 function n3lo_vt_rc(z)
        real*8 z
        n3lo_vt_rc=3.0d0*ga**4/(256.0d0*pi*fpi**4)
     +   *(5.0d0*mpi**2+2.0d0*normq(z)**2)*afunc(z)
        return
      end function

        real*8 function n3lo_wt_rc(z)
        real*8 z
        n3lo_wt_rc=ga**2/(128.0d0*pi*fpi**4)
     +  *(ga**2*(3.0d0*mpi**2+normq(z)**2)-wfunc(z)**2)
     +  *afunc(z)
        return 
        end function

        real*8 function n3lo_vs_rc(z)
        real*8 z
        n3lo_vs_rc=-normq(z)**2*n3lo_vt_rc(z)
        return
      end function

        real*8 function n3lo_ws_rc(z)
        real*8 z
        n3lo_ws_rc=-normq(z)**2*n3lo_wt_rc(z)
        return
        end function

        real*8 function n3lo_vls_rc(z)
        real*8 z
        n3lo_vls_rc=3.0d0*ga**4/(32.0d0*pi*fpi**4)
     +  *(2.0d0*mpi**2+normq(z)**2)*afunc(z)
        return
      end function

        real*8 function n3lo_wls_rc(z)
        real*8 z
        n3lo_wls_rc=ga**2*(1.0d0-ga**2)/(32.0d0*pi*
     +   fpi**4)*wfunc(z)**2*afunc(z)
        return
       end function

c     spectral functions,name as n3lo_imv(mu)
      
        real*8 function n3lo_imvc(mu)
        real*8 mu
        n3lo_imvc=3*ga**4*(2.0d0*mpi**2-mu**2)/(pi*mu
     +  *(4.0d0*fpi)**6)*((mpi**2-2.0d0*mu**2)*(2.0d0*mpi
     +  +(2.0d0*mpi**2-mu**2)/(2.0d0*mu)*dlog((mu+2.0d0*mpi)
     + /(mu-2.0d0*mpi)))+4.0d0*ga**2*mpi*(2.0d0*mpi**2-mu**2))
        return
      end function

        real*8 function n3lo_imws(mu)
        real*8 mu
        n3lo_imws=ga**4*(4*mpi**2-mu**2)/(pi*(4.0d0*fpi)**6)
     +   *((mpi**2-mu**2/4.0d0)*dlog((mu+2.0d0*mpi)/(mu-2.0d0*mpi))
     +  +(1.0d0+2.0d0*ga**2)*mu*mpi)
        return
        end function

        real*8 function n3lo_imwt(mu)
        real*8 mu
        n3lo_imwt=n3lo_imws(mu)/mu**2
        return
      end function

        real*8 function n3lo_imvs(mu)
        implicit none
        real*8 mu
        real*8 rk  !rootk
        real*8 integral,xr  !x regulation
        real*8 wt(96),ct(96)
        integer i
        call fset(ct,wt,0.0d0,1.0d0,12)
        integral=0.0d0
        rk=dsqrt(mu**2/4.0d0-mpi**2)
        do i=1,12
        xr=ct(i)
        integral=integral+wt(i)*((1.0d0-xr**2)*(1.0d0/6.0d0
     +  -mpi**2/(rk**2*xr**2)+(1.0d0+mpi**2/(rk**2*xr**2))**(1.5d0)
     +   *dlog((rk*xr+dsqrt(mpi**2+rk**2*xr**2))/mpi)))  
        end do
        n3lo_imvs=ga**2*mu*rk**3/(8.0d0*pi*fpi**4)*d15m14
     +  +2.0d0*ga**6*mu*rk**3/(8.0d0*pi*fpi**2)**3*integral
        return
        end function

        real*8 function n3lo_imvt(mu)
        real*8 mu
        n3lo_imvt=n3lo_imvs(mu)/mu**2
        return
      end function
        
        real*8 function n3lo_imwc(mu)
        implicit none
        real*8 mu
        real*8 rk,wt(96),ct(96),xr
        real*8 term1,term21,term22,term23,term24,term
        real*8 integral
        integer ::i
        integral=0.0d0
        call fset(ct,wt,0.0d0,1.0d0,12)
        rk=dsqrt(mu**2/4.0d0-mpi**2)
        do i=1,12
        xr=ct(i)
        term1=ga**2*(mu**2-2.0d0*mpi**2)+2.0d0*
     +   (1.0d0-ga**2)*rk**2*xr**2
        term21=96.0d0*pi**2*fpi**2*((2.0d0*mpi**2-mu**2)
     +   *d1p2-2.0d0*rk**2*xr**2*d3+4.0d0*mpi**2*d5)
        term22=(4.0d0*mpi**2*(1.0d0+2.0d0*ga**2)-mu**2
     +   *(1.0d0+5.0d0*ga**2))*rk/mu*dlog((mu+2.0d0*rk)
     +  /(2.0d0*mpi))+mu**2/12.0d0*(5.0d0+13.0d0*ga**2)
     +  -2.0d0*mpi**2*(1.0d0+2.0*ga**2)
        term23=-3.0d0*rk**2*xr**2+6.0d0*rk*xr*dsqrt(mpi**2
     +   +rk**2*xr**2)*dlog((rk*xr+dsqrt(mpi**2+rk**2*xr**2))/mpi)
        term24=ga**4*(mu**2-2.0d0*rk**2*xr**2-2.0d0*mpi**2)
     +   *(5.0d0/6.0d0+mpi**2/(rk**2*xr**2)-(1.0d0+mpi**2
     +   /(rk**2*xr**2))**(1.5d0)*dlog((rk*xr+dsqrt(mpi**2
     +   +rk**2*xr**2))/mpi))
        term=term1*(term21+term22+term23+term24)
        integral=integral+wt(i)*term
        end do
        n3lo_imwc=(2.0d0*rk)/(3.0d0*mu*(8.0d0*pi*fpi**2)**3)*integral
        return
        end function
         
c     two loop contributions(tl)

c     v_{C,S}=-2q^2/pi*integrate_{2mpi,Lambda}(imvcs(i mu)/(mu^5(mu^2+q^2)))
        real*8 function vcstl(func,z)
c      input
        real*8 z
        real*8 func(96)
c     local 
        real*8 wt(96),ct(96)
        real*8 xlb
        real*8 integral
        integer i
        xlb=2.0d0*mpi
        call fset(ct,wt,xlb,tidelambda,96)
c     the integral
        integral=0.0d0
        do i=1,96
         integral=integral+wt(i)*(func(i)/
     +    (ct(i)**5*(ct(i)**2+normq(z)**2)))
        end do
        vcstl=-2.0d0*normq(z)**6/pi*integral
        return
      end function

c     v_{T,LS}
        real*8 function vtlstl(func,z)
c      input
        real*8 z
        real*8 func(96)
c     local 
        real*8 wt(96),ct(96)
        real*8 xlb
        real*8 integral
        integer i
        xlb=2.0d0*mpi
        call fset(ct,wt,xlb,tidelambda,96)
c     the integral
        integral=0.0d0
        do i=1,96
         integral=integral+wt(i)*(func(i)/
     +    (ct(i)**3*(ct(i)**2+normq(z)**2)))
        end do
        vtlstl=2.0d0*normq(z)**4/pi*integral
        return
      end function

      real*8 function n3lo_vc_tl(z)
      real*8 z
      n3lo_vc_tl=vcstl(imvc,z)
      return
      end function
      
      real*8 function n3lo_ws_tl(z)
      real*8 z
      n3lo_ws_tl=vcstl(imws,z)
      return
      end function 

      real*8 function n3lo_wt_tl(z)
      real*8 z
      n3lo_wt_tl=vtlstl(imwt,z)
      return
      end function

      real*8 function n3lo_vs_tl(z)
      real*8 z
      n3lo_vs_tl=vcstl(imvs,z)
      return
      end function 

      real*8 function n3lo_vt_tl(z)
      real*8 z
      n3lo_vt_tl=vtlstl(imvt,z)
      return
      end function

      real*8 function n3lo_wc_tl(z)
      real*8 z
      n3lo_wc_tl=vcstl(imwc,z)
      return
      end function

c    ci/M contributions('cM')
      real*8 function n3lo_vc_cM(z)
      real*8 z
      n3lo_vc_cM=-ga**2*lfunc(z)/(32.0d0*pi**2*fpi**4)
     +  *((c2-6.0d0*c3)*normq(z)**4+4.0d0*(6.0d0*c1+c2-3.0d0*c3)
     +  *normq(z)**2*mpi**2+6.0d0*(c2-2.0d0*c3)*mpi**4
     +  +24.0d0*(2.0d0*c1+c3)*mpi**6/wfunc(z)**2)
      return
      end function

      real*8 function n3lo_wc_cM(z)
      real*8 z
      n3lo_wc_cM=-c4*normq(z)**2*lfunc(z)/(192.0d0*pi**2*fpi**4)
     + *(ga**2*(8.0d0*mpi**2+5.0d0*normq(z)**2)+wfunc(z)**2)
      return
      end function
      
      real*8 function n3lo_wt_cM(z)
      real*8 z
      n3lo_wt_cM=-c4*lfunc(z)/(192.0d0*pi**2*fpi**4)*
     + (ga**2*(16.0d0*mpi**2+7.0d0*normq(z)**2)-wfunc(z)**2)
      return
      end function

      real*8 function n3lo_ws_cM(z)
      real*8 z
      n3lo_ws_cM=-n3lo_wt_cM(z)*normq(z)**2
      return
      end function

      real*8 function n3lo_vls_cM(z)
      real*8 z
      n3lo_vls_cM=c2*ga**2/(8.0d0*pi**2*fpi**4)
     + *wfunc(z)**2*lfunc(z)
      return
      end function

      real*8 function n3lo_wls_cM(z)
      real*8 z
      n3lo_wls_cM=-c4*lfunc(z)/(48.0d0*pi**2*fpi**4)*
     + (ga**2*(8.0d0*mpi**2+5.0d0*normq(z)**2)+wfunc(z)**2)
      return
      end function

c   contact terms(ct)
      real*8 function n3lo_vc_ct(z)
      real*8 z
       n3lo_vc_ct=c(10)*normq(z)**4+c(11)*normk(z)**4
     + +c(12)*normq(z)**2*normk(z)**2+c(13)*kcrossq2(z)
       return
      end function

      real*8 function n3lo_vs_ct(z)
      real*8 z
      n3lo_vs_ct=c(14)*normq(z)**4+c(15)*normk(z)**4
     + +c(16)*normq(z)**2*normk(z)**2+c(17)*kcrossq2(z)
       return
       end function

       real*8 function n3lo_vls_ct(z)
       real*8 z
       n3lo_vls_ct=c(18)*normq(z)**2+c(19)*normk(z)**2 
       return
      end function

       real*8 function n3lo_vt_ct(z)
       real*8 z
       n3lo_vt_ct=c(20)*normq(z)**2+c(21)*normk(z)**2
       return
       end function

       real*8 function n3lo_vsk_ct(z)
       real*8 z
       n3lo_vsk_ct=c(22)*normq(z)**2+c(23)*normk(z)**2
       return
      end function

       real*8 function n3lo_vslsl_ct(z)
       real*8 z
       n3lo_vslsl_ct=c(24)
       return
       end function


c    pi-gamma (charge dependent) (wt)
        real*8 function pi_gamma(z)
        real*8 z
        real*8 beta
        beta=normq(z)/mpi
        pi_gamma=-ga**2/(4.0d0*fpi**2*mpi**2*pi*alpha)
     +  *(-1.0d0/4.0d0*(1.0d0-beta**2)**2/(2.0d0*beta**4
     +  *(1.0d0+beta**2))*dlog(1.0d0+beta**2)+1.0d0/
     +  (2.0d0*beta**2))
        return
      end function
      end module

      


c   we write another legendre 
        function legendre(x,j)
c
c
c   x is the independent variable
c
        real*8  x,legendrem1,a,b
        real*8  legendre
        integer j
c
c
c
c        compute legendre polynom for j equals zero
c
c
        if (j.gt.0) go to 1
        legendre=1.d0
        legendrem1=0.d0
        if (j.lt.0) legendre=0.d0
        return
c
c
c
c        compute legendre polynoms for j equals one
c
c
c
    1 legendre=x
        legendrem1=1.d0
        if (j.eq.1) return
c
c
c
c        compute legendre polynom for j greater or equal two
c
c
c
        do i=2,j
        a=x*legendre
        b=a-legendrem1
        legendrem1=legendre
        legendre=-b/dfloat(i)+b+a
        end do
c

        return
        end

        function alj(func,l,j)
c   function alj calculate the integral a^j(l), l is obit-angular
c   momentum and j is total angular momentum,func is the function 
c   being integrated

c   ALJ^J(q',q)=pi*INTEGRATE[fuc(q',q,z) z^l P_J(z)，-1,1]
        implicit none
c    input        
        real*8,external :: func
        real*8,external :: legendre
        integer :: l,j
c    output        
        real*8 :: alj
c    local
        real*8 :: pi   
        real*8 ::ct(96),wt(96)  
        integer ::i  

        alj=0.0d0
        pi=3.141592653589793d0
        call fset(ct,wt,-1.0d0,1.0d0,20)
        do i=1,20
        alj=alj+wt(i)*(legendre(ct(i),j)*func(ct(i))*ct(i)**l*pi)
        end do
        return
        end

c   this function calculate the cutoff
c   lambda is cutoff energy ,n adjust the sharp degree 

        function cutoff(lambda,n)
         use potential_global,only:xmev,ymev
        implicit real*8 (a-h,o-z)
        real*8 lambda,t,expo
        real*8 cutoff
        integer n
        t=dfloat(n)
        expo=(xmev/lambda)**(2.0d0*t)+(ymev/lambda)**(2.0d0*t)
        cutoff=dexp(-expo)  
        end
        module decompose
c   this module calculate lsj decomposition
c   we use pot(6) to record them
c   0v(singlet), 1v(uncoupled triplet), v++, v--, v+-, v-+ (coupled) 
         implicit none
c   all the subroutines are public                 

         contains
      subroutine lsjvcentral(pot,vcentral,j)

c        input         
         real*8,external :: vcentral
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         
         pot=0.0d0
         if (j .eq. 0) then
         pot(1)=2.0d0*alj(vcentral,0,j)
         pot(3)=2.0d0*alj(vcentral,0,j+1)
         else 
         pot(1)=2.0d0*alj(vcentral,0,j)
         pot(2)=2.0d0*alj(vcentral,0,j)
         pot(3)=2.0d0*alj(vcentral,0,j+1)
         pot(4)=2.0d0*alj(vcentral,0,j-1)
         pot(5)=0.0d0
         pot(6)=0.0d0
         end if
        end   

        subroutine lsjvspinspin(pot,vspinspin,j)

c        input         
         real*8,external :: vspinspin
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj

         pot=0.0d0        
         if (j .eq. 0)then
            pot(1)=-6.0d0*alj(vspinspin,0,j)
            pot(3)=2.0d0*alj(vspinspin,0,j+1)
            else 
            pot(1)=-6.0d0*alj(vspinspin,0,j)
            pot(2)=2.0d0*alj(vspinspin,0,j)
            pot(3)=2.0d0*alj(vspinspin,0,j+1)
            pot(4)=2.0d0*alj(vspinspin,0,j-1)
            pot(5)=0.0d0
            pot(6)=0.0d0
            end if   
        end   
       
      subroutine lsjvtensor(pot,vtensor,j)

c        global(x,y)
         use paravari,only:x,y,x2,y2
                
c        input         
         real*8,external :: vtensor
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
         if (j .eq. 0) then
        pot(1)=2.0d0*(-(x2+y2)*alj(vtensor,0,j)
     1  +2.0d0*x*y*alj(vtensor,1,j))
        pot(3)=2.0d0/jd2p1*(-(x2+y2)*alj(vtensor,0,j+1)
     1  +2.0d0*x*y*alj(vtensor,0,j))
        else
        pot(1)=2.0d0*(-(x2+y2)*alj(vtensor,0,j)
     1  +2.0d0*x*y*alj(vtensor,1,j))
        pot(2)=2.0d0*((x2+y2)*alj(vtensor,0,j)
     1  -2.0d0*x*y/jd2p1*(jd*alj(vtensor,0,j+1)
     2  +jdp1*alj(vtensor,0,j-1)))
        pot(3)=2.0d0/jd2p1*(-(x2+y2)*alj(vtensor,0,j+1)
     1  +2.0d0*x*y*alj(vtensor,0,j))
        pot(4)=2.0d0/jd2p1*((x2+y2)*alj(vtensor,0,j-1)
     1  -2.0d0*x*y*alj(vtensor,0,j))
        pot(5)=-4.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vtensor,0,j+1)
     1  +x2*alj(vtensor,0,j-1)-2.0d0*x*y*alj(vtensor,0,j) )   
        pot(6)=-4.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vtensor,0,j-1)
     1  +x2*alj(vtensor,0,j+1)-2.0d0*x*y*alj(vtensor,0,j) )
        end if 
        end 

        subroutine lsjvspinobit(pot,vspinobit,j)
c        global(x,y)
         use paravari,only:x,y
                
c        input         
         real*8,external :: vspinobit
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
        if (j .eq. 0) then
        pot(1)=pot(1)+0.0d0
        pot(3)=pot(3)+2.0d0*x*y*(jd+2.0d0)/(2.0d0*jd+3.0d0)
     1  *(alj(vspinobit,0,j+2)-alj(vspinobit,0,j))
        else
        pot(1)=pot(1)+0.0d0
        pot(2)=pot(2)+2.0d0*x*y/jd2p1*(alj(vspinobit,0,j+1)
     1  -alj(vspinobit,0,j-1))
        pot(3)=pot(3)+2.0d0*x*y*(jd+2.0d0)/(2.0d0*jd+3.0d0)
     1  *(alj(vspinobit,0,j+2)-alj(vspinobit,0,j))
        pot(4)=pot(4)+2.0d0*x*y*(jd-1.0d0)/(2.0d0*jd-1.0d0)
     1  *(alj(vspinobit,0,j-2)-alj(vspinobit,0,j))
        pot(5)=pot(5)+0.0d0
        pot(6)=pot(6)+0.0d0
        end if
       end subroutine

       subroutine lsjvsigmal(pot,vsigmaL,j)
c        global(x,y,x2,y2)
         use paravari,only:x2,y2
                
c        input         
         real*8,external :: vsigmaL
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
         if (j .eq. 0) then
         pot(1)=pot(1)+2.0d0*x2*y2*(alj(vsigmaL,2,j)-alj(vsigmaL,0,j))
         pot(3)=pot(3)+2.0d0*x2*y2*((2.0d0*jd+3.0d0)/jd2p1
     1  *alj(vsigmaL,0,j+1)-2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j+1))
         else
         pot(1)=pot(1)+2.0d0*x2*y2*(alj(vsigmaL,2,j)-alj(vsigmaL,0,j))
         pot(2)=pot(2)+2.0d0*x2*y2*(-alj(vsigmaL,0,j)+(jd-1.0d0)/jd2p1
     1  *alj(vsigmaL,1,j+1)+(jd+2.0d0)/jd2p1*alj(vsigmaL,1,j-1))
         pot(3)=pot(3)+2.0d0*x2*y2*((2.0d0*jd+3.0d0)/jd2p1
     1  *alj(vsigmaL,0,j+1)-2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j+1))
         pot(4)=pot(4)+2.0d0*x2*y2*((2*jd-1.0d0)/jd2p1
     1  *alj(vsigmaL,0,j-1)+2.0d0/jd2p1*alj(vsigmaL,1,j)
     2  -alj(vsigmaL,2,j-1))
         pot(5)=pot(5)-4.0d0*x2*y2*dsqrt(jd*jdp1)/(jd2p1)**2*
     1  (alj(vsigmaL,0,j+1)-alj(vsigmaL,0,j-1))
         pot(6)=pot(5)
         end if 
         return
         end subroutine

         subroutine lsjvsigmak(pot,vsigmak,j)
c        global(x,y,x2,y2)
         use paravari,only:x,y,x2,y2
                
c        input         
         real*8,external :: vsigmak
         integer j
c        output
         real*8 pot(6)
c        local
         real*8,external :: alj
         real*8 jd,jdp1,jd2p1
         
         pot=0.0d0
         jd=dfloat(j)
         jdp1=jd+1.0d0
         jd2p1=2.0d0*jd+1.0d0
         if (j .eq. 0) then
        pot(1)=pot(1)+0.5d0*(-(x2+y2)*alj(vsigmak,0,j)
     1  -2.0d0*x*y*alj(vsigmak,1,j))
        pot(3)=pot(3)+0.5d0/jd2p1*(-(x2+y2)*alj(vsigmak,0,j+1)
     1  -2.0d0*x*y*alj(vsigmak,0,j))
        else
        pot(1)=pot(1)+0.5d0*(-(x2+y2)*alj(vsigmak,0,j)
     1  -2.0d0*x*y*alj(vsigmak,1,j))
        pot(2)=pot(2)+0.5d0*((x2+y2)*alj(vsigmak,0,j)
     1  +2.0d0*x*y/jd2p1*(jd*alj(vsigmak,0,j+1)
     2  +jdp1*alj(vsigmak,0,j-1)))
        pot(3)=pot(3)+0.5d0/jd2p1*(-(x2+y2)*alj(vsigmak,0,j+1)
     1  -2.0d0*x*y*alj(vsigmak,0,j))
        pot(4)=pot(4)+0.5d0/jd2p1*((x2+y2)*alj(vsigmak,0,j-1)
     1  +2.0d0*x*y*alj(vsigmak,0,j))
        pot(5)=pot(5)-1.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vsigmak,0,j+1)
     1  +x2*alj(vsigmak,0,j-1)+2.0d0*x*y*alj(vsigmak,0,j) )   
        pot(6)=pot(6)-1.0d0*dsqrt(jd*jdp1)/jd2p1*(y2*alj(vsigmak,0,j-1)
     1  +x2*alj(vsigmak,0,j+1)+2.0d0*x*y*alj(vsigmak,0,j) )
        end if
        return
       end subroutine
    
        subroutine isospindependent(pot1,j,pot2)
c       this subroutine contains the tau1 dot tau2 terms factor

c       input 
         real*8 pot1(6)
         integer j

c       output
         real*8 pot2(6)
         pot2=0.0d0
         if (mod(j,2).eq.0)then
            pot2(1)=pot1(1)
            pot2(2)=-3.0d0*pot1(2)
            pot2(3:6)=pot1(3:6)
         else
            pot2(1)=-3.0d0*pot1(1)
            pot2(2)=pot1(2)
            pot2(3:6)=-3.0d0*pot1(3:6)
         end if                
         end
      end module
        subroutine lsjdecomposition(pot,j)
        use addterms
        use potential_global,only:lambda
        use n3lopot
        use lopot
        use decompose
        implicit none
        real*8 pot(6),jd,jdp1,jd2p1,temp1(6),temp2(6),vcct(6)
        real*8 nlowcel(6),nlovtel(6)
        real*8 n2lovcel(6),n2lowsel(6),n2lowtel(6)
        real*8,external :: cutoff
        real*8,external :: alj
        real*8 ::ex,ey,ree,expexp2,expexp3,expexp4
        integer ::j,i
        call ini_paravari
        pot=0.0d0
        jd=dfloat(j)
        jdp1=jd+1.0d0
        jd2p1=2.0d0*jd+1.0d0

c   central force part
c   when j=0,there is only two terms do not equal to zero
            
        ex=dsqrt(1.0d0+x*x)
        ey=dsqrt(1.0d0+y*y)
        ree=dsqrt(ex*ey)
        expexp2=cutoff(lambda,2)
        expexp3=cutoff(lambda,3)
        expexp4=cutoff(lambda,4)
c       lo 

c       contact terms
        call lo_ct%init
        call lo_onepi%init
        call lsjvcentral(lo_ct%vc,lo_vc_ct,j)
        call lsjvspinspin(lo_ct%vss,lo_vs_ct,j)
        call lo_ct%add
        lo_ct%sum=lo_ct%sum*expexp3

c       one-pion terms
        call lsjvtensor(temp1,onepii0,j)
        call lsjvtensor(temp2,onepii1,j)
        if(mod(j,2) .eq. 1)then 
         lo_onepi%vt(1)=temp1(1)
         lo_onepi%vt(2)=temp2(2)
         lo_onepi%vt(3:6)=temp1(3:6)
         else
         lo_onepi%vt(1)=temp2(1)
         lo_onepi%vt(2)=temp1(2)
         lo_onepi%vt(3:6)=temp2(3:6)
         end if
         call lo_onepi%add
         lo_onepi%sum=lo_onepi%sum*expexp4
c        nlo

c        contact terms
         call nlo_ct%init
         call lsjvcentral(nlo_ct%vc,nlo_vc_ct,j)
         call lsjvspinspin(nlo_ct%vss,nlo_vs_ct,j)
         call lsjvspinobit(nlo_ct%vls,nlo_vls_ct,j)
         call lsjvtensor(nlo_ct%vt,nlo_vt_ct,j)
         call lsjvsigmak(nlo_ct%vsk,nlo_vsk_ct,j)
         call nlo_ct%add
c        the C_{3P1} term n=3         
         if(j .eq. 1)then
            nlo_ct%sum(2)=nlo_ct%sum(2)*expexp3
            nlo_ct%sum(1)=nlo_ct%sum(1)*expexp2
            nlo_ct%sum(3:6)=nlo_ct%sum(3:6)*expexp2
         else
            nlo_ct%sum=nlo_ct%sum*expexp2
         end if

c       2-pi terms
        call nlo_tp%init 
        call lsjvcentral(temp1,nlowc,j)
        call isospindependent(temp1,j,nlo_tp%wc)
        call lsjvspinspin(nlo_tp%vss,nlovss,j)
        call lsjvtensor(nlo_tp%vt,nlovt,j)
        call nlo_tp%add
        nlo_tp%sum=nlo_tp%sum*expexp2


        call lsjvcentral(n2lovcel,n2lovc,j)
        pot=pot+n2lovcel

c   spin-spin force part
       call lsjvspinspin(temp1,vspinspin,j)
c       pot=pot+temp1

       call lsjvspinspin(temp1,n2lows,j)
       call isospindependent(temp1,j,n2lowsel)
       pot=pot+n2lowsel   

c   spin-obit force part
       call lsjvspinobit(temp1,vspinobit,j)
c       pot=pot+temp1
c   sigmaL force part
        
       call lsjvsigmal(temp1,vsigmaL,j)
c       pot=pot+temp1
        
c   tensor force part
        call lsjvtensor(temp1,vtensor,j)
c        pot=pot+temp1



        call lsjvtensor(temp1,n2lowt,j)
        call isospindependent(temp1,j,n2lowtel)
        pot=pot+n2lowtel
c   sigmak force part

        call lsjvsigmak(temp1,vsigmak,j)
c        pot=pot+temp1
        call lsjvtensor(temp1,pi_gamma,j)
        call isospindependent(temp1,j,temp2)
c        pot=pot+temp2
        pot=pot*expexp2
c   n3lo
        
c      contact terms
        call n3lo_ct%init
        call lsjvcentral(n3lo_ct%vc,n3lo_vc_ct,j)
        call lsjvspinspin(n3lo_ct%vss,n3lo_vs_ct,j)
        call lsjvspinobit(n3lo_ct%vls,n3lo_vls_ct,j)
        call lsjvtensor(n3lo_ct%vt,n3lo_vt_ct,j)
        call lsjvsigmak(n3lo_ct%vsk,n3lo_vsk_ct,j)
        call lsjvsigmal(n3lo_ct%vslsl,n3lo_vslsl_ct,j)
        call n3lo_ct%add
        if(j .eq.0)then
c       D_{3P0}         
         n3lo_ct%sum(3)=n3lo_ct%sum(3)*expexp3
         n3lo_ct%sum(1:2)=n3lo_ct%sum(1:2)*expexp2
         n3lo_ct%sum(4:6)=n3lo_ct%sum(4:6)*expexp2
        else if(j.eq.1)then
c       D_{3P1}         
         n3lo_ct%sum(2)=n3lo_ct%sum(2)*expexp3
         n3lo_ct%sum(1)=n3lo_ct%sum(1)*expexp2
         n3lo_ct%sum(3:6)=n3lo_ct%sum(3:6)*expexp2
        else if(j.eq.2)then
c       D_{1D2}         
         n3lo_ct%sum(1)=n3lo_ct%sum(1)*expexp3
c       D_{3D2}         
         n3lo_ct%sum(2)=n3lo_ct%sum(2)*expexp3
c       D_{3PF2}         
         n3lo_ct%sum(5:6)=n3lo_ct%sum(5:6)*expexp4
         n3lo_ct%sum(3:4)=n3lo_ct%sum(3:4)*expexp2
        else
         n3lo_ct%sum=n3lo_ct%sum*expexp2
        end if



c      pi exchange terms
        call n3lo_rc%init
        call n3lo_fd%init
        call n3lo_tl%init
        call n3lo_cM%init
        call lsjvcentral(n3lo_rc%vc,n3lo_vc_rc,j)
        call lsjvtensor(n3lo_rc%vt,n3lo_vt_rc,j)
        call lsjvspinspin(n3lo_rc%vss,n3lo_vs_rc,j)
        call lsjvspinobit(n3lo_rc%vls,n3lo_vls_rc,j)
        call lsjvcentral(temp1,n3lo_wc_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wc)
        call lsjvtensor(temp1,n3lo_wt_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wt)
        call lsjvspinspin(temp1,n3lo_ws_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wss)
        call lsjvspinobit(temp1,n3lo_wls_rc,j)
        call isospindependent(temp1,j,n3lo_rc%wls)
        call n3lo_rc%add
        n3lo_rc%sum=n3lo_rc%sum*expexp2
        call lsjvcentral(n3lo_fd%vc,n3lo_vc_fd,j)
        call lsjvtensor(temp1,n3lo_wt_fd,j)
        call isospindependent(temp1,j,n3lo_fd%wt)
        call lsjvspinspin(temp1,n3lo_ws_fd,j)
        call isospindependent(temp1,j,n3lo_fd%wss)
        call n3lo_fd%add
        n3lo_fd%sum=n3lo_fd%sum*expexp2
        call lsjvcentral(n3lo_tl%vc,n3lo_vc_tl,j)
        call lsjvspinspin(temp1,n3lo_ws_tl,j)
        call isospindependent(temp1,j,n3lo_tl%wss)
        call lsjvtensor(temp1,n3lo_wt_tl,j)
        call isospindependent(temp1,j,n3lo_tl%wt)
        call lsjvspinspin(n3lo_tl%vss,n3lo_vs_tl,j)
        call lsjvtensor(n3lo_tl%vt,n3lo_vt_tl,j)
        call lsjvcentral(temp1,n3lo_wc_tl,j)
        call isospindependent(temp1,j,n3lo_tl%wc)
        call n3lo_tl%add 
        n3lo_tl%sum=n3lo_tl%sum*expexp2
        call lsjvcentral(n3lo_cM%vc,n3lo_vc_cM,j)
        call lsjvcentral(temp1,n3lo_wc_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wc)
        call lsjvtensor(temp1,n3lo_wt_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wt)
        call lsjvspinspin(temp1,n3lo_ws_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wss)
        call lsjvspinobit(n3lo_cM%vls,n3lo_vls_cM,j)
        call lsjvspinobit(temp1,n3lo_wls_cM,j)
        call isospindependent(temp1,j,n3lo_cM%wls)
        call n3lo_cM%add
        n3lo_cM%sum=n3lo_cM%sum*expexp2
         
c     sum all
        pot=lo_ct%sum+lo_onepi%sum+nlo_ct%sum+nlo_tp%sum+pot
     + +n3lo_tl%sum+n3lo_ct%sum+n3lo_fd%sum+n3lo_cM%sum+n3lo_rc%sum  
        do i=1,6
        pot(i)=pot(i)/(2.0d0*pi)**3*dwnq/ree
        end do
        
        return
        end

c   this function write the 3P2 contact terms 
         subroutine n3lo500new
         use potential_global,only:v,j
         call lsjdecomposition(v,j)
         return
         end subroutine

      subroutine fset(ct,wt,xlb,xub,n)
c     the intergral is sun( wt(i)*f(ct(i))
c       the integrate low bound and up bound        
        real*8 ::xlb,xub
c     the gauss point number        
        integer ::n
c       output  
        real*8 ::ct(96),wt(96)        
c     local
        real*8 ::x(273),a(273)
        integer :: i
c     n=8
      data x(16)/0.960289856497536 d0/, a(16)/0.101228536290376 d0/
      data x(17)/0.796666477413627 d0/, a(17)/0.222381034453374 d0/
      data x(18)/0.525532409916329 d0/, a(18)/0.313706645877887 d0/
      data x(19)/0.183434642495650 d0/, a(19)/0.362683783378362 d0/
c     n=12
      data x(36)/0.981560634246719 d0/, a(36)/0.047175336386512 d0/
      data x(37)/0.904117256370475 d0/, a(37)/0.106939325995318 d0/
      data x(38)/0.769902674194305 d0/, a(38)/0.160078328543346 d0/
      data x(39)/0.587317954286617 d0/, a(39)/0.203167426723066 d0/
      data x(40)/0.367831498998180 d0/, a(40)/0.233492536538355 d0/
      data x(41)/0.125233408511469 d0/, a(41)/0.249147045813403 d0/
c     n=16
      data x(64)/0.989400934991650 d0/, a(64)/0.027152459411754 d0/
      data x(65)/0.944575023073233 d0/, a(65)/0.062253523938648 d0/
      data x(66)/0.865631202387832 d0/, a(66)/0.095158511682493 d0/
      data x(67)/0.755404408355003 d0/, a(67)/0.124628971255534 d0/
      data x(68)/0.617876244402644 d0/, a(68)/0.149595988816577 d0/
      data x(69)/0.458016777657227 d0/, a(69)/0.169156519395003 d0/
      data x(70)/0.281603550779259 d0/, a(70)/0.182603415044924 d0/
      data x(71)/0.095012509837637 d0/, a(71)/0.189450610455069 d0/      
c      n=20
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
c       n=24       
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
c       n=96        
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
        wt=0.0d0 
        ct=0.0d0
        select case(n)
        case(8)
         do  i=1,n/2
         wt(i)=a(15+i)   
         ct(i)=-x(15+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do  
      case(12)
         do  i=1,n/2
         wt(i)=a(35+i)   
         ct(i)=-x(35+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do  
      case(16)
         do  i=1,n/2
         wt(i)=a(63+i)   
         ct(i)=-x(63+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      case(20)
        do  i=1,n/2
        wt(i)=a(71+i)   
        ct(i)=-x(71+i)
        end do
        do  i=n,n/2,-1
        ct(i)=-ct(n+1-i)
        wt(i)=wt(n+1-i)
        end do
      case(24)
         do  i=1,n/2
         wt(i)=a(81+i)   
         ct(i)=-x(81+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      case(64)
         do  i=1,n/2
         wt(i)=a(153+i)   
         ct(i)=-x(153+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      case(96)
         do  i=1,n/2
         wt(i)=a(225+i)   
         ct(i)=-x(225+i)
         end do
         do  i=n,n/2,-1
         ct(i)=-ct(n+1-i)
         wt(i)=wt(n+1-i)
         end do
      end select
        do i=1,96
         ct(i)=xlb+(ct(i)+1.0d0)*(xub-xlb)/2.0d0
        end do
        wt=wt*(xub-xlb)/2.0d0
        return             
      end subroutine