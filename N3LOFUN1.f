      module phasecal
        real*8,save,allocatable :: phaseshifts(:,:)
        integer,save :: phasenum
        integer,save,allocatable :: phasename(:,:)
      end module
      program main
             use phasecal
             implicit real*8 (a-h,o-z)
            
c     we use common conta(2)to pass on this para
           

c     we use phaseshifts to pass out the useful phaseshifts
c     to be specific,it's  
             real*8 elab(41)
             parameter (n=64)
             parameter (n6=6*n,na=2*(n+1)*2*(n+1),naa=6*n/2*(n+1))
             dimension vv(n6),s(n),u(n),a(na),b(na),aa(naa),qq(n),eq(n)
             character*2 filenum
             external n3lo500new
             common /alpha/ melab
             common /einject/ elab
             common /crdwrt/ kread,kwrite,kpunch,kda(9)          
             call getarg(1,filenum)
            open (unit=10,file='input'//filenum)
            open (unit=11,file='output'//filenum)
            kread=10
            kwrite=11
            read(kread,*) phasenum
            allocate(phasename(phasenum,2))
            do i=1,phasenum
               read(kread,*) phasename(i,:)
            end do
            call phases (n3lo500new,vv,s,u,a,b,aa,qq,eq)
10055 format (f6.1,20f24.20)
            do i=1,melab
               write(kwrite,10055) elab(i),phaseshifts(i,:)
            end do
            deallocate(phasename)
            deallocate(phaseshifts)
            close(11)
            close(10)
            end
