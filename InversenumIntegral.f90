      subroutine InversenumIntegral(calcpntfl,dwmhgrdfl,gravgrdfl,knd,dr)
      implicit none
	character*800::calcpntfl,dwmhgrdfl,gravgrdfl
 	character*800::line,str,astr
      integer knd,sn,len,astat(8)
	integer i,j,nlon,nlat,kk
	real*8::dr,hd(6),hd1(6),hd2(6),rec(800),grad,ksi0
	real*8::GRS(6),BLH(3),XYZ(3),NFD(5),rr,gr
	real*8::CGrdPntD2,RaGradientBLH,rst(3),mr
	real*8,allocatable::kai(:,:),nta(:,:),dwm(:,:)
	integer::status=0
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378136.3d0; GRS(3)=1.082636277388d-3
      GRS(4) = 7.292115d-5; GRS(5) = 1.d0/298.25641153d0
      mr=datan(1.d0)/45.d0/3.6d3
      open(unit=8,file=gravgrdfl,status="old",iostat=status)
      if(status/=0)goto 902
      read(8,'(a)') line
      call PickRecord(line,len,rec,sn)
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/real(nlon)
	hd(6)=(hd(4)-hd(3))/real(nlat)
 	allocate(kai(nlat,nlon), stat=astat(1))
 	allocate(nta(nlat,nlon), stat=astat(2))
 	allocate(dwm(nlat,nlon), stat=astat(3))
	if (sum(astat(1:3)) /= 0) then
          close(8);goto 902
	endif
 	do i=1,nlat
	   read(8,*,end=903)(kai(i,j),j=1,nlon)
      enddo
      if(knd<2)goto 903
 	do i=1,nlat
	   read(8,*,end=903)(nta(i,j),j=1,nlon)
      enddo
      kai=kai*mr;nta=nta*mr
903   close(8)
      open(unit=10,file=dwmhgrdfl,status="old",iostat=status)
      if(status/=0)goto 904
      read(10,'(a)') line
      call PickRecord(line,len,rec,sn)
      hd1(1:6)=rec(1:6)
      if(sum(hd1-hd)>1.d-5)then !格网规格不同The grid specifications are different
         close(10);goto 904
      endif
 	do i=1,nlat
	   read(10,*,end=905)(dwm(i,j),j=1,nlon)
      enddo
905   close(10)
      open(unit=8,file=calcpntfl,status="old",iostat=status)
      if(status/=0)goto 904
      open(unit=10,file='reslt.txt',status="replace")
      read(8,'(a)') line  !读取头文件read the file header
      write(10,101)trim(line)
      kk=0
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickRecord(line,len,rec,sn)
         if(sn<4)goto 906
         kk=kk+1;BLH(2)=rec(2);BLH(1)=rec(3);rst=0.d0
         !内插等位面上计算点的大地高
         !interpolate the ellipsoidal height of the calculation point on the equipotential surface
         BLH(3)=CGrdPntD2(BLH(2),BLH(1),dwm,nlat,nlon,hd)
         if(knd<2)then
           call BLH_XYZ(GRS,BLH,XYZ);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
           rr=dsqrt(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
           !内插等位面上计算点的残差高程异常
           !interpolate the residual height anomaly at the calculation point on the equipotential surface
           ksi0=CGrdPntD2(BLH(2),BLH(1),kai,nlat,nlon,hd)
           !calculate the radial gradient of residual height anomaly
           grad=RaGradientBLH(BLH,kai,dwm,nlat,nlon,hd,dr,GRS)!计算高程异常径向梯度
           if(knd==0)rst(1)=-gr*grad
           if(knd==1)rst(1)=-gr*grad-ksi0*gr/2.d0/rr
           !输出计算点大地高和残差场元（mGal，knd=0空间异常，knd=1扰动重力）
           write(10,101)trim(line),BLH(3),rst(1)*1.d5
         else
           call AnVMdeflBLH(BLH,kai,nta,dwm,nlat,nlon,hd,dr,rst,GRS)
           !输出计算点大地高、残差高程异常(m)、残差空间异常(mGal)和残差扰动重力(mGal）
           write(10,101)trim(line),BLH(3),rst(1),rst(2:3)*1.d5
         endif
         if(kk/200*200==kk)write(*, '(a,i9)'), '    Calculated point number: ',kk
      enddo
906   close(8); close(10)
904   deallocate(kai,nta,dwm)
902   continue
101   format(a,40F12.4)
      write (*,*)'  Complete the computation! The results are saved in the file reslt.txt.'
      end
