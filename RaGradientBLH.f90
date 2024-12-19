      real*8 function RaGradientBLH(BLH,ksi,hgt,nlat,nlon,hd,dr,GRS)
      !按严密球面积分公式计算重力场元的径向梯度
      !dr-积分半径m
!-------------------------------------------------------------
      implicit none
	integer::i,j,nlat,nlon,i0,j0,ni,nj,astat(5)
	real*8::dr,ksi(nlat,nlon),hgt(nlat,nlon)
	real*8::hd(6),gp,pi,RAD,ds,mdr,tt,rr,r1,rst,tmp
	real*8::GRS(6),BLH(3),XYZ(3),rln(3),BLH1(3),XYZ1(3),rln1(3)
	real*8 CGrdPntD2,RaGradSng,ksi0,L1
!-----------------------------------------------------------------
      rst=0.d0;pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      ksi0=CGrdPntD2(BLH(2),BLH(1),ksi,nlat,nlon,hd)
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      mdr=rr*hd(5)*RAD*dcos(rln(2)*RAD)/4.d0 !奇异点判断
      ni=nint(dr/rr/RAD/hd(6)+1.d0) !积分半径dr对应的地面格网数
      nj=nint(dr/rr/RAD/hd(5)/dcos(rln(2)*RAD)+1.d0)
	i0=nint((BLH(1)-hd(3))/hd(6)+0.5d0)
	j0=nint((BLH(2)-hd(1))/hd(5)+0.5d0)!计算点所在的地面格网i0,j0
	do i=i0-ni,i0+ni
	  if(i<1.or.i>nlat)goto 9100
        BLH1(1)=hd(3)+(real(i)-0.5d0)*hd(6)
	  do j=j0-nj,j0+nj
	    if(j<1.or.j>nlon)goto 9101
	    BLH1(2)=hd(1)+(real(j)-0.5d0)*hd(5)
          BLH1(3)=hgt(i,j);call BLH_XYZ(GRS,BLH1,XYZ1)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          if(L1>dr)goto 9101
          if(L1<mdr)goto 9101!!!!!!!!!!!!!!!!!!!!
          if(L1<mdr)then!计算奇异积分
             tmp=RaGradSng(BLH,ksi,hgt,nlat,nlon,hd,i,j,2,GRS)
             rst=rst+tmp
             goto 9101 
          endif
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)*r1**2
          rst=rst+(ksi(i,j)-ksi0)*ds/pi/L1**3/2.d0
9101      continue
	  enddo
9100    continue
	enddo
	RaGradientBLH=rst
9002	return
      end
!********************************************************************************
      real*8 function RaGradSng(BLH,ksi,hgt,nlat,nlon,hd,i0,j0,m,GRS)
      !计算BLH点的重力场元的径向梯度的奇异积分
!-------------------------------------------------------------
      implicit none
	integer::nlat,nlon,i0,j0,m,i,j
	real*8::CGrdPntD2,ksi(nlat,nlon),hgt(nlat,nlon),ksi1,rv,mdr,rst,ksi0,L1
	real*8::hd(6),pi,RAD,ds,r0,r1,GRS(6),BLH(3),BLH1(3),XYZ(3),XYZ1(3),rln1(3),lat,lon
!-----------------------------------------------------------------
      rv=hd(5)/dble(m);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      lat=hd(3)+real(i0-1)*hd(6);lon=hd(1)+real(j0-1)*hd(5)!格网左下角经纬度
      ksi0=CGrdPntD2(BLH(2),BLH(1),ksi,nlat,nlon,hd)
      call BLH_XYZ(GRS,BLH,XYZ);r0=dsqrt(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
      rst=0.d0;mdr=r0*rv*RAD*dcos(BLH(1)*RAD)/dble(m)/4.d0  !奇异点判断
	do i=1,m
        BLH1(1)=lat+(real(i)-0.5d0)*rv
	  do j=1,m
	    BLH1(2)=lon+(real(j)-0.5d0)*rv
          BLH1(3)=CGrdPntD2(BLH1(2),BLH1(1),hgt,nlat,nlon,hd)
          ksi1=CGrdPntD2(BLH1(2),BLH1(1),ksi,nlat,nlon,hd)
          call BLH_XYZ(GRS,BLH1,XYZ1)
          L1=dsqrt((XYZ1(1)-XYZ(1))**2+(XYZ1(2)-XYZ(2))**2+(XYZ1(3)-XYZ(3))**2)
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          ds=rv**2*RAD**2*dcos(rln1(2)*RAD)*r1**2
          if(L1<mdr)L1=mdr
          rst=rst+(ksi1-ksi0)*ds/pi/L1**3/2.d0
	  enddo
9100    continue
	enddo
      RaGradSng=rst
9002	return
      end
