      subroutine AnVMdeflBLH(BLH,kai,nta,hgt,nlat,nlon,hd,dr,rst,GRS)
      !按严密球面积分公式计算垂线偏差逆远算积分
      !dr-积分半径m
!-------------------------------------------------------------
      implicit none
	integer::knd,i,j,nlat,nlon,i0,j0,ni,nj,astat(5)
	real*8::dr,kai(nlat,nlon),nta(nlat,nlon),hgt(nlat,nlon)
	real*8::GRS(6),hd(6),pi,RAD,ds,mdr,gr,tt,rr,r1,rst(3),tmp,qp
	real*8::BLH(3),XYZ(3),rln(3),NFD(5),BLH1(3),XYZ1(3),rln1(3)
	real*8 rlon,rlat,rlon1,rlat1,dlon,dlat,sgn(3)
	real*8 sin2f,cos2f,sinf,sina,cosa,CGrdPntD2,L1
!-----------------------------------------------------------------
      rst=0.d0;pi=datan(1.d0)*4.d0;RAD=pi/180.d0;sgn=0.d0
      call BLH_XYZ(GRS,BLH,XYZ);call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
	call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      rlat=rln(2)*RAD;rlon=rln(3)*RAD
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
          call BLH_RLAT(GRS,BLH1,rln1);r1=rln1(1)
          ds=hd(5)*hd(6)*RAD**2*dcos(rln1(2)*RAD)
          if(L1<mdr)then!计算奇异积分
             !call AnVMSng(BLH,kai,nta,hgt,nlat,nlon,hd,sgn,GRS)
             !rst=rst+sgn; 
             goto 9101 
          endif
          tt=1.d0-2.d0*(L1/r1/2.d0)**2
          rlat1=rln1(2)*RAD;rlon1=rln1(3)*RAD
          sin2f=dsqrt((1.d0-tt)/2.d0);cos2f=dsqrt(1.d0-sin2f**2);sinf=2.d0*sin2f*cos2f
	    cosa=(dcos(rlat)*dsin(rlat1)-dsin(rlat)*dcos(rlat1)*dcos(rlon1-rlon))/sinf
	    sina=dcos(rlat1)*dsin(rlon1-rlon)/sinf
          qp=-(kai(i,j)*cosa+nta(i,j)*sina)*ds/pi/4.d0
          rst(1)=rst(1)-qp*r1*cos2f/sin2f
	    tmp=3.d0/sinf-1.d0/sinf/sin2f-sin2f/cos2f-2.d0*cos2f/sin2f
          rst(2)=rst(2)+qp*gr*tmp
          rst(3)=rst(3)+qp*gr*(3.d0/sinf-1.d0/sinf/sin2f-sin2f/cos2f)
9101      continue
	  enddo
9100    continue
	enddo
9002	return
      end
!--------------------------------------------------------------------------------
      subroutine AnVMSng(BLH,kai,nta,hgt,nlat,nlon,hd,rst,GRS)
      !计算BLH点的垂线偏差逆远算的奇异积分
!-------------------------------------------------------------
      implicit none
	integer::nlat,nlon
	real*8::kai(nlat,nlon),nta(nlat,nlon),hgt(nlat,nlon),kax(5),ntx(5)
	real*8::hd(6),pi,RAD,ds,rr,GRS(6),BLH(3),rln(3),NFD(5),gr,rst(3)
!-----------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      call HGradientOne(BLH,kai,hgt,nlat,nlon,hd,kax,GRS)
      call HGradientOne(BLH,nta,hgt,nlat,nlon,hd,ntx,GRS)
      call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
      call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      ds=hd(5)*hd(6)*RAD**2*dcos(rln(2)*RAD)*rr**2
      rst(1)=ds*(kax(2)+ntx(1))/4.d0/pi
      rst(2)=-dsqrt(ds/pi)*(kax(2)+ntx(1))/4.d0*gr
      rst(3)=-(dsqrt(ds*pi)+ds/rr)*(kax(2)+ntx(1))/2.d0/pi*gr
9002	return
      end
