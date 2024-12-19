      subroutine HGradientOne(BLH0,rga,hgt,nlat,nlon,hd,tx,GRS)
      !BLH点rga的1阶水平梯度估计。x轴指向北,y轴指向东
      implicit none
	integer nlon,nlat,i,j,ki,kj,si,sj
	real*8::hd(6),rga(nlat,nlon),hgt(nlat,nlon),CGrdPntD2,tx(5)
	real*8::RAD,pi,rr,rlat,lon,rlat0,lon0,dlat,dlon,dx,dy,obs,dl,tmp
	real*8::GRS(6),BLH(3),BLH0(3),rln(3),BPB(2,2),BB(2),BPL(2),xx(2)
!---------------------------------------------------------------------
   	pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      BLH0(3)=CGrdPntD2(BLH0(2),BLH0(1),hgt,nlat,nlon,hd)
      call BLH_RLAT(GRS,BLH0,rln)
	rlat0=rln(2)*RAD;lon0=rln(3)*RAD
	rr=rln(1);BPB=0.d0;BPL=0.d0;xx=0.d0
	i=nint((BLH0(1)-hd(3))/hd(6)+0.5d0)
	j=nint((BLH0(2)-hd(1))/hd(5)+0.5d0)
	do ki=i-1,i+1
	   BLH(1)=hd(3)+(real(ki)-0.5d0)*hd(6)
	   do kj=j-1,j+1
	      BLH(2)=hd(1)+(real(kj)-0.5d0)*hd(5)
	      if(ki<1.or.ki>nlat.or.kj<1.or.kj>nlon)goto 9200
	      if(ki==i.and.kj==j)goto 9200
		  BLH(3)=hgt(ki,kj)
            call BLH_RLAT(GRS,BLH,rln)
	      rlat=rln(2)*RAD;lon=rln(3)*RAD
	      dlon=lon-lon0;dlat=rlat-rlat0
		  dx=2.d0*dsin(dlat/2.d0)*rr
		  dy=2.d0*dsin(dlon/2.d0)*rr*dcos((rlat+rlat0)/2.d0)
	      dl=dsqrt(dx**2+dy**2)
            if(dl<1.d-3)goto 9200
            obs=rga(ki,kj)-rga(i,j)
            BB(1)=dx;BB(2)=dy
	      do si=1,2
	         BPL(si)=BPL(si)+BB(si)*obs/dl
	         do sj=1,2
	            BPB(si,sj)=BPB(si,sj)+BB(si)*BB(sj)/dl
	         enddo
	       enddo
9200	       continue
	    enddo
	enddo
	tmp=0.d0
	do si=1,2
	   do sj=1,si-1
	      BPB(sj,si)=BPB(si,sj)
	   enddo
         tmp=tmp+BPB(si,si)**2
	enddo
	tmp=dsqrt(tmp/2.d0)
	do si=1,2
	   BPB(si,si)=BPB(si,si)+tmp*1.d-5
	enddo
	call EqJordon(BPB,xx,2,BPL)
      tx(1)=xx(1);tx(2)=xx(2)
      end
