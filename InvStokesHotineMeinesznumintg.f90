!  InvStokesHotineMeinesznumintg.f90 
!
!  FUNCTIONS:
!  InvStokesHotineMeinesznumintg - Entry point of console application.
!
!****************************************************************************
      program InvStokesHotineMeinesznumintg
      implicit none
	character*800::calcpntfl,dwmhgrdfl,gravgrdfl
	integer knd
	real*8::dr
!---------------------------------------------------------------------
      !�����Ŷ������������������knd
      !knd=0 Inverse Stokes, =1 Inverse Hotine =2 Inverse Vening-Meinesz
      knd=1!knd=0-��Stokes,1-��Hotine,2-��Vening-Meinesz
      !���������ļ���������Ĭ�ϼ�����ڵ�λ�߽�����
      !Input the calculation point file on the equipotential boundary surface.
      write(calcpntfl,*)'calcpnt.txt'
      !�����λ�߽����ظ߸����ļ���
      !Input the ellipsoidal height grid file of the equipotential surface.
      write(dwmhgrdfl,*)'landgeoidhgt.dat'
      !�����λ���ϲвԪ�����ļ�����knd=0��1ʱ����в�߳��쳣������knd=2ʱ����в��ƫ������������
      !Input residual height anomaly (m, when knd=0 or 1) orvertical deflection vector
      !(��, SW, when knd=2) grid file on the equipotential surface.
      !����Ҫ���λ���ظ߸����������ϲвԪ��������ͬ�ĸ������
      !The same grid specifications required for the ellipsoidal height grid of the equipotential
      !surface and residual field element grid on the surface.
      write(gravgrdfl,*)'resGMlgeoid541_1800.ksi'
      !������ְ뾶(m)
      dr=150.d3!Integral radius (m)
      write(*, *)"    Begin compulation......"
      call InversenumIntegral(calcpntfl,dwmhgrdfl,gravgrdfl,knd,dr)
      pause
      end
