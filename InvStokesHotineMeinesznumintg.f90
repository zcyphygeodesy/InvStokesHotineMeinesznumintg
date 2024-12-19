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
      !输入扰动重力逆运算积分类型knd
      !knd=0 Inverse Stokes, =1 Inverse Hotine =2 Inverse Vening-Meinesz
      knd=1!knd=0-逆Stokes,1-逆Hotine,2-逆Vening-Meinesz
      !输入计算点文件名，程序默认计算点在等位边界面上
      !Input the calculation point file on the equipotential boundary surface.
      write(calcpntfl,*)'calcpnt.txt'
      !输入等位边界面大地高格网文件名
      !Input the ellipsoidal height grid file of the equipotential surface.
      write(dwmhgrdfl,*)'landgeoidhgt.dat'
      !输入等位面上残差场元格网文件名。knd=0或1时输入残差高程异常格网，knd=2时输入残差垂线偏差向量格网。
      !Input residual height anomaly (m, when knd=0 or 1) orvertical deflection vector
      !(″, SW, when knd=2) grid file on the equipotential surface.
      !程序要求等位面大地高格网及其面上残差场元格网有相同的格网规格
      !The same grid specifications required for the ellipsoidal height grid of the equipotential
      !surface and residual field element grid on the surface.
      write(gravgrdfl,*)'resGMlgeoid541_1800.ksi'
      !输入积分半径(m)
      dr=150.d3!Integral radius (m)
      write(*, *)"    Begin compulation......"
      call InversenumIntegral(calcpntfl,dwmhgrdfl,gravgrdfl,knd,dr)
      pause
      end
