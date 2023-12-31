! Author: Sarthak Gupta !
! Please contact sg207@rice.edu for any questions!


! This code analyze the time-dependent properties of the viral uptake !
! It takes the trajectory data generated by the Trig_March.f90 and !
! calculates the rate of wrapping from virus and cell surface perspective. !


module virus_uptake
implicit none
	integer :: TM,V,Layer,ACE2,CL,NE,NT,NSN,PosV,SVN,POS_ACE2,NVF,NN,NE_Comp,Lx,Comp_Spike,NE_Spike,SS,N_Spike,Vim_Edge
	Real*8  :: space
	REAL*8,allocatable,dimension(:)    :: RX,RY,RZ,RX0,RY0,RZ0
	REAL*8,allocatable,dimension(:)    :: Height
	Integer,allocatable,dimension(:)   :: E1,E2,Spike_occupy,N_Vim_Sp,RC,E_S1,E_S2
	integer :: istep,step,tstep,nf,NBP,Particle_No,INC
	integer :: VE
	integer :: Sp,Sp_Edge,Comp,Comp_Edge,Comp_Edge_complete
	Integer :: Comp_Edge_Sphere
	integer :: interval,cyclee,ll,oo
	integer :: Flag_connect,K_Vim,K_Spike,Flag_Spike_occupied
	real*8  :: dis,xr,yr,zr
	Real*8 :: xcm,ycm,zcm
	Real*8 :: xdis,ydis,zdis,Rg_sqr,SPOCP,r_vim_covid,Z_Avg
	Real*8 :: hx,hy,hz,h_dis,RMS_Rough_Mem,RMS_Rough_Virus
	Real*8 :: Wrap_Boundary,Wrap_distance,Wrap_distance_mem
	Real*8 :: Sigma_VIM
	integer :: num_args
	integer,allocatable,dimension(:) :: Flag_Bowl
	character(len=20), dimension(:), allocatable :: args

	Real*8  :: probe_dis

	Real*8 :: High_Diff,High_Diff_2,Hrs
	Integer :: Count_Inside

	Real*8  :: Pie

	Real*8 :: VX_COM,VY_COM,VZ_COM
	Real*8 :: Virus_Raidus,Check_Raidus
	Real*8 :: dis_check
	integer :: count_in
	Integer :: In_Circle



	Integer :: Ben,Sp_no,VD,Ens
	Character*20 :: Ben_ch,Sp_no_ch,VD_ch,Ens_ch
	character*50    :: filename1,filename2,filename3,filename4,filename5,filename14
	character*50    :: filename6,filename7,filename8,filename110,filename13
	character*50    :: filename15,filename16,filename17

	integer,allocatable,dimension(:,:) :: Neigh,Neigh_right
	integer,allocatable,dimension(:) :: Flag_6
	integer :: mem_point_count
	integer :: PC1,PC2
	integer,allocatable,dimension(:) :: Flag_in_area
	real*8 :: dyara
	real*8,allocatable,dimension(:) :: guass_sheet

	integer :: Total_Virus_Particle
	Real*8  :: Flat_Sheet_Area
	
End module virus_uptake




program endocytosis
use virus_uptake
implicit none
integer :: i
  istep = 1     	! Initial Step you want to start analysing data
  step  = 1	     	! Interval between two steps
  tstep = 2000	! Final step of simulation
  PIE=4.D0*DATAN(1.D0)    ! Just pie, you know.
  Sigma_VIM=1.d0 


num_args = command_argument_count()
allocate(args(num_args))  

do i = 1, num_args
  call get_command_argument(i,args(i))
end do

READ(args(1),*)Ben   !then, convert them to REALs
READ(args(2),*)Sp_no	
READ(args(3),*)VD
READ(args(4),*)Ens

Call Make_File_Names
Read(40,*)TM,V,Layer,ACE2,CL,NE,NT,space,NSN,PosV,SVN,POS_ACE2,NVF,NN,NE_Comp,Lx,Comp_Spike,NE_Spike,SS,N_Spike,Vim_Edge

open(unit=31,file="Sphere_Data.dat")
read(31,*) Sp,Sp_Edge
NBP=(sqrt(real(TM))*2) + (sqrt(real(TM))-2)*2
 CLOSE(31)

Comp=Sp+CL
Comp_Edge=NE+Vim_Edge
Comp_Edge_Sphere=Comp_Edge+Sp_Edge
Comp_Edge_complete=Comp_Edge_Sphere+NE_Spike
allocate(RX(Comp_Spike),RY(Comp_Spike),RZ(Comp_Spike))
allocate(RX0(Comp_Spike),RY0(Comp_Spike),RZ0(Comp_Spike))
allocate(Spike_occupy(N_Spike))
allocate(N_Vim_Sp(N_Spike))
allocate(RC(TM),Height(TM))
allocate(E1(Comp_Edge_complete),E2(Comp_Edge_complete))
allocate(E_S1(NSN),E_S2(NSN))
allocate(Flag_in_area(TM))

allocate(Neigh(TM,6),Neigh_right(TM,6))
allocate(Flag_6(TM))
allocate(guass_sheet(TM))
allocate(Flag_Bowl(TM))

Do I=1,Comp_Spike
	Read(10,*) RX0(I),RY0(I),RZ0(I)
End Do


Do I=1,Comp_Edge_complete
	Read(20,*) E1(I),E2(I)
End Do

Do I=1,NSN
	Read(60,*) E_S1(I),E_S2(I)
End Do

nf=tstep
dyara=1000.d0
call cut_membrane_ini

Total_Virus_Particle=Comp-Cl


!______ main loop __________!
	Do step=1,nf
		call Read_Data
			Call Bulge_Size
			Call Virus_Wrapping
			Call Velocity
	End do
!___________________________!

End Program endocytosis



Subroutine Velocity
use virus_uptake
implicit none
real*8 :: xr_ini,yr_ini,zr_ini
real*8 :: cos_angle,angle

	zr=rz0(Comp)-rz(Comp)
	write(533,*) step,Ben,VD,ZR
	
return
End Subroutine Velocity







! Candy Wrapper vs Taco Shape
Subroutine Bulge_Size
use virus_uptake
implicit none
integer :: i,j,k
real*8,dimension(TM)  :: H
integer,dimension(TM) :: Band_ID
real*8 :: Max_Rad,Band_Width
integer :: N_Band
Real*8 :: Total_Angle,Increment
integer :: N_Section
integer,dimension(TM) :: Section_ID
Real*8 :: Surface_Angle
Real*8 :: H_Cut_Off
integer :: Flag_Exit,Point_Count
Real*8 :: Membrane_Surface_Area
integer :: p1,p2,p3
real*8  :: Vx_p1p2,Vy_p1p2,Vz_p1p2
real*8  :: Vx_p1p3,Vy_p1p3,Vz_p1p3
real*8  :: Ax,Ay,Az,Area
Real*8 :: Whole_Surface_Area

	
!Find height of each particle for reference Z
open(unit=444,file="Height.dat")

	Do I=1,TM
		H(I)=RZ(I)-RZ(POS_ACE2)
	End Do
	
! Band Disect	
Max_Rad=3.0*(5.0)!sqrt(360.d0*360.d0 + 180.d0*180.d0)
Band_Width=0.1!1.5 ! Decides how far to look and make boundary tight
N_Band=int(Max_Rad/Band_Width)

! Angle Disect
Total_Angle=360.d0
Increment=0.1!1.d0		
N_Section=Total_Angle/Increment
	
Do I=1,TM
	Band_ID(I)=0
	Section_ID(I)=0
	Flag_Bowl(I)=0
End Do

! Giving each particle a band and section ID

Do I=1,TM	
	XR=RX(I)-RX(POS_ACE2)
	YR=RY(I)-RY(POS_ACE2)
	
	Dis=sqrt(XR*XR + YR*YR)
	
	Do J=1,N_Band
		IF((DIS.gt.(J-1)*Band_Width).and.(DIS.le.(J*Band_Width))) then	
			Band_ID(I)=J
		End IF
	End Do
	
	Surface_Angle=acos(XR/Real(Dis))
	Surface_Angle=Surface_Angle*(180.d0/PIE)
	If (Surface_Angle.lt.(0.d0)) Surface_Angle=Surface_Angle+360.d0
	If (Surface_Angle.ge.(360.d0)) Surface_Angle=Surface_Angle-360.d0
	If (YR .lt. 0.d0) Surface_Angle=360.d0-Surface_Angle

	Do K=1,N_Section
		IF((Surface_Angle.GT.((K-1)*Increment)).AND.(Surface_Angle.LE.(K*Increment))) then
			Section_ID(I)=K
		End IF
	End Do
End Do


H_Cut_Off=6.d0
Point_Count=0
	Do K=1,N_Section
		Flag_Exit=0
		Do J=1,N_Band
			I=1
			Do While((Flag_Exit==0).AND.(I.LE.TM))
				IF((Section_ID(I)==K).AND.(Band_ID(I)==J)) Then
					IF((H(I).LE. H_Cut_Off)) then
						Point_Count=Point_Count+1
						Flag_Bowl(I)=1
					Else
						Flag_Exit=1	
					End IF
				End IF
			I=I+1
			End Do
		End Do
	End Do

Membrane_Surface_Area=0.d0
Whole_Surface_Area=0.d0	
Do I=1,TM	
	IF((Flag_Bowl(I)==1).and.(Flag_6(I)==0)) then
		Do J=1,6
			p1=I
			p2=Neigh_right(I,J)
			if(j.ne.6)	p3=Neigh_right(I,J+1)

			IF(j==6) p3=Neigh_right(I,1)

			Vx_p1p2=rx(p2)-rx(p1)
			Vy_p1p2=ry(p2)-ry(p1)
			Vz_p1p2=rz(p2)-rz(p1)

			Vx_p1p3=rx(p3)-rx(p1)
			Vy_p1p3=ry(p3)-ry(p1)
			Vz_p1p3=rz(p3)-rz(p1)

			Ax=Vy_p1p2*Vz_p1p3 - Vy_p1p3*Vz_p1p2
			Ay=(Vx_p1p2*Vz_p1p3 - Vx_p1p3*Vz_p1p2)*(-1)
			Az=Vx_p1p2*Vy_p1p3-Vx_p1p3*Vy_p1p2

			Area=0.5*sqrt(Ax*Ax + Ay*Ay + Az*Az)
			Membrane_Surface_Area=Membrane_Surface_Area+Area
		End Do
	End IF	
	
	
	IF(Flag_6(I)==0) then
		Do J=1,6
			p1=I
			p2=Neigh_right(I,J)
			if(j.ne.6)	p3=Neigh_right(I,J+1)

			IF(j==6) p3=Neigh_right(I,1)

			Vx_p1p2=rx(p2)-rx(p1)
			Vy_p1p2=ry(p2)-ry(p1)
			Vz_p1p2=rz(p2)-rz(p1)

			Vx_p1p3=rx(p3)-rx(p1)
			Vy_p1p3=ry(p3)-ry(p1)
			Vz_p1p3=rz(p3)-rz(p1)

			Ax=Vy_p1p2*Vz_p1p3 - Vy_p1p3*Vz_p1p2
			Ay=(Vx_p1p2*Vz_p1p3 - Vx_p1p3*Vz_p1p2)*(-1)
			Az=Vx_p1p2*Vy_p1p3-Vx_p1p3*Vy_p1p2

			Area=0.5*sqrt(Ax*Ax + Ay*Ay + Az*Az)
			Whole_Surface_Area=Whole_Surface_Area+Area
		End Do
	End IF	
	
End Do	
	
write(333,*)step,Ben,VD,Point_Count,Point_Count/Real(TM),Membrane_Surface_Area/Whole_Surface_Area
Return
End Subroutine Bulge_Size


subroutine cut_membrane_ini
use virus_uptake
implicit none
integer :: i

		Virus_Raidus=5.d0
		Do I=1,TM
			Flag_in_area(I)=0
		End Do

		Do I=1,TM
			XR=RX0(I)-RX0(layer+1)
			YR=RY0(I)-RY0(layer+1)
			dis=sqrt(xr*xr + yr*yr)
			IF(dis.le. (dyara*Virus_Raidus)) then
				Flag_in_area(I)=1
			End IF
		End Do
return
end subroutine cut_membrane_ini



subroutine cut_membrane
use virus_uptake
implicit none
integer :: i
		Virus_Raidus=5.d0
		Do I=1,TM
			Flag_in_area(I)=0
		End Do

		Do I=1,TM
		
			XR=RX(I)-RX(layer+1)
			YR=RY(I)-RY(layer+1)
			
			dis=sqrt(xr*xr + yr*yr)
			IF(dis.le. (dyara*Virus_Raidus)) then
				Flag_in_area(I)=1
			End IF
		End Do
return
end subroutine cut_membrane




Subroutine Virus_Wrapping
use virus_uptake
implicit none
integer :: i,j,k,jj,kk
integer,dimension(Comp_Spike) :: Wrap_Z,Wrap_Virus

REAL*8,dimension(Total_Virus_Particle)    :: CX,CY,CZ
Real*8,dimension(Total_Virus_Particle)   :: Theta,Phi
Real*8 :: CX_center,CY_center,CZ_center 
Real*8 :: Raidus_Org,R_2D,Raidus_Org_Avg

Real*8 :: Theta_Max,Theta_min,Theta_Inc
Real*8 :: Phi_Max,Phi_min,Phi_Inc
integer :: Theta_Sec,Phi_Sec,Total_Sec

real*8,allocatable,dimension(:,:) :: Area_Mesh
real*8 :: Total_SA
real*8 :: phi_low,phi_up,theta_low,theta_up	
integer :: Count_Stick
real*8 :: detla_theta,detla_phi

r_vim_covid=4.d0 ! Spike and ECC connnected

!__Calculate the number of Spikes-Vimentin connections __!
		Do i=1,N_Spike
			Spike_occupy(i)=0
		End do

		Do K_Vim=1,NVF
			IF(K_Vim==NN) cycle
			i=TM+(K_Vim*SVN)
			Do K_Spike=1,N_Spike
				j=Cl+Sp+(K_spike*SS)
				IF(j==(comp+ss)) cycle

				xr=Rx(i)-Rx(j)          
				yr=Ry(i)-Ry(j)
				zr=Rz(i)-Rz(j)          

				dis= xr*xr + yr*yr + zr*zr

				If(dis.lt.r_vim_covid) then
					Spike_occupy(K_Spike)=1
				End if	
			End do
		End do
!_________________________________________________________!


! Trace back the spikes to the virus surface and this what is covered !
!___________Find particles on the virus surface______________!

	Do I=1,Comp_Spike
		Wrap_Z(I)=0
	End Do
Count_Stick=0
	Do K_Spike=1,N_Spike
		IF(Spike_occupy(K_Spike)==1) then
		Count_Stick=Count_Stick+1
			jj=Cl+Sp+(K_spike*SS)
				Do I=Comp_Edge_Sphere+1,Comp_Edge_complete
					IF(JJ==E2(I))	then
						kk=E1(I-1)
						Wrap_Z(KK)=1
					End IF	
				 End Do

		End IF
	End Do	
	
!_________________________________________________________!


! Find the center of mass of the virus surface !
	 CX_center=0.d0
	 CY_center=0.d0
	 CZ_center=0.d0

	K=1
	Do I=CL+1,Comp
			CX(K)=RX(I)		
			CY(K)=RY(I)
			CZ(K)=RZ(I)
		
			Wrap_Virus(K)=Wrap_Z(I)
		
			CX_center=CX_center+CX(K)
			CY_center=CY_center+CY(K)
			CZ_center=CZ_center+CZ(K)

			K=K+1
	End Do
	
	
	
 CX_center=CX_center/Real(Total_Virus_Particle)
 CY_center=CY_center/Real(Total_Virus_Particle)
 CZ_center=CZ_center/Real(Total_Virus_Particle)

!_________________________________________________________!


! Give all points a theta and phi value !
! Move the center of mass to the origin !
	Do I=1,Total_Virus_Particle
		 CX(I)=CX(I)-CX_center
		 CY(I)=CY(I)-CY_center
		 CZ(I)=CZ(I)-CZ_center
	End Do
	
		Do I=1,Total_Virus_Particle
			Theta(I)=0.d0
			Phi(I)=0.d0
		End Do

		Raidus_Org_Avg=0.d0
		DO I=1,Total_Virus_Particle
			Raidus_Org= sqrt(cx(I)*cx(I)+cy(I)*cy(I)+cz(I)*cz(I))
			Phi(I)   = acos(CZ(I)/Raidus_Org)
			
			Raidus_Org_Avg=Raidus_Org_Avg+Raidus_Org
				
			R_2D=sqrt(CX(I)*CX(I) + CY(I)*CY(I))

			If((CX(I).GT.(0.d0)).AND.(CY(I).GE.(0.d0))) Then
				Theta(I) = asin(CY(I)/R_2D)
			ElseIf((CX(I).LE.(0.d0)).AND.(CY(I).GT.(0.d0))) Then
				Theta(I) =PIE- asin(CY(I)/R_2D)
			ElseIf((CX(I).LT.(0.d0)).AND.(CY(I).LE.(0.d0))) Then
				Theta(I) =PIE- asin(CY(I)/R_2D)
			ElseIf((CX(I).GE.(0.d0)).AND.(CY(I).LT.(0.d0))) Then
				Theta(I) =(2.d0*PIE)+ asin(CY(I)/R_2D)
			End IF
		End Do
		
		Raidus_Org_Avg=Raidus_Org_Avg/Real(Total_Virus_Particle)
!_________________________________________________________!

Total_SA=0.d0
		Do K_Spike=1,N_Spike
			jj=Cl+Sp+(K_spike*SS)
			Do I=Comp_Edge_Sphere+1,Comp_Edge_complete
				IF(JJ==E2(I))	then
					kk=E1(I-1)
					kk=kk-Cl
				if((kk.ge.6).and.(kk.le.996)) then
 detla_theta=(theta(kk+2)+ 0.5*(theta(kk+3)-theta(kk+2)))-(theta(kk-2)+0.5*(theta(kk-2)-theta(kk-3)))
 detla_phi=cos(phi(kk+2)+ 0.5*(phi(kk+3)-phi(kk+2)))-cos(phi(kk-2)+0.5*(phi(kk-2)-phi(kk-3)))
					 Total_SA=Total_SA+ detla_theta*detla_phi*(-1)*Raidus_Org_Avg*Raidus_Org_Avg
					 end if
				End IF	
			End Do
		End Do	


! Section the whole sphere in small rectangular patches!
! Find the weighted area !
! Find which virus points comes under the window!
! Add all those area patches!


	Theta_Max=360.d0
	Theta_Min=0.d0
	Theta_Inc=7.25*2.255


	Phi_Max=180.d0
	Phi_Min=0.d0
	Phi_Inc=7.25*2.255
	
	Theta_Sec=nint(Theta_Max/Theta_Inc)
	Phi_Sec=nint(Phi_Max/Phi_Inc)
	Total_Sec=Theta_Sec*Phi_Sec

allocate(Area_Mesh(Phi_Sec,Theta_Sec))

	
Total_SA=0.d0
K=1
Do While(K.lt.997)
	IF (Wrap_Virus(K)==1) then

Do I=1,Phi_Sec ! Phi loop
		Phi_low=(I-1)*Phi_Inc*(PIE/180.d0)
		Phi_Up=(I)*Phi_Inc*(PIE/180.d0)
			Do J=1,Theta_Sec ! Theta loop
			Theta_low=(J-1)*Theta_Inc*(PIE/180.d0)
			Theta_Up=(J)*Theta_Inc*(PIE/180.d0)
			
					IF((Phi(K).GE.Phi_Low).and.(Phi(K).LT.Phi_Up)) Then
						IF((Theta(K).GE.Theta_Low).and.(Theta(K).LT.Theta_Up)) Then
						Area_Mesh(I,J)=( (Theta_Up-Theta_low)*(cos(Phi_Up)-cos(Phi_Low))*(-1.d0) )*Raidus_Org_Avg*Raidus_Org_Avg
						Total_SA=Total_SA+Area_Mesh(I,J)
						End IF
					End IF
			
			End Do
		End Do
		End IF		
k=k+5	
	End Do
	write(120,*)step,Ben,VD,Total_SA/(4*pie*Raidus_Org_Avg*Raidus_Org_Avg),Count_Stick
return
End Subroutine Virus_Wrapping

!_________________________________________________________!





Subroutine Read_Data
use virus_uptake
implicit none
integer :: I

	DO i = 1,Comp_Spike
		READ(100) RX(i),RY(i),RZ(i)
	END DO 

return
end subroutine Read_Data


Subroutine Make_File_Names
use virus_uptake
implicit none


	write(Ben_Ch,'(i0)') Ben
	write(Sp_no_Ch,'(i0)') Sp_no
	write(VD_Ch,'(i0)') VD
	write(Ens_Ch,'(i0)') Ens

	filename1='Ini_seed'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename2='Tri_Nodes'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename3='Tri_Edges'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename4='Tri_Faces'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename5='Triangle_data'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename6='Tri_Neigh'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename7='Tri_Neigh_Sec'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename8='Edges'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename110='Mem'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename13='Time_Wrapping'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename14='Time_Membrane_Area'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	filename16='Time_Velocity'//trim(adjustl(Ben_ch))//'_'//trim(adjustl(Sp_no_Ch))//'_' &
	   //trim(adjustl(VD_Ch))//'_'//trim(adjustl(Ens_Ch))//'.dat'
	
	
	Open(unit=45,file=filename1)

	Open(unit=10,file=filename2)
	Open(unit=20,file=filename3)
	Open(unit=30,file=filename4)
	Open(unit=40,file=filename5)
	Open(unit=50,file=filename6)
	Open(unit=60,file=filename7)

	Open(unit=70,file=filename8)

	open(unit=100,file=filename110,form="unformatted")
	open(unit=120,file=filename13)
	open(unit=333,file=filename14)
	open(unit=533,file=filename16)
	
return
end subroutine Make_File_Names

