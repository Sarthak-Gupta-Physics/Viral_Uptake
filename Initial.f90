! Author: Sarthak Gupta !
! Please contact sg207@rice.edu for any questions!

! This code generates the initial conditions of the viral uptake simulations !

module filenames
implicit none
character*50    :: filename1,filename2,filename3,filename4,filename5
character*50    :: filename6,filename7,filename8
Real*8  :: Ben,Sp_No,VD,Ens
end module filenames

! Initilazition of Triangulated Mesh !
! Gives: 
! Number of nodes
! Number of Edges
! Number of Triangles
! Triangle Neighbours

! Lx, Ly cannot be smaller than space^2
! or whatever is the distance^2 between
! two points on x-axis only


Program Membrane_Initial
use filenames
implicit none
! Keep Lx/space an even number, it calculates the POS_ACE2 correctly.
integer,parameter :: Lx=57     ! Dimension of membrane in x-direction.
integer,parameter :: Ly=Lx     ! Dimension of membrane in y-direction.
Real*8,parameter :: space=1.5  ! 
Real*8,allocatable,Dimension(:) :: X,Y,Z
Real*8 :: ran2

Integer :: Seed
Real*8  :: Mx,My
Integer :: SVN,NE_Comp,Comp_Spike
Integer :: PosV,NVF,TM,Vim_Edge
Integer :: Interval_VIM,NN
Integer :: V,ACE2,POS_ACE2

Integer :: Layer,CL,Comp,Sp,Sp_Edge,SE1,SE2
Real*8	:: Sigma_VIM
Real*8 :: xr,yr
Real*8 :: dis,error
Real*8 :: E1,E2,E3
Real*8 :: s,Area,AI
Real*8 :: dx1,dx2,dx3,dy1,dy2,dy3
Real*8 :: Pie

Integer :: I,J,K
Integer :: V1,V2,V3
Integer :: NT,NE,NSN,NE_Spike
Integer,allocatable,Dimension(:) :: P1,P2,P3,N1,N2,N3

Integer,allocatable,Dimension(:) :: Flag_Vim,Shuffle_1,Shuffle_2

real*8  :: Sx,Sy,Sz,R_2D
real*8,allocatable,Dimension(:) :: Theta,Phi
Integer :: NS,NM_Spike

Integer :: N_Spike,Interval_Spike,SS,Raidus
Real*8  :: Raidus_Org

Real*8  :: Temp,Frac,Temp_Seed

Integer :: Flag_Vim_Num


Real*8  :: Up


integer :: num_args
character(len=12), dimension(:), allocatable :: args


Up=12.d0 ! Shifting Z to up
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


Call init_random_seed()
Call Random_number(Temp_Seed)
seed=Temp_Seed*(-1000000)
write(45,*)seed,Ben,Sp_no,VD,Ens

TM = int(Lx*Lx/(space*space)) ! Total Number of Membrane monomer
Mx=int(Lx/Space)		      ! Number of monomers in one X-axis.
My=int(Ly/Space)		      ! Number of monomers in one Y-axis.
SVN  = 4			      ! No. of monomers in One Vimentin
ACE2 = 4			      ! No. of monomers in ACE2 receptor
POS_ACE2 = (Mx/2)*(My-1)	      ! Position of ACE at the center of Membrane
Interval_VIM=1	 		      ! 1/frequency of Vimentin 
!Basically number of spots on X-axis*No. of total rows.
Frac=VD/100
NVF=(Frac*(TM))


V=NVF*SVN	 		      ! Total number of vimentin monomers
PosV=Interval_VIM		      ! Frequency of Position of Vimentin
Layer=TM+V			      ! Membrane + Vimentin
 CL = Layer + ACE2		      ! Membrane + Vimentin + Ace2

open(unit=1000,file="Sphere_Data.dat")
read(1000,*) Sp,Sp_Edge
 Comp= CL + Sp			! Membrane + Vimentin + Ace2 + Sphere

SS=2			  ! Size of Spike.

N_spike=Sp_no!Sp/Interval_Spike ! Number of spike protein.
Interval_Spike=Sp/N_spike	  ! Interval of spikes on the virus.


Comp_Spike=Comp+(N_Spike*SS) ! Membrane + Vimentin + Ace2 + Sphere + Spikes



allocate(X(Comp_Spike),Y(Comp_Spike),Z(Comp_Spike))
allocate(Theta(Comp_Spike),Phi(Comp_Spike))
allocate(Flag_Vim(TM),Shuffle_1(TM))

PIE=4.D0*DATAN(1.D0)    ! Just pie, you know.


!Initilazation of arrays!
Do I=1,Comp_Spike
	X(I)=0.d0
	Y(I)=0.d0
	Z(I)=0.d0	
End do
!!!!!!!!!!!!!!!!!!!!!!!! 






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Giving values to arrays      !
! Position of nodes in membrane   !
! This gives you Up-Down triangle !
!    Equillateral triangles       !

Do J=1,int(Ly/space)
	Do I=1,int(Lx/space)
		Y((int(Ly/space)*(J-1))+I)=(J-1)*(space*sin(Pie/3.d0)) + 1
		If(mod(J,2).ne.0) then
			X((int(Lx/space)*(J-1))+I)=space*(I-1)+1
		Else
			X((int(Lx/space)*(J-1))+I)=space*(I-1)+1+(space/2.d0)
		End if	
	End do
End do



Do J=1,int(Ly/space)
	IF(J==1) Then
		Do I=(1+((J-1)*int(Ly/space)))+1,(J*int(Ly/space))-1
!			Write(120,*)I
		End Do
	End IF
	IF(J==int(Ly/space)) Then
		Do I=(1+((J-1)*int(Ly/space)))+1,(J*int(Ly/space))-1
!			Write(120,*)I
		End Do
	End IF
End do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Vimentin Brush      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Fisher–Yates shuffle   !

Do I=1,TM
	Flag_Vim(I)=0
	Shuffle_1(I)=I
End Do
	Flag_Vim(POS_ACE2)=0


Do I=TM,2,-1
	J=(ran2(seed))*I + 1
	IF((I.Eq.POS_ACE2).or.(J.Eq.POS_ACE2)) cycle
	Temp=Shuffle_1(I)
	Shuffle_1(I)=Shuffle_1(J)
	Shuffle_1(J)=Temp
End Do

Do I=1,NVF
	Flag_Vim(Shuffle_1(I))=1
End do
Flag_Vim(POS_ACE2)=0

!Print*,Flag_Vim(POS_ACE2),POS_ACE2,"dafs"

Flag_Vim_Num=0

Do J=1,TM
IF(Flag_Vim(J)==1) Then	
	Flag_Vim_Num=Flag_Vim_Num+1	
	Do I=TM+1+((Flag_Vim_Num-1)*SVN),TM+(Flag_Vim_Num*SVN)
		X(I)=X(J)
		Y(I)=Y(J)
		Z(I)=Z(J)+(I-(TM+(Flag_Vim_Num-1)*SVN))*Sigma_VIM
	End Do
End IF
End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ace2 Receptor	    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Do J=1,1
	Do I=Layer+1+((J-1)*ACE2),Layer+(J*ACE2)

		K=POS_ACE2
			! print*,I
			X(I)=X(K)
			Y(I)=Y(K)
			Z(I)=Z(K)+(I-(Layer+(J-1)*ACE2))*Sigma_VIM

	End Do
End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!



open(unit=2000,file="sphere_fibonacci_grid_n1000.xyz")



Do I=1,CL
	write(10,*) X(I),Y(I),Z(I)
End Do


DO I=Cl+1,Comp
	Read(2000,*) Sx,Sy,Sz
	x(I)=Sx
	y(I)=Sy
	z(I)=Sz
	Write(10,*)  Sx+X(POS_ACE2),Sy+Y(POS_ACE2),Sz+Up
End Do

!#########################################!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For putting Spikes on The Coronavirus             !
! Find all the angle where one wants to keep Spikes !
! 		Angles in Radians		      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!https://arxiv.org/pdf/1512.00424.pdf!

Do I=1,Comp_Spike
Theta(I)=0.d0
Phi(I)=0.d0
End Do

DO I=Cl+1,comp
	Raidus_Org= sqrt(x(I)*x(I)+y(I)*y(I)+z(I)*z(I))
	Phi(I)   = acos(Z(I)/Raidus_Org)

	R_2D=sqrt(X(I)*X(I) + Y(I)*Y(I))

	If((X(I).GT.(0.d0)).AND.(Y(I).GE.(0.d0))) Then
		Theta(I) = asin(Y(I)/R_2D)
	ElseIf((X(I).LE.(0.d0)).AND.(Y(I).GT.(0.d0))) Then
		Theta(I) =PIE- asin(Y(I)/R_2D)
	ElseIf((X(I).LT.(0.d0)).AND.(Y(I).LE.(0.d0))) Then
		Theta(I) =PIE- asin(Y(I)/R_2D)
	ElseIf((X(I).GE.(0.d0)).AND.(Y(I).LT.(0.d0))) Then
		Theta(I) =(2.d0*PIE)+ asin(Y(I)/R_2D)
	End IF

End Do


Do K=1,N_Spike ! Going over no. of Spikes
I=Cl+1+((K-1)*Interval_Spike) ! Position of Spike on the Sphere
	Do J=Cl+Sp+1+((K-1)*SS),Cl+Sp+(K*SS) ! Position of Spike particles

		Raidus=Raidus_Org+(J-(Cl+Sp+((K-1)*SS)))*Sigma_VIM
		X(J)=Raidus*Sin(Phi(I))*cos(Theta(I))
		Y(J)=Raidus*Sin(Phi(I))*sin(Theta(I))
		Z(J)=Raidus*cos(Phi(I))
		
		Write(10,*)  X(J)+X(POS_ACE2),Y(J)+Y(POS_ACE2),Z(J)+Up
	End Do
End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Edge Finder       !
!____ Membrane Edges ______!
error=0.01
NE=0  ! No of Edges
NE_Comp=0
Do I=1,TM-1
	Do J=I+1,I+int(Lx/space)+1
		xr=x(j)-x(i)
		yr=y(j)-y(i)

		dis=sqrt( xr*xr + yr*yr)

		IF(abs((Dis-Real(space))).le.error) then
			NE=NE+1
			write(20,*) I,J
			!Print*,I,J
			NE_Comp=NE_Comp+1
			write(70,*) I,J
		End IF
	End do
End do


!____ Vimentin Edges ______!
Vim_Edge=0
Flag_Vim_Num=0

Do J=1,TM
	IF(Flag_Vim(J)==1) Then	
		Flag_Vim_Num=Flag_Vim_Num+1	
		 I=TM+1+((Flag_Vim_Num-1)*SVN)
			
			NE_Comp=NE_Comp+1
			Vim_Edge=Vim_Edge+1
			write(70,*) J,I	
			write(20,*) J,I	

		Do I=TM+1+((Flag_Vim_Num-1)*SVN),TM+(Flag_Vim_Num*SVN)-1	
			 NE_Comp=NE_Comp+1
		  	 Vim_Edge=Vim_Edge+1
			 write(70,*) I,I+1
			 write(20,*) I,I+1	
		End Do
	End IF
End Do


!____ ACE2 Receptor ______!
i=POS_ACE2
j=Layer+1
NE_Comp=NE_Comp+1
write(70,*) I,J

Do i=Layer+1,CL-1
	j=I+1
	!print*,i,j
	NE_Comp=NE_Comp+1
	write(70,*) I,J
end do


!____ Sphere Edges ______!
Open(unit=3000,file="Sphere_Edge.dat")

	Do I=1,Sp_Edge
		Read(3000,*) SE1,SE2
		Write(20,*)  SE1+CL,SE2+CL
		NE_Comp=NE_Comp+1
		write(70,*) SE1+CL,SE2+CL
	End Do




!____ Spike Edges ______!
NE_Spike=0
!Finding Edges of Spikes
Do K=1,N_Spike ! Going over no. of Spikes
I=Cl+1+((K-1)*Interval_Spike) ! Position of Spike on the Sphere

	Do J=Cl+Sp+1+((K-1)*SS),Cl+Sp+1+((K-1)*SS) ! Position of Spike particles
		write(20,*) I,J
		NE_Spike=NE_Spike+1
		NE_Comp=NE_Comp+1
		!Print*,K,I,J
		write(70,*) I,J
	End Do

	Do J=Cl+Sp+1+((K-1)*SS),Cl+Sp+(K*SS)-1 ! Position of Spike particles
		write(20,*) J,J+1
		NE_Spike=NE_Spike+1
		NE_Comp=NE_Comp+1
		!Print*,K,J,J+1
		write(70,*) J,J+1
	End Do
End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Second Nearest Neighbour !
error=0.01
NSN=0  ! No of Edges
Do I=1,TM-1
	Do J=I+1,TM!I+int(Lx/space)+1

	xr=x(j)-x(i)
	yr=y(j)-y(i)

	dis=sqrt( xr*xr + yr*yr)

	IF(abs(Dis-(2.d0*space*sin(Pie/3.d0))).le.error) then
		NSN=NSN+1
		write(60,*) I,J
		!Print*,I,J

	End IF

	End do
End do
!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Triangle finder !!!!!!!
AI=sqrt(3.d0)*space*space/4.d0
NT=0 ! Number of Triangles
Do V1=1,TM-2
	Do V2=V1+1,TM-1
		Do V3=V2+1,TM
		!print*,V1,V2,V3
		dx1 = X(V1)-X(V2)
		dx2 = X(V2)-X(V3)
		dx3 = X(V3)-X(V1)
		
		dy1 = Y(V1)-Y(V2)
		dy2 = Y(V2)-Y(V3)
		dy3 = Y(V3)-Y(V1)

		E1= sqrt (dx1*dx1 + dy1*dy1)
		E2= sqrt (dx2*dx2 + dy2*dy2)
		E3= sqrt (dx3*dx3 + dy3*dy3)


		If( ((E1+E2+E3)-(3.d0*space)).le.error ) Then
			s=(E1+E2+E3)/2.d0
			Area=sqrt(s*(s-E1)*(s-E2)*(s-E3))
			IF(abs(Area-AI).le.error) then
				!print*,V1,V2,V3
				write(30,*) V1,V2,V3
				NT=NT+1
			End if
		End IF
		End Do
	End Do
End Do

 Close(30)
!!!!!!!!!!!!!!!!!!!!

write(40,*)TM,V,Layer,ACE2,CL,NE,NT,space,NSN,PosV,SVN,POS_ACE2,NVF,NN,NE_Comp,Lx,Comp_Spike,NE_Spike,SS,N_Spike,Vim_Edge

!! Triangle Neighbour Finder !!


allocate(P1(NT),P2(NT),P3(NT))
allocate(N1(NT),N2(NT),N3(NT))
Open(unit=30,file=filename4)

Do I=1,NT
	Read(30,*) P1(I),P2(I),P3(I)

	N1(I)=0
	N2(I)=0
	N3(I)=0
End Do




Do I=1,NT
	Do J=1,NT
	IF((J.gt.0).and.(J.le.NT).and.(I.ne.J)) Then
IF(((P1(I)==P1(J)).and.(P2(I)==P2(J))).or.((P1(I)==P2(J)).and.(P2(I)==P3(J))).or.((P1(I)==P1(J)).and.(P2(I)==P3(J)))) N1(I)=J

IF(((P2(I)==P1(J)).and.(P3(I)==P2(J))).or.((P2(I)==P2(J)).and.(P3(I)==P3(J))).or.((P2(I)==P1(J)).and.(P3(I)==P3(J)))) N2(I)=J

IF(((P1(I)==P1(J)).and.(P3(I)==P2(J))).or.((P1(I)==P2(J)).and.(P3(I)==P3(J))).or.((P1(I)==P1(J)).and.(P3(I)==P3(J)))) N3(I)=J
	End IF

	End do
End Do




Do I=1,NT
	write(50,*) N1(I),N2(I),N3(I)
End Do



End Program Membrane_Initial


Subroutine Make_File_Names
use filenames
implicit none

if( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i1,a1,i1,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"



elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i1,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i1,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i1,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i1,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i1,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i1,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i1,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i1,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i1,a1,i2,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i2,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i2,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i2,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i2,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i2,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i2,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i2,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i2,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"

elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i1,a1,i3,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<10).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i1,a1,i3,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i1,a1,i3,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i1,a1,i3,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i1,a1,i3,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i1,a1,i3,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i1,a1,i3,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i1,a1,i3,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i1,a1,i3,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
!___________________________________________________________________________________________________________!

elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i2,a1,i1,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"

elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i1,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i1,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i1,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i1,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i1,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i1,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i1,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i1,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i2,a1,i2,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i2,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i2,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i2,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i2,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i2,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i2,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i2,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i2,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i2,a1,i3,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"

elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"

elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<100).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i2,a1,i3,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i2,a1,i3,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i2,a1,i3,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i2,a1,i3,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i2,a1,i3,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i2,a1,i3,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i2,a1,i3,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i2,a1,i3,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
!______________________________________________________________________________________________________________!

elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i3,a1,i1,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"

elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<10).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i1,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i1,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i1,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i1,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i1,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i1,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i1,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i1,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"



elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i3,a1,i2,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"

elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<100).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i2,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i2,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i2,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i2,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i2,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i2,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i2,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i2,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"

elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<10) )then          ! Setting name for files.
WRITE(filename1,'(a8,i3,a1,i3,a1,i1,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i1,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i1,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i1,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i1,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i1,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i1,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i1,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i1,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i1,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i1,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i1,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i1,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i1,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i1,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i1,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<10).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i1,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i1,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i1,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i1,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i1,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i1,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i1,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i1,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i2,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i2,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i2,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i2,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i2,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i2,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i2,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i2,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i2,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i2,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i2,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i2,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i2,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i2,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i2,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i2,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<100).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i2,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i2,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i2,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i2,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i2,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i2,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i2,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i2,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"


elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<10) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i3,a1,i1,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i3,a1,i1,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i3,a1,i1,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i3,a1,i1,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i3,a1,i1,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i3,a1,i1,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i3,a1,i1,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i3,a1,i1,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<100) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i3,a1,i2,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i3,a1,i2,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i3,a1,i2,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i3,a1,i2,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i3,a1,i2,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i3,a1,i2,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i3,a1,i2,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i3,a1,i2,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
elseif( (int(Ben)<1000).and.(int(Sp_No)<1000).and.(int(VD)<1000).and.(int(Ens)<1000) )then          
WRITE(filename1,'(a8,i3,a1,i3,a1,i3,a1,i3,a4)')"Ini_seed",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename2,'(a9,i3,a1,i3,a1,i3,a1,i3,a4)')"Tri_Nodes",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename3,'(a9,i3,a1,i3,a1,i3,a1,i3,a4)')"Tri_Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename4,'(a9,i3,a1,i3,a1,i3,a1,i3,a4)')"Tri_Faces",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename5,'(a13,i3,a1,i3,a1,i3,a1,i3,a4)')"Triangle_data",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename6,'(a9,i3,a1,i3,a1,i3,a1,i3,a4)')"Tri_Neigh",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename7,'(a13,i3,a1,i3,a1,i3,a1,i3,a4)')"Tri_Neigh_Sec",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
WRITE(filename8,'(a5,i3,a1,i3,a1,i3,a1,i3,a4)')"Edges",int(Ben),"_",int(Sp_No),"_",int(VD),"_",int(Ens),".dat"
end if
!__________________________________________________________________________________________________!


Open(unit=45,file=filename1)
Open(unit=10,file=filename2)
Open(unit=20,file=filename3)
Open(unit=30,file=filename4)
Open(unit=40,file=filename5)
Open(unit=50,file=filename6)
Open(unit=60,file=filename7)
Open(unit=70,file=filename8)

return
end subroutine Make_File_Names



SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)

            DEALLOCATE(seed)
END SUBROUTINE



FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1

          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
END

