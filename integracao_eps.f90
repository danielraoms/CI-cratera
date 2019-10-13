module integracao_eps
	use parametros
	implicit none

contains

!subroutine para evolução das partículas MÓVEIS 
subroutine integracao_verlet()
	integer :: j !índice para arrays

	!não evolua as partículas das paredes
	do j = paredes+1, N

		!translacional - horizontal (Verlet quarta ordem)
		xnewer(j) = 2.0d0*xnew(j) - xold(j) + (forcax(j)/m(j))*(dt**2.0d0)!xold(j) + vxold(j)*dt + (forcax(j)/m(j))*(dt**2.0d0)/2.0d0
		vxnew(j) = (xnewer(j) - xnew(j))/dt !vxold(j) + (forcax(j)/m(j))*dt

		xold(j) = xnew(j)
		xnew(j) = xnewer(j)
		vxold(j) = vxnew(j)

		!translacional - vertical (Verlet quarta ordem)
		ynewer(j) = 2.0d0*ynew(j) - yold(j) + (forcay(j)/m(j))*(dt**2.0d0) !yold(j) + vyold(j)*dt + (forcay(j)/m(j))*(dt**2.0d0)/2.0d0
		vynew(j) = (ynewer(j) - ynew(j))/dt!vyold(j) + (forcay(j)/m(j))*dt

		yold(j) = ynew(j)
		ynew(j) = ynewer(j)
		vyold(j) = vynew(j)

		!tangencial (Euler)
		omega_new(j) = omega_old(j) + (torque(j)/inertia(j))*dt
		theta_new(j) = theta_old(j) + omega_new(j)*dt
	
		theta_old(j) = theta_new(j)
		omega_old(j) = omega_new(j)

	end do

	return
end subroutine

!subroutine para evolução das partículas MÓVEIS 
subroutine integracao_verlet_cratera()
	integer :: j !índice para arrays

	!não evolua as partículas das paredes
	do j = paredes+1, N
		if ((flag_dig(j) .eq. 1)) then
		!translacional - horizontal (Verlet quarta ordem)
		xnewer(j) = 2.0d0*xnew(j) - xold(j) + (forcax(j)/m(j))*(dt**2.0d0)!xold(j) + vxold(j)*dt + (forcax(j)/m(j))*(dt**2.0d0)/2.0d0
		vxnew(j) = (xnewer(j) - xnew(j))/dt !vxold(j) + (forcax(j)/m(j))*dt

		xold(j) = xnew(j)
		xnew(j) = xnewer(j)
		vxold(j) = vxnew(j)

		!translacional - vertical (Verlet quarta ordem)
		ynewer(j) = 2.0d0*ynew(j) - yold(j) + (forcay(j)/m(j))*(dt**2.0d0) !yold(j) + vyold(j)*dt + (forcay(j)/m(j))*(dt**2.0d0)/2.0d0
		vynew(j) = (ynewer(j) - ynew(j))/dt!vyold(j) + (forcay(j)/m(j))*dt

		yold(j) = ynew(j)
		ynew(j) = ynewer(j)
		vyold(j) = vynew(j)

		!tangencial (Euler)
		omega_new(j) = omega_old(j) + (torque(j)/inertia(j))*dt
		theta_new(j) = theta_old(j) + omega_new(j)*dt
	
		theta_old(j) = theta_new(j)
		omega_old(j) = omega_new(j)

		end if
	end do

	return
end subroutine

!subroutine para gerar interface gráfica .eps da configuração atual do SISTEMA no decorrer da SIMULAÇÃO
subroutine salva_eps(ttwrite,yplaca,parts,parede,part_raio,part_pos_x,part_pos_y,arq,part_ang,force,vel) !cria arquivo eps com as partículas desenhadas
  integer, intent (in) :: ttwrite, parts, parede
  double precision, intent(in) :: yplaca
  double precision, dimension (1:N), intent (in) :: part_raio
  double precision, dimension (1:N), intent (in) :: part_pos_x
  double precision, dimension (1:N), intent (in) :: part_pos_y
  double precision, dimension (1:N), intent (in) :: part_ang
  double precision, dimension (1:N), intent (in) :: force
  double precision, dimension (2,N), intent (in) :: vel
  double precision, dimension (1:N) :: arrow_size
  integer, intent (in) :: arq
  integer :: l, j, k
  double precision :: max_force, max_vel, nx, ny
  double precision, dimension (N,3) :: rgb
  double precision, dimension (1:N) :: linewidth
  real*8 :: scale = 400.0d0 !para raiomed = 0.1

  character (len = 20) :: filename1

  !create output value - number entre 000 e 999!!!
  if(ttwrite .ne. 0) then !Se diferente de 0, arquivo de CI
    if(arq .eq. 0) then
        filename1 = "outputsxy-0000.eps"
    else
        filename1 = "outputixy-000.eps"
    end if
  filename1(14:14)=CHAR(48+mod(ttwrite,10))
  filename1(13:13)=CHAR(48+(mod(ttwrite,100)/10))
  filename1(12:12)=CHAR(48+(mod(ttwrite,1000)/100))
  filename1(11:11)  =CHAR(48+ttwrite/1000)
  else
     if(arq.eq.0) then
        filename1="outputsxy-000.eps"
     else
        filename1="outputixy-000.eps"
     end if
  end if

  open(unit=210,file=filename1,status='unknown')

  !Writes initial head of eps
  write(210,90)
  write(210,91)
  write(210,92) filename1
  write(210,93)
  write(210,94)
  write(210,95)
  write(210,96) 0, 0, 3500, 2500
  write(210,97)
  write(210,98)
  write(210,99)
  write(210,100)
  write(210,101)
  write(210,102)
  write(210,103)

  !identificando a força máxima
  max_force = 0.0d0
  do l = 1, N
	max_force = max(max_force,force(l))
  end do

  !identificando a velocidade máxima
  max_vel = 0.0d0
  do l = 1, N
	max_vel = max(max_vel,norm2(vel(:,l)))
  end do


  !normalização das forças
  do l = 1, N
  	rgb(l,1) = 1.0d0
	if (max_force .gt. 0.0d0) then
		rgb(l,2) = 1.0d0 - (force(l)/max_force)
	else
		rgb(l,2) = 1.0d0
	end if
	rgb(l,3) = 0.0d0
	linewidth(l) = 1.0d0 + 2.5d0*(force(l)/max_force)

	!write(*,*) "eps", l, force(l), max_force, rgb(l,2)
	!read(*,*)
  end do

!  !desenhe as partículas da parede
!   do l = 1, parede
!     write(210,104) nint(scale*part_pos_x(l)), nint(scale*part_pos_y(l)), nint(scale*part_raio(l))
!     write(210,*) 'gsave'
!     write(210,*) '0.5 setgray'
!     write(210,*) 'fill'
!     write(210,*) 'grestore'
!     write(210,*) 'stroke'
!  end do

  !parte do scheme de validar atrito
  !desenha a placa plana inferior
 ! write(210,*) 'newpath'
 ! write(210,108) int(scale*0.5d0), int(scale*(0.5d0 + 2.0d0*0.015d0))
 ! write(210,109) 1815, 0 
 ! write(210,*) 'stroke'

  !desenhe as partículas móveis
  do l=1,parts!parede+1,parts
     write(210,*) 'gsave'
     write(210,*) 0.08, 'setlinewidth' 
     !desenha a borda do círculo de uma partícula
     write(210,104) nint(scale*part_pos_x(l)), nint(scale*part_pos_y(l)), nint(scale*part_raio(l))
     if (l .gt. parede) then
        write(210,108) int(scale*part_pos_x(l)), int(scale*part_pos_y(l)) 
     	write(210,109) int(scale*part_raio(l)*cos(part_ang(l))), int(scale*part_raio(l)*sin(part_ang(l)))
     end if
     !salva o desenho do círculo
     write(210,*) 'gsave'
     write(210,*) 'stroke'
  end do

  !desenhe o eixo das ordenadas
     write(210,108) int(scale*part_pos_x(1)), int(scale*part_pos_y(1))
     write(210,*) 3.0, 'setlinewidth'
     write(210,109) 0, int(scale*part_pos_y(parede))

     write(210,*) 'gsave'

  !desenhe o eixo das abscissas
     write(210,108) int(scale*part_pos_x(1)), int(scale*part_pos_y(1))
     write(210,*) 3.0, 'setlinewidth'
     write(210,109) int(scale*part_pos_x(parede)), 0

    write(210,*) 'gsave'
    write(210,*) 'stroke'

  !desenhe o grid
    !grids horizontais
    l = 1
    do while (k .lt. (bw + lw))
    	 k = (bw + 4) + 5*(l-1) !começa em bw+4, incrementa a cada 5 partículas 
    	 write(210,108) int(scale*part_pos_x(1)), int(scale*part_pos_y(k))
 	 write(210,109) int(scale*part_pos_x(parede)), 0
	 write(210,*) 0.4, 'setlinewidth'
         write(210,110) 0.58, 0.58, 0.58
	
	 write(210,*) 'gsave'
	 write(210,*) 'stroke'

	 l = l + 1
    end do


    !grids verticais
    l = 1
    do while (k .lt. bw)
    	 k = 5 + 5*(l-1) !começa em 5, incrementa a cada 5 partículas 
    	 write(210,108) int (scale*part_pos_x(k)), int(scale*part_pos_y(1))
	 write(210,109) 0, int(scale*part_pos_y(parede))
	 write(210,*) 0.4, 'setlinewidth'
	 write(210,110) 0.58, 0.58, 0.58

         write(210,*) 'gsave'
         write(210,*) 'stroke'
	 l = l + 1
    end do

    write(210,*) 'stroke'


  write(210,105)
  write(210,106)
  write(210,107)
  close(unit=210)

90 format('%%!PS-Adobe-3.0 EPSF-3.0')
91 format('%%Document-Fonts: Times-Roman')
92 format('%%Title: ', A15)
93 format('%%Creator: YDS')
94 format('%%CreationDate: unknown')
95 format('%%Pages: 1')
96 format('%%BoundingBox: ',4I5)
97 format('%%LanguageLevel: 1')
98 format('%%EndComments')
99 format('%%BeginProlog')
100 format('%%EndProlog')
101 format('0.0000 0.0000 0.0000 setrgbcolor')
102 format('%% Page:     1    1')
103 format('save')
104 format(3I8,'  0   360  arc closepath')
105 format('restore showpage')
106 format('%%Trailer')
107 format('%%EOF')
108 format(2I8,' moveto')
109 format(2I8, ' rlineto')
110 format(3F10.3,' setrgbcolor')
111 format(1I8,'setlinewidth')


end subroutine salva_eps

!subroutine para gerar interface gráfica .eps da configuração atual do SISTEMA no decorrer da SIMULAÇÃO
subroutine salva_eps_cratera(ttwrite,yplaca,parts,parede,part_raio,part_pos_x,part_pos_y,arq,part_ang,force,vel,bandeira_dig) !cria arquivo eps com as partículas desenhadas
  integer, intent (in) :: ttwrite, parts, parede
  double precision, intent(in) :: yplaca
  integer, dimension (1:N), intent (in) :: bandeira_dig
  double precision, dimension (1:N), intent (in) :: part_raio
  double precision, dimension (1:N), intent (in) :: part_pos_x
  double precision, dimension (1:N), intent (in) :: part_pos_y
  double precision, dimension (1:N), intent (in) :: part_ang
  double precision, dimension (1:N), intent (in) :: force
  double precision, dimension (2,N), intent (in) :: vel
  double precision, dimension (1:N) :: arrow_size
  integer, intent (in) :: arq
  integer :: l, j
  double precision :: max_force, max_vel, nx, ny
  double precision, dimension (N,3) :: rgb
  double precision, dimension (1:N) :: linewidth
  real*8 :: scale = 800.0d0 !para raiomed = 0.1

  character (len = 20) :: filename1

  !create output value - number entre 000 e 999!!!
  if(ttwrite .ne. 0) then !Se diferente de 0, arquivo de CI
    if(arq .eq. 0) then
        filename1 = "outputcxy-0000.eps"
    else
        filename1 = "outputixy-000.eps"
    end if
  filename1(14:14)=CHAR(48+mod(ttwrite,10))
  filename1(13:13)=CHAR(48+(mod(ttwrite,100)/10))
  filename1(12:12)=CHAR(48+(mod(ttwrite,1000)/100))
  filename1(11:11)  =CHAR(48+ttwrite/1000)
  else
     if(arq.eq.0) then
        filename1="outputsxy-000.eps"
     else
        filename1="outputixy-000.eps"
     end if
  end if

  open(unit=210,file=filename1,status='unknown')

  !Writes initial head of eps
  write(210,90)
  write(210,91)
  write(210,92) filename1
  write(210,93)
  write(210,94)
  write(210,95)
  write(210,96) 0, 0, 3000, 2500
  write(210,97)
  write(210,98)
  write(210,99)
  write(210,100)
  write(210,101)
  write(210,102)
  write(210,103)

  !identificando a força máxima
  max_force = 0.0d0
  do l = 1, N
	max_force = max(max_force,force(l))
  end do

  !identificando a velocidade máxima
  max_vel = 0.0d0
  do l = 1, N
	max_vel = max(max_vel,norm2(vel(:,l)))
  end do


  !normalização das forças
  do l = 1, N
  	rgb(l,1) = 1.0d0
	if (max_force .gt. 0.0d0) then
		rgb(l,2) = 1.0d0 - (force(l)/max_force)
	else
		rgb(l,2) = 1.0d0
	end if
	rgb(l,3) = 0.0d0
	linewidth(l) = 1.0d0 + 2.5d0*(force(l)/max_force)

	!write(*,*) "eps", l, force(l), max_force, rgb(l,2)
	!read(*,*)
  end do

!  !desenhe as partículas da parede
!   do l = 1, parede
!     write(210,104) nint(scale*part_pos_x(l)), nint(scale*part_pos_y(l)), nint(scale*part_raio(l))
!     write(210,*) 'gsave'
!     write(210,*) '0.5 setgray'
!     write(210,*) 'fill'
!     write(210,*) 'grestore'
!     write(210,*) 'stroke'
!  end do

  !parte do scheme de validar atrito
  !desenha a placa plana inferior
 ! write(210,*) 'newpath'
 ! write(210,108) int(scale*0.5d0), int(scale*(0.5d0 + 2.0d0*0.015d0))
 ! write(210,109) 1815, 0 
 ! write(210,*) 'stroke'

  !desenhe as partículas móveis
  do l=1,parts!parede+1,parts
     if  (bandeira_dig(l) .eq. 1) then
     write(210,*) 'gsave'
     write(210,*) 0.08, 'setlinewidth' 
     !desenha a borda do círculo de uma partícula
     write(210,104) nint(scale*part_pos_x(l)), nint(scale*part_pos_y(l)), nint(scale*part_raio(l))
     if (l .gt. parede) then
        write(210,108) int(scale*part_pos_x(l)), int(scale*part_pos_y(l)) 
     	write(210,109) int(scale*part_raio(l)*cos(part_ang(l))), int(scale*part_raio(l)*sin(part_ang(l)))
     end if
     !salva o desenho do círculo
     write(210,*) 'gsave'
     write(210,*) 'stroke'
     end if
  end do

  write(210,105)
  write(210,106)
  write(210,107)
  close(unit=210)

90 format('%%!PS-Adobe-3.0 EPSF-3.0')
91 format('%%Document-Fonts: Times-Roman')
92 format('%%Title: ', A15)
93 format('%%Creator: YDS')
94 format('%%CreationDate: unknown')
95 format('%%Pages: 1')
96 format('%%BoundingBox: ',4I5)
97 format('%%LanguageLevel: 1')
98 format('%%EndComments')
99 format('%%BeginProlog')
100 format('%%EndProlog')
101 format('0.0000 0.0000 0.0000 setrgbcolor')
102 format('%% Page:     1    1')
103 format('save')
104 format(3I8,'  0   360  arc closepath')
105 format('restore showpage')
106 format('%%Trailer')
107 format('%%EOF')
108 format(2I8,' moveto')
109 format(2I8, ' rlineto')
110 format(3F10.3,' setrgbcolor')
111 format(1I8,'setlinewidth')


end subroutine salva_eps_cratera

end module integracao_eps
