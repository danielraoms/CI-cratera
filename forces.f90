!calcula as forças elástica, viscosa e de atrito de uma colisão entre duas partículas MÓVEIS
subroutine all_forces(gama_n_dummy,mi_dummy,gama_s_dummy,mi_roll_dummy,E1,E2,vpois1,vpois2,G1,G2,m1,m2,a1,a2,x1,x2,&
		        y1,y2,vx1,vx2,vy1,vy2,o_1,o_2,Fx_elastica_parcial,Fy_elastica_parcial,Fx_viscosa_parcial,Fy_viscosa_parcial,&
		        Fs_x_dummy,Fs_y_dummy,Fs_dummy,T_rolling_dummy,normaldirectionx,normaldirectiony,dx_acum_x,dx_acum_y,sinal_v_rot) !,flag_dummy)
	use parametros, only: N, g, dt
	implicit none
	double precision, intent(in) :: gama_n_dummy, mi_dummy, gama_s_dummy, mi_roll_dummy !coeficientes elástico, viscoso, de atrito, regularização de atrito e rolling friction
	double precision, intent(in) :: E1, E2, vpois1, vpois2, G1, G2, m1, m2, a1, a2, x1, x2, y1, y2, vx1, vx2, vy1, vy2, o_1, o_2 !módulos de effective shear modulus, massas, raios, posições e velocidades nas componentes x e y das partículas p e q e velocidades angulares
	double precision, intent(out) :: Fx_elastica_parcial, Fy_elastica_parcial !componentes da força elástica
	double precision, intent(out) :: Fx_viscosa_parcial, Fy_viscosa_parcial !componentes da força viscosa
	double precision, intent(out) :: Fs_x_dummy, Fs_y_dummy !força de atrito
	double precision, intent(out) :: Fs_dummy 	!contribuição da força de atrito ao torque
	double precision, intent(out) :: T_rolling_dummy !torque de atrito de rolamento (rolling friction)
	double precision, intent(out) :: normaldirectionx, normaldirectiony !componentes da direção normal do contato
	double precision, intent(out) :: sinal_v_rot
	double precision, intent(inout) :: dx_acum_x, dx_acum_y
	double precision :: dx_t_trial_x, dx_t_trial_y
	double precision :: d 
	double precision :: d_con
	double precision :: vrelx, vrely, v_rel_normal
	double precision :: Fx_normal, Fy_normal, F_normal
	double precision :: v_rel_t, vt_x, vt_y, sinalx_v_rot, sinaly_v_rot, reg, nor
	double precision :: reg_x, reg_y, nor_x, nor_y
	double precision :: Fs_dummy_norma
	double precision :: csi, a1_c, a2_c, vx_rel_surf, vy_rel_surf, vx_rel_t, vy_rel_t
	double precision :: sinal_elasticax, sinal_elasticay, sinal_normalx, sinal_normaly
	double precision :: E_star, G_star, m_star, a_star, k_at, d_at, csi_ponto_t 
	double precision :: dx_trial_x_old, dx_trial_y_old
	double precision :: dx_r_x, dx_r_y, dx_trial_x, dx_trial_y, F_at_trial_x, F_at_trial_y, F_at_trial
	double precision :: dx_x, dx_y, t_trial_x, t_trial_y
	double precision :: normal_aux

	!distância entre as duas partículas
	d = sqrt((x2 - x1)**2.0d0 + (y2 - y1)**2.0d0)

	!direção normal do par de contato
	normaldirectionx = (x2 - x1)/d
	normaldirectiony = (y2 - y1)/d

	!componentes x e y das velocidades relativas das partículas
	vrelx = (vx2 - vx1)
	vrely = (vy2 - vy1)

	!velocidade normal relativa entre a partículas
	v_rel_normal = vrelx*normaldirectionx + vrely*normaldirectiony

	!interpenetração virtual
	csi = dabs(d - (a1 + a2))

	!raios corrigidos relativos ao ponto de contato
	a1_c = a1 - dabs(csi)/2.0d0
	a2_c = a2 - dabs(csi)/2.0d0


	!CAMPELLO - Tangential elasticity
	!calculate coefficients

	!effective Young's modulus of the contacting pair
	E_star = (E1*E2)/(E2*(1.0d0 - vpois1**2.0d0) + E1*(1.0d0 - vpois2**2.0d0))

	!effective shear modulus of the contacting pair
	G_star = (G1*G2)/(G1 + G2)
	
	!effective radius of the contacting pair	
	a_star = (a1*a2)/(a1 + a2)

	!effective mass of the contacting pair
	m_star = (m1*m2)/(m1 + m2)

	!stiffness of the spring
	k_at = 8.0d0*G_star*dsqrt(a_star)*dsqrt(abs(csi))
	
	!damping constant of dashpot
	d_at = 2.0d0*gama_s_dummy*dsqrt(m_star*k_at)


	!**********FORÇA ELÁSTICA*************
	!cálculo das componentes x e y da força elástica
	Fx_elastica_parcial = -(4.0d0/3.0d0)*dsqrt(a_star)*E_star*abs(csi)**(1.5d0)*(normaldirectionx)
  Fy_elastica_parcial = -(4.0d0/3.0d0)*dsqrt(a_star)*E_star*abs(csi)**(1.5d0)*(normaldirectiony)

	sinal_elasticax = Fx_elastica_parcial/abs(Fx_elastica_parcial)
	sinal_elasticax = Fy_elastica_parcial/abs(Fy_elastica_parcial)


	!**********FORÇA VISCOSA***************	
	!calculating contact damping constant, d_con (Campello's model for normal contacts)
	d_con = 2.0d0*gama_n_dummy*dsqrt(2.0d0*E_star*m_star*dsqrt(a_star))*(abs(csi)**0.25d0)

	!cálculo das componentes x e y da força viscosa
 	Fx_viscosa_parcial = (d_con)*(vrelx*normaldirectionx + vrely*normaldirectiony)*normaldirectionx
	Fy_viscosa_parcial = (d_con)*(vrelx*normaldirectionx + vrely*normaldirectiony)*normaldirectiony

	normal_aux = -(4.0d0/3.0d0)*dsqrt(a_star)*E_star*abs(csi)**(1.5d0) - (d_con*(vrelx*normaldirectionx + vrely*normaldirectiony))

	!artifício para que o amortecimento não supere a repulsão elástica
	if (normal_aux .le. 0.0d0) then

		!**********FORÇA NORMAL RESULTANTE******
		Fx_normal = Fx_elastica_parcial + Fx_viscosa_parcial
		Fy_normal = Fy_elastica_parcial + Fy_viscosa_parcial

		F_normal = dsqrt(Fx_normal**2.0d0 + Fy_normal**2.0d0)

		!sinais da força normal
		sinal_normalx = Fx_normal/abs(F_normal)
		sinal_normaly = Fy_normal/abs(F_normal)
	else
		Fx_normal = 0.0d0
		Fy_normal = 0.0d0
		
		F_normal = 0.0d0

		sinal_normalx = 0.0d0
		sinal_normaly = 0.0d0
	end if

	!velocidade relativa tangencial
	v_rel_t = vrelx*normaldirectiony - vrely*normaldirectionx + a1_c*o_1 + a2_c*o_2


!step 1: calculate trial friction force
	!trial incremental elongation is the sum of the elongation of the reference configuration projected on the tangential
	!direction of the contact plane and the trial incremental elongation of the pair
	dx_trial_x = (dx_acum_x*normaldirectiony) + (v_rel_t*normaldirectiony)*(dt)
	dx_trial_y = (dx_acum_y*(-normaldirectionx)) + (v_rel_t*(-normaldirectionx))*(dt)

	!trial friction force
	F_at_trial_x = -k_at*dx_trial_x + d_at*v_rel_t*(normaldirectiony)
	F_at_trial_y = -k_at*dx_trial_y + d_at*v_rel_t*(-normaldirectionx)

	F_at_trial = dsqrt(F_at_trial_x**2.0d0 + F_at_trial_y**2.0d0)

	!friction force
	Fs_x_dummy = F_at_trial_x
	Fs_y_dummy = F_at_trial_y

	Fs_dummy = F_at_trial

	!direction of the trial friction force
	t_trial_x = F_at_trial_x/abs(F_at_trial)
	t_trial_y = F_at_trial_y/abs(F_at_trial) 


!step 1: assume STICKING 
	if ((Fs_dummy) .le. (mi_dummy*F_normal)) then
		!sticking occurs

		!elongation of the spring is the trial elongation
		dx_acum_x = dx_trial_x
		dx_acum_y = dx_trial_y

		!direction of the tangential direction of friction force is the same direction as the trial friction force
		Fs_x_dummy = abs(F_at_trial_x)!*t_trial_x
		Fs_y_dummy = abs(F_at_trial_y)!*t_trial_y

		Fs_dummy = dsqrt(Fs_x_dummy**2.0d0 + Fs_y_dummy**2.0d0)

!step 2: check for SLIDING
	else if (Fs_dummy .gt. (mi_dummy*F_normal)) then
		!sliding occurs

		!friction force
		Fs_x_dummy = mi_dummy*F_normal*t_trial_x
		Fs_y_dummy = mi_dummy*F_normal*t_trial_y

		Fs_dummy = dsqrt(Fs_x_dummy**2.0d0 + Fs_y_dummy**2.0d0)

		!elongation of the spring
		!dx_acum_x = (1.0d0/k_at)*(mi_dummy*F_normal*normaldirectionx - d_at*v_rel_t*normaldirectionx)
		!dx_acum_y = (1.0d0/k_at)*(mi_dummy*F_normal*(-normaldirectiony) - d_at*v_rel_t*(-normaldirectiony))
		dx_acum_x = (mi_dummy/k_at)*F_normal*t_trial_x
		dx_acum_y = (mi_dummy/k_at)*F_normal*t_trial_y
	end if

	!componentes da velocidade relativa tangencial
	!vt_x = v_rel_t*normaldirectiony
	!vt_y = v_rel_t*(-normaldirectionx)

	!sinal da velocidade de rotação no contato oblíquo
	sinal_v_rot = v_rel_t/abs(v_rel_t)


	!**********TORQUE DE ATRITO DE ROLAMENTO*************
	!cálculo do torque da rolling friction
	if ((o_1-o_2) .EQ. 0.0d0) then	
		T_rolling_dummy = 0.0d0
	else
		T_rolling_dummy = -mi_roll_dummy*a_star*abs(F_normal)*(o_1 - o_2)/abs(o_1 - o_2)
	end if	

	!sinaliza que os cálculos das quantidades para o par de contato (part_p, part_q) já foram feitos
	!flag_dummy(part_p, part_q) = 1
	!flag_dummy(part_q, part_p) = 1

end subroutine


!soma as energias de cada partícula no dado instante
subroutine energias()
	use parametros
	implicit none
	integer :: l
	double precision, dimension(N) :: velocidade_partial !velocidade de cada partícula do sistema
	double precision, dimension(N) :: cinetica_partial, potg_partial, rot_partial !energias cinética, potencial gravitacional e rotacional de cada partícula do sistema

	velocidade_partial(:) = dsqrt((vxold(:)**2.0d0)+(vyold(:)**2.0d0))

	cinetica_partial(:) = (m(:)/2.0d0)*(velocidade_partial(:)**2.0d0)
	potg_partial(:) = m(:)*g*(yold(:))
	rot_partial(:) = (inertia(:)/2.0d0)*(omega_old(:)**2.0d0) 

	kinetic = sum(cinetica_partial)
	rotational_en = sum(rot_partial)
	potentialg = sum(potg_partial)

	totalenergy = kinetic + potentialg + rotational_en

	return
end subroutine
