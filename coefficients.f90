subroutine coeficientes()
	use parametros
	implicit none
	double precision :: p1, p2, p3
	double precision :: m_eff

	!massa efetiva de um conjunto de duas partículas
	m_eff = (m(1) * m(2))/(m(1) + m(2))

	!coeficiente normal de dissipação viscosa
	gama_n = 0.33
	
	!coeficiente de atrito dinâmico
	mi_t = 0.6d0

	!coeficiente de regularização do atrito dinâmico
	gama_s = 30000.0d0

	!coeficiente de atrito de rolamento (rolling friction)
	!os modelos referentes a estes coeficientes são baseados no artigo "Rolling friction in the dynamic simulation of sandpile formation", por Y.C. Zhou, B.D. Wright, R.Y. Yang, B.H. Xu, A.B. Yu 
	!(página 4, equações (3) e (4), respectivamente)
	!UPDATE: modelos baseados na página 83 de "Um modelo computacional para o estudo de materiais granulares" by Eduardo Campello
	mi_roll_1 = 0.1d0 !Modelo 1
	mi_roll_2 = 0.1d0 !1.0d0 !Modelo 2

	return
end subroutine coeficientes

