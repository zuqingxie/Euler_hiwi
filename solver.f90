subroutine solver
    implicit none

    call solver_kfvs
end subroutine solver

subroutine solver_kfvs
    use mod_solver_kfvs
    use mod_struct_to_array
    use mod_read_gmsh,     only: fname
    implicit none
    integer :: n

    call struct_to_array

    call donnee_initiale

    if (fname(1:3) == 'c31') then
        call allocate_vardummy_multi_elem_airfoil
        call conditions_aux_limites_multi_elem_airfoil
        call assign_lr_cell_multi_elem_airfoil
        call calcul_derived_quantities_multi_elem_airfoil
    else
        call allocate_vardummy
        call conditions_aux_limites
        call assign_lr_cell
        call calcul_derived_quantities
    endif

    call calcul_conservative_vector
        
    call timestep

    !--- temps de simulation et parametres de sorties
    tmax= 1.0d0 !0.3828823925d-2
    nmax=floor(tmax/dt)
    print *,'the maximum step is',nmax
    !--- evolution
    do n=1,nmax

        !-- calcul des flux
        if (fname(1:3) == 'c31') then
            call calcul_flux_multi_elem_airfoil
        else
            call calcul_flux
        endif

        !-- iteration en temps
        call calcul_rhs
        call euler_time_iteration

        !-- vérifier si convergé
        call chk_converge(n)

        !-- mise a jour
        vect_u=vect_unew

        !-- calcul de rho,ux,uy,t
        call calcul_rho_ux_uy_t
!            print *, rho(:,3)
        !-- mise a jour cl
        if (fname(1:3) == 'c31') then
            call conditions_aux_limites_multi_elem_airfoil
        else
            call conditions_aux_limites
        endif

        !-- mise à jour des quantités dérivées
        if (fname(1:3) == 'c31') then
            call calcul_derived_quantities_multi_elem_airfoil
        else
            call calcul_derived_quantities
        endif

        !-- sauvegarde resultats (format vtk, lisible par Paraview)
        if (mod(n,10000) == 0) then
            write(*,*) 'Writing solution file at iteration ', n, '...'
            !call write_solution_vtk(n)
            !call write_pressure_coefficient(n)
        endif
    enddo

end subroutine solver_kfvs

