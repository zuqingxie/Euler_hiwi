
subroutine pre_processing
    use mod_objects
    implicit none

    real :: startTime,endTime

    call cpu_time(startTime)

    call list_cell%read_msh()
    call list_cell%assign_id_nodes()

    call list_cell%build_cells()
    call list_cell%find_neighbor()

    visualizer = obj_msh_vis(list_cell)                           ! replace a pointer assignment with a Fortran constructor.
    call visualizer%write_vtk()

    call cpu_time(endTime)
    write(*,*) 'Time for treating mesh:', endTime - startTime
    
    call visualizer%write_tecplot()

    call list_cell%calcul_vol()
    call list_cell%calcul_face_info()
    call list_cell%assign_face_id()

    faces_fvm = obj_fvm_face_2D(list_cell)                        ! replace a pointer assignment with a Fortran constructor.
    call faces_fvm%alloc()
    call faces_fvm%get()

end subroutine pre_processing

