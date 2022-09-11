program Main
use setup_routines
use io_routines
use stiffness_routines
use post_process_routines

implicit none

type(global_params) :: globals
type(task_params) :: task
type(dof_matrices) :: dof_mats
type(mesh_info) :: mesh
type(post_process_info) :: post_proc
real, dimension(:,:), allocatable :: global_stiffness
real, dimension(:), allocatable :: inner_force, outer_force, global_force, global_reactions

! for sgetrf and sgetrs routines
integer, dimension(:), allocatable :: ipiv
integer :: gsd, info_int

call globals%set_globals(2,5,4)
call allocate_task_matrices(task, globals)
call allocate_dof_matrices(dof_mats, globals)

task%prop = reshape([1, 1, 1, 1, 1, &
                    200000, 200000, 200000, 200000, 200000, &
                    100, 200, 100, 200, 100, &
                    300000, 300000, 300000, 300000, 300000, &
                    100, 200, 100, 200, 100, &
                    100, 200, 100, 200, 100, &
                    100, 200, 100, 200, 100],shape(task%prop))

task%ien = reshape([1, 1, 2, 2, 3, &
                    2, 3, 3, 4, 4],shape(task%ien))

task%idb = reshape([1, 0, 0, 1, &
                    1, 0, 0, 1, &
                    0, 0, 0, 0, &
                    0, 0, 0, 0, &
                    0, 0, 0, 0, &
                    0, 0, 0, 0, &
                    0, 0, 0, 0], shape(task%idb))

task%xn = reshape([0.0, 4.e3, 4.e3, 8.e3, &
                   0.0, 0.0, 3.e3, 3.e3],shape(task%xn))
task%pu = 0
task%pu(3,2) = -9000.0

task%ds = 0
task%ds(1,1) = -4.0

call create_dof_matrices(dof_mats, task, globals)
call mesh%setup_table(task,globals)
call mesh%allocate_all_types(globals)
call mesh%setup_stiffness_mats(task, globals)

global_stiffness = setup_global_mat(mesh, dof_mats, task, globals)

call mesh%setup_forces(dof_mats, task, globals)

inner_force = setup_global_force(mesh, dof_mats, task, globals)
outer_force = mat_to_vec(task%pu,dof_mats%idu)
global_force = outer_force - inner_force
gsd = size(global_stiffness,dim=1)

allocate(ipiv(gsd))

call sgetrf(gsd,gsd, global_stiffness, gsd, ipiv, info_int)
call sgetrs('N', gsd, 1, global_stiffness, gsd, ipiv, global_force, gsd, info_int)
call post_proc%allocate_post_mats(mesh, task)

post_proc%global_displacements = vec_to_mat(global_force, dof_mats%idu) + task%ds
call post_proc%find_element_forces(mesh, dof_mats, task, globals)

global_reactions = find_reactions(post_proc, mesh, dof_mats, task, globals)
call write_mat_real(vec_to_mat(global_reactions,dof_mats%ids))

end program Main
