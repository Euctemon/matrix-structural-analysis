module post_process_routines
use stiffness_routines
use setup_routines

implicit none
    
type :: post_process_info
    real, dimension(:,:), allocatable :: global_displacements
    real, dimension(:,:), allocatable :: truss_forces
    real, dimension(:,:), allocatable :: truss_displacements
    real, dimension(:,:), allocatable :: beam_forces
    real, dimension(:,:), allocatable :: beam_displacements
    real, dimension(:,:), allocatable :: frame_forces
    real, dimension(:,:), allocatable :: frame_displacements
contains
    procedure :: allocate_post_mats
    procedure :: find_element_forces
    procedure :: find_reactions
end type

contains
    pure subroutine allocate_post_mats(self, mesh, task)
        class(post_process_info), intent(inout) :: self
        class(mesh_info), intent(in) :: mesh
        class(task_params), intent(in) :: task

        allocate(self%global_displacements, mold = task%ds)
        
        if (mesh%num_of_trusses .gt. 0) then
            allocate(self%truss_forces, mold = mesh%trusses%local_forces)
            allocate(self%truss_displacements, mold = mesh%trusses%local_displacements)
        end if
        
        if (mesh%num_of_beams .gt. 0) then
            allocate(self%beam_forces, mold = mesh%beams%local_forces)
            allocate(self%beam_displacements, mold = mesh%beams%local_displacements)
        end if
        
        if (mesh%num_of_frames .gt. 0) then
            allocate(self%frame_forces, mold = mesh%frames%local_forces)
            allocate(self%frame_displacements, mold = mesh%frames%local_displacements)
        end if
    end subroutine allocate_post_mats

    pure subroutine find_element_forces(self, mesh, dof_mats, task, globals)
        class(post_process_info), intent(inout) :: self
        class(mesh_info), intent(in) :: mesh
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals    
        class(dof_matrices), intent(in) :: dof_mats

        integer :: elem, i, j, local_index, p, n
    
        do elem=1,globals%nel
            local_index = mesh%gtl_table(2,elem)
            do n=1,dof_mats%nen(elem)
                do i=1,globals%ndf
                    if (dof_mats%ied(elem,i) .gt. 0) then
                        p = (n-1)*dof_mats%ned(elem) + dof_mats%ied(elem,i)
                        select case (mesh%gtl_table(1,elem))
                            case (1)
                                self%truss_displacements(local_index,p) = self%global_displacements(task%ien(elem,n),i)
                            case (2)
                                self%beam_displacements(local_index,p) = self%global_displacements(task%ien(elem,n),i)
                            case (3)
                                self%frame_displacements(local_index,p) = self%global_displacements(task%ien(elem,n),i)
                        end select
                    end if
                end do
            end do
            select case (mesh%gtl_table(1,elem))
                case (1)
                    self%truss_forces(local_index,:)=matmul(mesh%trusses%rotated_stiff_mats(local_index,:,:),self%truss_displacements(local_index,:))
                case (2)
                    self%beam_forces(local_index,:)=matmul(mesh%beams%rotated_stiff_mats(local_index,:,:),self%beam_displacements(local_index,:))
                case (3)
                    self%frame_forces(local_index,:)=matmul(mesh%frames%rotated_stiff_mats(local_index,:,:),self%frame_displacements(local_index,:))
            end select
        end do
    end subroutine find_element_forces

    function find_reactions(self, mesh, dof_mats, task, globals) result(global_force)
        class(post_process_info), intent(in), target :: self
        class(mesh_info), intent(in) :: mesh
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals
        class(dof_matrices), intent(in) :: dof_mats

        integer :: elem, i, n, local_index, p, local_p
        real, dimension(:), pointer :: current_force_vec
        real, dimension(:), allocatable :: global_force

        allocate(global_force(maxval(dof_mats%ids)))
        global_force = 0
    
        do elem=1, globals%nel
            local_index = mesh%gtl_table(2,elem)
            select case (mesh%gtl_table(1,elem))
                case (1)
                    current_force_vec => self%truss_forces(local_index,:)
                case (2)
                    current_force_vec => self%beam_forces(local_index,:)
                case (3)
                    current_force_vec => self%frame_forces(local_index,:)
            end select

            do n=1,dof_mats%nen(elem)
                do i=1,globals%ndf
                    p = dof_mats%ids(task%ien(elem,n),i)
                    if (and_select(p,dof_mats%ied(elem,i)) .eq. 1) then
                        local_p = (n-1)*dof_mats%ned(elem) + dof_mats%ied(elem,i)
                        global_force(p) = global_force(p) + current_force_vec(local_p)
                    end if
                end do
            end do
        end do
    end function find_reactions
end module post_process_routines