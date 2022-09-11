module stiffness_routines
use setup_routines
use io_routines
implicit none

! Declaration of types

type :: truss_elems
    real, dimension(:,:,:), allocatable :: stiff_mats
    real, dimension(:,:,:), allocatable :: rot_mats
    real, dimension(:,:,:), allocatable :: rotated_stiff_mats
    real, dimension(:,:), allocatable :: local_displacements
    real, dimension(:,:), allocatable :: local_forces
end type

type :: beam_elems
    real, dimension(:,:,:), allocatable :: stiff_mats
    real, dimension(:,:,:), allocatable :: rot_mats
    real, dimension(:,:,:), allocatable :: rotated_stiff_mats
    real, dimension(:,:), allocatable :: local_displacements
    real, dimension(:,:), allocatable :: local_forces
end type

type :: frame_elems
    real, dimension(:,:,:), allocatable :: stiff_mats
    real, dimension(:,:,:), allocatable :: rot_mats
    real, dimension(:,:,:), allocatable :: rotated_stiff_mats
    real, dimension(:,:), allocatable :: local_displacements
    real, dimension(:,:), allocatable :: local_forces    
end type

type :: mesh_info
    integer, dimension(:,:), allocatable :: gtl_table       ! global to local indices table
    integer :: num_of_trusses = 0
    integer :: num_of_beams = 0
    integer :: num_of_frames = 0
    type(truss_elems) :: trusses
    type(beam_elems) :: beams
    type(frame_elems) :: frames
    
contains
    procedure :: setup_table
    procedure :: allocate_all_types
    procedure :: setup_stiffness_mats
    procedure :: write_stiffness_mats
    procedure :: setup_forces
end type

! Routines for elements

contains
    pure subroutine setup_table(self, task, globals)
        class(mesh_info), intent(inout) :: self
        class(global_params), intent(in) :: globals
        class(task_params), intent(in) :: task
        integer :: elem, elem_type

        allocate(self%gtl_table(2,globals%nel))
        do elem=1,globals%nel
            elem_type = task%prop(elem,1)
            self%gtl_table(1,elem)= elem_type
            select case (elem_type)
                case (1)
                    self%num_of_trusses = self%num_of_trusses + 1
                    self%gtl_table(2,elem) = self%num_of_trusses
                case (2)
                    self%num_of_beams = self%num_of_beams + 1
                    self%gtl_table(2,elem) = self%num_of_beams
                case (3)
                    self%num_of_frames = self%num_of_frames + 1
                    self%gtl_table(2,elem) = self%num_of_frames
            end select
        end do
    end subroutine setup_table

    pure subroutine allocate_truss_elems(self, globals)
        class(mesh_info), intent(inout) :: self
        class(global_params), intent(in) :: globals
        integer :: mat_size

        select case (globals%nsd)
            case (1)
                mat_size = 2
            case (2)
                mat_size = 4
            case (3)
                mat_size = 6
        end select
        allocate(self%trusses%stiff_mats(self%num_of_trusses,2,2))
        allocate(self%trusses%rot_mats(self%num_of_trusses,2,mat_size))
        allocate(self%trusses%rotated_stiff_mats(self%num_of_trusses,mat_size,mat_size))
        allocate(self%trusses%local_displacements(self%num_of_trusses,mat_size))
        allocate(self%trusses%local_forces(self%num_of_trusses,mat_size))

    end subroutine allocate_truss_elems
    
    pure subroutine allocate_beam_elems(self, globals)
        class(mesh_info), intent(inout) :: self
        class(global_params), intent(in) :: globals
        integer :: mat_size

        if (globals%nsd .eq. 1) then
            mat_size = 4
        else
            mat_size = 6
        end if
        allocate(self%beams%stiff_mats(self%num_of_beams, mat_size, mat_size))
        allocate(self%beams%rot_mats(self%num_of_beams, mat_size, mat_size))
        allocate(self%beams%rotated_stiff_mats(self%num_of_beams,mat_size,mat_size))
        allocate(self%beams%local_displacements(self%num_of_beams,mat_size))
        allocate(self%beams%local_forces(self%num_of_beams,mat_size))
        
    end subroutine allocate_beam_elems
    
    pure subroutine allocate_frame_elems(self, globals)
        class(mesh_info), intent(inout) :: self
        class(global_params), intent(in) :: globals
        integer :: mat_size

        if (globals%nsd .eq. 3) then
            mat_size = 12
        else
            mat_size = 6
        end if
        allocate(self%frames%stiff_mats(self%num_of_frames, mat_size, mat_size))
        allocate(self%frames%rot_mats(self%num_of_frames, mat_size, mat_size))
        allocate(self%frames%rotated_stiff_mats(self%num_of_frames,mat_size,mat_size))
        allocate(self%frames%local_displacements(self%num_of_frames,mat_size))
        allocate(self%frames%local_forces(self%num_of_frames,mat_size))

    end subroutine allocate_frame_elems
    
    pure subroutine allocate_all_types(self, globals)
        class(mesh_info), intent(inout) :: self
        class(global_params), intent(in) :: globals
        
        if (self%num_of_trusses .ne. 0) then
            call allocate_truss_elems(self,globals)
        end if
        
        if (self%num_of_beams .ne. 0) then
            call allocate_beam_elems(self,globals)
        end if
    
        if (self%num_of_frames .ne. 0) then
            call allocate_frame_elems(self,globals)
        end if
    end subroutine allocate_all_types

    pure subroutine compute_truss_stiffness(stiff_mat, rot_mat, xn_mat, prop_vec, ien_vec)
        real, dimension(:,:), intent(out) :: stiff_mat, rot_mat
        integer, dimension(:), intent(in) :: prop_vec, ien_vec
        real, dimension(:,:), intent(in) :: xn_mat
        
        real, dimension(:), allocatable :: nx       ! orientation vector of the frame element
        integer :: E, A, first_node, second_node, nx_size
        real :: L
    
        stiff_mat = reshape([1, -1, -1, 1], shape(stiff_mat))
        rot_mat = 0

        E = prop_vec(2)
        A = prop_vec(3)
        
        first_node = ien_vec(1)
        second_node = ien_vec(2)
        nx = xn_mat(second_node,:)-xn_mat(first_node,:)
        L = norm2(nx)
        nx = nx/L
        
        stiff_mat = (E*A/L)*stiff_mat
        nx_size = size(nx)
        rot_mat(1,1:nx_size)=nx
        rot_mat(2,nx_size+1:2*nx_size)=nx
    end subroutine compute_truss_stiffness
    
    pure subroutine compute_beam_stiffness(stiff_mat, rot_mat, xn_mat, prop_vec, ien_vec)
        real, dimension(:,:), intent(out) :: stiff_mat, rot_mat
        integer, dimension(:), intent(in) :: prop_vec, ien_vec        
        real, dimension(:,:), intent(in) :: xn_mat
        
        real, dimension(3,3) :: qe
        real, dimension(:), allocatable :: nx       ! orientation vector of the truss element
        integer :: E, Ix, J, v, first_node, second_node
        real :: L, G, k1, k2, k3, k4
    

        E = prop_vec(2)
        Ix = prop_vec(4)
        J = prop_vec(5)
        v = prop_vec(6)
        G = E/(2+2*v)
        
        first_node = ien_vec(1)
        second_node = ien_vec(2)
        nx = xn_mat(second_node,:)-xn_mat(first_node,:)
        L = norm2(nx)
        nx = nx/L
        
        k1 = 12*E*Ix/L**3
        k2 = 6*E*Ix/L**2
        k3 = 4*E*Ix/L
        k4 = G*J/L

        if (size(nx) .eq. 1) then
            stiff_mat = reshape([k1, k2, -k1, k2, &
                                    k2, k3, -k2, k3/2, &
                                    -k1, -k2, k1, -k2, &
                                    k2, k3/2, -k2, k3], shape(stiff_mat))
            rot_mat = 0.0
            rot_mat(1,1) = nx(1)
            rot_mat(2,2) = 1.0
            rot_mat(3,3) = nx(1)
            rot_mat(4,4) = 1.0
        else
            stiff_mat = reshape([k1, 0.0, -k2, -k1, 0.0, -k2, &
                                0.0, k4, 0.0, 0.0, -k4, 0.0, &
                                -k2, 0.0, k3, k2, 0.0, k3/2, &
                                -k1, 0.0, k2, k1, 0.0, k2, &
                                0.0, -k4, 0.0, 0.0, k4, 0.0, &
                                -k2, 0.0, k3/2, k2, 0.0, k3], shape(stiff_mat))
            rot_mat = 0.0
            qe = reshape([1.0, 0.0, 0.0, &
                            0.0, nx(1), nx(2), &
                            0.0, -nx(2), nx(1)], shape(qe))
            rot_mat(1:3,1:3) = qe
            rot_mat(4:6, 4:6) = qe
        end if
    end subroutine compute_beam_stiffness
    
    pure subroutine compute_frame_stiffness(stiff_mat, rot_mat, xn_mat, prop_vec, ien_vec)
        real, dimension(:,:), intent(out) :: stiff_mat, rot_mat
        integer, dimension(:), intent(in) :: prop_vec, ien_vec
        real, dimension(:,:), intent(in) :: xn_mat
        
        real, dimension(3,3) :: qe
        real, dimension(:), allocatable :: nx,ny       ! orientation vector of the truss element
        integer :: E, A, Iz, Iy, J, v, first_node, second_node
        real :: L, G, k1, k2, k3, k4, k5, k6, k7, k8
    
    
        E = prop_vec(1)
        A = prop_vec(2)
        Iz = prop_vec(3)
        Iy = prop_vec(4)
        J = prop_vec(5)
        v = prop_vec(6)
        G = E/(2+2*v)
        
        first_node = ien_vec(1)
        second_node = ien_vec(2)
        nx = xn_mat(second_node,:)-xn_mat(first_node,:)
        L = norm2(nx)
        nx = nx/L
        
        k1 = E*A/L
        k2 = 12*E*Iz/l**3
        k3 = 6*E*Iz/L**2
        k4 = 4*E*Iz/L
        k5 = 12*E*Iy/L**3
        k6 = 6*E*Iy/L**2
        k7 = 4*E*Iy/L
        k8 = G*J/L
        
        if (size(nx) .eq. 3) then
            stiff_mat = reshape([k1, 0.0, 0.0, 0.0, 0.0, 0.0, -k1, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                 0.0, k2, 0.0, 0.0, 0.0, k3, 0.0, -k2, 0.0, 0.0, 0.0, k3, &
                                 0.0, 0.0, k5, 0.0, -k6, 0.0, 0.0, 0.0, -k5, 0.0, -k6, 0.0, &
                                 0.0, 0.0, 0.0, k8, 0.0, 0.0, 0.0, 0.0, 0.0, -k8, 0.0, 0.0, &
                                 0.0, 0.0, -k6, 0.0, k7, 0.0, 0.0, 0.0, k6, 0.0, k7/2, 0.0, &
                                 0.0, k3, 0.0, 0.0, 0.0, k4, 0.0, -k3, 0.0, 0.0, 0.0, k4/2, &
                                 -k1, 0.0, 0.0, 0.0, 0.0, 0.0, k1, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                 0.0, -k2, 0.0, 0.0, 0.0, -k3, 0.0, k2, 0.0, 0.0, 0.0, -k3, &
                                 0.0, 0.0, -k5, 0.0, k6, 0.0, 0.0, 0.0, k5, 0.0, k6, 0.0, &
                                 0.0, 0.0, 0.0, -k8, 0.0, 0.0, 0.0, 0.0, 0.0, k8, 0.0, 0.0, &
                                 0.0, 0.0, -k6, 0.0, k7/2, 0.0, 0.0, 0.0, k6, 0.0, k7, 0.0, &
                                 0.0, k3, 0.0, 0.0, 0.0, k4/2, 0.0, -k3, 0.0, 0.0, 0.0, k4], shape(stiff_mat))
            if (nx(3) .eq. 1.0) then
                ny = [0.0, 1.0, 0.0]
            else
                ny = [-nx(2), nx(1), 0.0]
                ny = ny/norm2(nx(1:2))
            end if
            qe(1,:) = nx
            qe(2,:) = ny
            qe(3,:) = cross(nx,ny)
            rot_mat(1:3,1:3) = qe
            rot_mat(4:6,4:6) = qe
            rot_mat(7:9,7:9) = qe
            rot_mat(10:12,10:12) = qe
        else
            stiff_mat = reshape([k1, 0.0, 0.0, -k1, 0.0, 0.0, &
                                 0.0, k2, k3, 0.0, -k2, k3, &
                                 0.0, k3, k4, 0.0, -k3, k4/2, &
                                 -k1, 0.0, 0.0, k1, 0.0, 0.0, &
                                 0.0, -k2, -k3, 0.0, k2, -k3, &
                                 0.0, k3, k4/2, 0.0, -k3, k4], shape(stiff_mat))
            if (size(nx) .eq. 1) then
                qe = reshape([nx, 0.0, 0.0, &
                              0.0, nx, 0.0, &
                              0.0, 0.0, 1.0], shape(qe))
            else
                qe = reshape([nx(1), nx(2), 0.0, &
                              -nx(2), nx(1), 0.0, &
                              0.0, 0.0, 1.0], shape(qe))
            end if
            rot_mat(1:3,1:3) = qe
            rot_mat(4:6,4:6) = qe
        end if
    end subroutine compute_frame_stiffness

    pure subroutine setup_stiffness_mats(self,task,globals)
        class(mesh_info), intent(inout) :: self
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals
        real, dimension(:,:), allocatable :: right_mul
        integer :: elem, local_index
    
        do elem=1,globals%nel
            local_index = self%gtl_table(2,elem)
            select case (self%gtl_table(1,elem))
            case (1)
                call compute_truss_stiffness(self%trusses%stiff_mats(local_index,:,:),self%trusses%rot_mats(local_index,:,:), task%xn, task%prop(elem,:), task%ien(elem,:))
                right_mul = matmul(self%trusses%stiff_mats(local_index,:,:),self%trusses%rot_mats(local_index,:,:))
                self%trusses%rotated_stiff_mats(local_index,:,:) = matmul(transpose(self%trusses%rot_mats(local_index,:,:)),right_mul)
            case (2)
                call compute_beam_stiffness(self%beams%stiff_mats(local_index,:,:),self%beams%rot_mats(local_index,:,:), task%xn, task%prop(elem,:), task%ien(elem,:))
                right_mul = matmul(self%beams%stiff_mats(local_index,:,:),self%beams%rot_mats(local_index,:,:))
                self%beams%rotated_stiff_mats(local_index,:,:) = matmul(transpose(self%beams%rot_mats(local_index,:,:)),right_mul)
            case (3)
                call compute_frame_stiffness(self%frames%stiff_mats(local_index,:,:),self%frames%rot_mats(local_index,:,:), task%xn, task%prop(elem,:), task%ien(elem,:))
                right_mul = matmul(self%frames%stiff_mats(local_index,:,:),self%frames%rot_mats(local_index,:,:))
                self%frames%rotated_stiff_mats(local_index,:,:) = matmul(transpose(self%frames%rot_mats(local_index,:,:)),right_mul)
            end select
        end do
    end subroutine setup_stiffness_mats
    
    subroutine write_stiffness_mats(self, globals)
        class(mesh_info), intent(inout) :: self
        class(global_params), intent(in) :: globals
        integer :: elem, local_index
    
        do elem=1,globals%nel
            local_index = self%gtl_table(2,elem)
            select case (self%gtl_table(1,elem))
            case (1)
                print *,'truss'
                call write_mat_real(self%trusses%rotated_stiff_mats(local_index,:,:))
                print *,''
            case (2)
                print *,'beam'
                call write_mat_real(self%beams%rotated_stiff_mats(local_index,:,:))
                print *,''
            case (3)
                print *,'frame'
                call write_mat_real(self%frames%rotated_stiff_mats(local_index,:,:))
                print *,''
            end select
            
        end do
    end subroutine write_stiffness_mats
    
    ! Routines for global stiffness and forces
    
    function setup_global_mat(self, dof_mats, task, globals) result(global_mat)
        class(mesh_info), intent(in), target :: self
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals    
        class(dof_matrices), intent(in) :: dof_mats
        
        real, dimension(:,:), allocatable :: global_mat
        real, dimension(:,:), pointer :: current_stiff_mat
        integer :: i, j, elem, global_p, global_q, local_p, local_q, m, n, local_index, neq

        neq = maxval(dof_mats%idu)        
        allocate(global_mat(neq,neq))
        global_mat = 0

        do elem = 1, globals%nel
            local_index = self%gtl_table(2,elem)
            select case (self%gtl_table(1,elem))
                case (1)
                    current_stiff_mat => self%trusses%rotated_stiff_mats(local_index,:,:)
                case (2)
                    current_stiff_mat => self%beams%rotated_stiff_mats(local_index,:,:)
                case (3)
                    current_stiff_mat => self%frames%rotated_stiff_mats(local_index,:,:)
            end select

            do m = 1,dof_mats%nen(elem)
                do i = 1,globals%ndf
                    global_p = dof_mats%idu(task%ien(elem,m),i)
                    if (and_select(global_p,dof_mats%ied(elem,i)) .eq. 1) then
                        do n=1,dof_mats%nen(elem)
                            do j=1,globals%ndf
                                global_q = dof_mats%idu(task%ien(elem,n),j)
                                if (and_select(global_q,dof_mats%ied(elem,j)) .eq. 1) then
                                    local_p = (m-1)*dof_mats%ned(elem) + dof_mats%ied(elem,i)
                                    local_q = (n-1)*dof_mats%ned(elem) + dof_mats%ied(elem,j)
                                    global_mat(global_p,global_q) = global_mat(global_p,global_q) + current_stiff_mat(local_p,local_q)
                                end if
                            end do
                        end do
                    end if
                end do
            end do
        end do
    end function setup_global_mat
    
    pure subroutine setup_forces(self, dof_mats, task, globals)
        class(mesh_info), intent(inout) :: self
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals    
        class(dof_matrices), intent(in) :: dof_mats

        integer :: elem, i, j, local_index, p, n
    
        do elem=1,globals%nel
            local_index = self%gtl_table(2,elem)
            do n=1,dof_mats%nen(elem)
                do i=1,globals%ndf
                    if (dof_mats%ied(elem,i) .gt. 0) then
                        p = (n-1)*dof_mats%ned(elem) + dof_mats%ied(elem,i)
                        select case (self%gtl_table(1,elem))
                            case (1)
                                self%trusses%local_displacements(local_index,p) = task%ds(task%ien(elem,n),i)
                            case (2)
                                self%beams%local_displacements(local_index,p) = task%ds(task%ien(elem,n),i)
                            case (3)
                                self%frames%local_displacements(local_index,p) = task%ds(task%ien(elem,n),i)
                        end select
                    end if
                end do
            end do
            select case (self%gtl_table(1,elem))
                case (1)
                    self%trusses%local_forces(local_index,:)=matmul(self%trusses%rotated_stiff_mats(local_index,:,:),self%trusses%local_displacements(local_index,:))
                case (2)
                    self%beams%local_forces(local_index,:)=matmul(self%beams%rotated_stiff_mats(local_index,:,:),self%beams%local_displacements(local_index,:))
                case (3)
                    self%frames%local_forces(local_index,:)=matmul(self%frames%rotated_stiff_mats(local_index,:,:),self%frames%local_displacements(local_index,:))
            end select
        end do
    end subroutine setup_forces
    
    function setup_global_force(self, dof_mats, task, globals) result(global_force)
        class(mesh_info), intent(in), target :: self
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals
        class(dof_matrices), intent(in) :: dof_mats

        integer :: elem, i, n, local_index, p, local_p
        real, dimension(:), pointer :: current_force_vec
        real, dimension(:), allocatable :: global_force

        allocate(global_force(maxval(dof_mats%idu)))
        global_force = 0
    
        do elem=1, globals%nel
            local_index = self%gtl_table(2,elem)
            select case (self%gtl_table(1,elem))
                case (1)
                    current_force_vec => self%trusses%local_forces(local_index,:)
                case (2)
                    current_force_vec => self%beams%local_forces(local_index,:)
                case (3)
                    current_force_vec => self%frames%local_forces(local_index,:)
            end select

            do n=1,dof_mats%nen(elem)
                do i=1,globals%ndf
                    p = dof_mats%idu(task%ien(elem,n),i)
                    if (and_select(P,dof_mats%ied(elem,i)) .eq. 1) then
                        local_p = (n-1)*dof_mats%ned(elem) + dof_mats%ied(elem,i)
                        global_force(p) = global_force(p) + current_force_vec(local_p)
                    end if
                end do
            end do
        end do
    end function setup_global_force

end module stiffness_routines