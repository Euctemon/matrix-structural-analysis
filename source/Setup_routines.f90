module setup_routines
implicit none

! first collumn TRUSS element DOFS, second collumn BEAM element DOFS, third collumn FRAME element DOFS
! from top to bottom u_x, u_y, u_z, phi_x, phi_y, phi_z, temperature

integer, parameter, dimension(3,7) :: one_dim_dofs = reshape( &
    [1, 0, 1, &
     0, 1, 1, &
     0, 0, 0, &
     0, 0, 0, &
     0, 0, 0, &
     0, 1, 1, &
     0, 0, 0], shape(one_dim_dofs))
    
integer, parameter, dimension(3,7) :: two_dim_dofs = reshape( &
    [1, 0, 1, &
     1, 0, 1, &
     0, 1, 0, &
     0, 1, 0, &
     0, 1, 0, &
     0, 0, 1, &
     0, 0, 0], shape(two_dim_dofs))
    
integer, parameter, dimension(3,7) :: three_dim_dofs = reshape( &
    [1, 0, 1, &
     1, 0, 1, &
     1, 0, 1, &
     0, 0, 1, &
     0, 0, 1, &
     0, 0, 1, &
     0, 0, 0], shape(three_dim_dofs),order=[2,1])

type :: global_params
    integer :: ndf = 7                  ! number of degrees of freedom
    integer :: nsd                      ! number of spatial dimensions
    integer :: nel                      ! number of elements
    integer :: nnp                      ! number of nodal points
contains
    procedure :: set_globals
end type

type :: task_params
    integer, dimension(:,:), allocatable :: prop        ! matrix of element properties
    integer, dimension(:,:), allocatable :: ien         ! matrix of element node indices
    integer, dimension(:,:), allocatable :: idb         ! matrix of supported DOFS, 1 is essential (prescribed) BC, 0 is natural BC
    real, dimension(:,:), allocatable :: xn             ! matrix of xyz nodal coordinates
    real, dimension(:,:), allocatable :: pu             ! matrix of forces at unrestrained DOFS
    real, dimension(:,:), allocatable :: ds             ! matrix of prescribed displacements
contains
    procedure :: allocate_task_matrices
end type

type :: dof_matrices
    integer, dimension(:,:), allocatable :: ied         ! matrix of element DOFS, 1 means the element has that DOF, 0 not
    integer, dimension(:,:), allocatable :: idt         ! matrix of node DOFS, 1 means the node has that DOF, 0 not
    integer, dimension(:,:), allocatable :: ids         ! matrix of supported DOFS indices for nodes
    integer, dimension(:,:), allocatable :: idu         ! matrix of unrestrained DOFS indices for nodes
    integer, dimension(:), allocatable :: nen           ! CONVENIENCE VECTOR, number of element nodes for each element
    integer, dimension(:), allocatable :: ned           ! CONVENIENCE VECTOR, number of element DOFS
contains
    procedure :: allocate_dof_matrices
    procedure :: create_info_vecs
    procedure :: create_ied_mat
    procedure :: create_idt_mat
end type

contains
    pure subroutine set_globals(self,number_spatial_dims, number_elements, number_nodal_points)
        class(global_params), intent(inout) :: self
        integer, intent(in) :: number_spatial_dims, number_elements, number_nodal_points
        
        self%nsd = number_spatial_dims
        self%nel = number_elements
        self%nnp = number_nodal_points
    end subroutine set_globals
    
    pure subroutine allocate_task_matrices(self, globals)
        class(task_params), intent(inout) :: self
        class(global_params), intent(in) :: globals
        
        allocate(self%prop(globals%nel, 7))
        allocate(self%ien(globals%nel, 2))
        allocate(self%idb(globals%nnp, globals%ndf))
        allocate(self%xn(globals%nnp, globals%nsd))
        allocate(self%pu(globals%nnp, globals%ndf))
        allocate(self%ds(globals%nnp, globals%ndf))
    end subroutine allocate_task_matrices
    
    pure subroutine allocate_dof_matrices(self, globals)
        class(dof_matrices), intent(inout) :: self
        class(global_params), intent(in) :: globals
        
        allocate(self%ied(globals%nel, globals%ndf))
        allocate(self%idt(globals%nnp, globals%ndf))
        allocate(self%ids(globals%nnp, globals%ndf))
        allocate(self%idu(globals%nnp, globals%ndf))
        allocate(self%nen(globals%nel))
        allocate(self%ned(globals%nel))
    end subroutine allocate_dof_matrices

    pure function count_index(input_mat,every_row) result(new_mat)
        integer, dimension(:,:), intent(in) :: input_mat
        integer, intent(in) :: every_row                    ! if 1, then every row counting starts with 1, otherwise continues numbering indefinitely
        integer, dimension(:,:), allocatable :: new_mat
        integer :: count, i, j
    
        allocate(new_mat, mold=input_mat)
        count = 0
        new_mat = 0

        do i = 1,size(input_mat,dim=1)
            if (every_row .eq. 1) then
                count = 0
            end if
            do j = 1,size(input_mat,dim=2)
                if (input_mat(i,j) .gt. 0) then
                    count = count + 1
                    new_mat(i,j) = count
                end if
            end do
        end do
    end function count_index

    elemental function or_select(val_one, val_two) result(or_val)
        integer, intent(in) :: val_one, val_two
        integer :: or_val
        
        if ((val_one + val_two) .gt. 0) then
            or_val = 1
        else
            or_val = 0
        end if
        
    end function or_select

    elemental function and_select(val_one, val_two) result(and_val)
        integer, intent(in) :: val_one, val_two
        integer :: and_val
        
        if ((val_one*val_two) .gt. 0) then
            and_val = 1
        else
            and_val = 0
        end if
    end function and_select
    
    pure function cross(a,b) result(cross_vec)
        real, dimension(3), intent(in) :: a,b
        real, dimension(3) :: cross_vec
        
        cross_vec(1) = a(2) * b(3) - a(3) * b(2)
        cross_vec(2) = a(3) * b(1) - a(1) * b(3)
        cross_vec(3) = a(1) * b(2) - a(2) * b(1)
    end function cross

    pure subroutine create_info_vecs(self, task, globals)
        class(dof_matrices), intent(inout) :: self
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals
        integer :: i

        do i = 1,globals%nel
            self%nen(i) = count(task%ien(i,:) .gt.0)
            self%ned(i) = count(self%ied(i,:) .gt.0)
        end do
    end subroutine create_info_vecs
    
    pure subroutine create_ied_mat(self, task, globals)
        class(dof_matrices), intent(inout) :: self
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals
        integer, dimension(:,:), allocatable :: dim_elem_mat
        integer :: elem, elem_type
    
        select case (globals%nsd)
            case (1)
                dim_elem_mat = one_dim_dofs
            case (2)
                dim_elem_mat = two_dim_dofs
            case (3)
                dim_elem_mat = three_dim_dofs
        end select

        do elem=1,globals%nel
            elem_type = task%prop(elem,1)
            self%ied(elem,:) = dim_elem_mat(elem_type,:)
        end do

        self%ied = count_index(self%ied,1)
    end subroutine create_ied_mat

    pure subroutine create_idt_mat(self, task, globals)
        class(dof_matrices), intent(inout) :: self
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals
        integer :: elem, i, node
        
        self%idt = 0
        do elem = 1, globals%nel
            do i = 1,self%nen(elem)
                node = task%ien(elem,i)
                self%idt(node,:) = or_select(self%idt(node,:), self%ied(elem,:))
            end do
        end do
    end subroutine create_idt_mat
    
    pure subroutine create_dof_matrices(dof_mats, task, globals)
        class(dof_matrices), intent(inout) :: dof_mats
        class(task_params), intent(in) :: task
        class(global_params), intent(in) :: globals
        
        call dof_mats%create_ied_mat(task, globals)
        call dof_mats%create_info_vecs(task, globals)
        call dof_mats%create_ied_mat(task, globals)
        call dof_mats%create_idt_mat(task, globals)
        
        dof_mats%ids = count_index(task%idb,0)
        dof_mats%idu = count_index(dof_mats%idt - task%idb,0)
    end subroutine create_dof_matrices

    pure function mat_to_vec(info_mat, shape_mat) result(info_vec)
        real, dimension(:,:), intent(in) :: info_mat
        integer, dimension(:,:), intent(in) :: shape_mat
        real, dimension(:), allocatable :: info_vec
        integer :: neq, i, j
        
        neq = maxval(shape_mat)
        allocate(info_vec(neq))
        info_vec = 0
        
        do i=1,size(shape_mat,dim=1)
            do j=1,size(shape_mat,dim=2)
                if (shape_mat(i,j) .gt. 0) then
                    info_vec(shape_mat(i,j)) = info_mat(i,j)
                end if
            end do
        end do
    end function mat_to_vec

    pure function vec_to_mat(info_vec, shape_mat) result(info_mat)
        real, dimension(:), intent(in) :: info_vec
        integer, dimension(:,:), intent(in) :: shape_mat
        real, dimension(:,:), allocatable :: info_mat
        integer :: neq, i, j
        
        allocate(info_mat(size(shape_mat,dim=1),size(shape_mat,dim=2)))
        info_mat = 0
        
        do i=1,size(shape_mat,dim=1)
            do j=1,size(shape_mat,dim=2)
                if (shape_mat(i,j) .gt. 0) then
                    info_mat(i,j) = info_vec(shape_mat(i,j))
                end if
            end do
        end do
    end function vec_to_mat
end module setup_routines