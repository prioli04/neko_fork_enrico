! Copyright (c) 2023-2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the `actuator_line_source_term_t` type.
module actuator_line_source_term
  use coefs, only : coef_t
  use device, only : device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use field, only : field_t
  use field_list, only : field_list_t
  use registry, only : neko_registry
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_lookup
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use source_term, only : source_term_t
  use utils, only : neko_error
  use time_state, only : time_state_t

  use global_interpolation, only: global_interpolation_t
  use matrix, only: matrix_t
  use vector, only: vector_t
  use math, only: pi
  use actuator_line_source_term_cpu, only : actuator_line_source_term_compute_cpu
  use actuator_line_source_term_device, only : actuator_line_source_term_compute_device
  implicit none
  private

  !> Actuator line id counter
  integer, save :: next_actuator_line_id = 0
  public :: make_registry_name

  !> Actuator line source term.
  type, public, extends(source_term_t) :: actuator_line_source_term_t
     private

     !> Actuator line id
     integer :: instance_id
     !> Number of actuator line elements
     integer :: N
     !> Wingspan
     real(kind=rp) :: b
     !> Gaussian kernel core size
     real(kind=rp) :: eps
     !> Actuator line center X coordinate
     real(kind=rp) :: x_center
     !> Actuator line center Y coordinate
     real(kind=rp) :: y_center
     !> Actuator line center Z coordinate
     real(kind=rp) :: z_center
     !> Actuator line elements' centroids
     type(matrix_t), pointer :: elem_center => null()
     !> Actuator line elements' nodal coordinates
     type(matrix_t), pointer :: nodes_coord => null()
     !> Gaussian kernel values (entire domain and for each actuator line elements)
     type(matrix_t), pointer :: kernel => null()
     !> Velocities interpolated at the actuator line (x-component)
     type(vector_t), pointer :: u_interp => null()
     !> Velocities interpolated at the actuator line (y-component)
     type(vector_t), pointer :: v_interp => null()
     !> Velocities interpolated at the actuator line (z-component)
     type(vector_t), pointer :: w_interp => null()
     !> Circulation vector
     real(kind=rp), allocatable :: gamma_vec(:)
     !> Sectional lift distribution
     real(kind=rp), allocatable :: lift_dist(:)
     !> Sectional drag distribution
     real(kind=rp), allocatable :: drag_dist(:)
     !> Resultant force vector
     type(vector_t), pointer :: resultant_force => null() 
     !> Interpolator object
     type(global_interpolation_t), pointer :: interpolator => null()
     !> Flag to indicate if the source term owns the interpolator or if it was provided externally
     logical :: owns_interpolator = .false.
     !> Aerodynamicc coefficients file id
     integer :: coefs_out_num 

   contains
     private

     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init => &
          actuator_line_source_term_init_from_json
     !> The constructor from type components.
     procedure, public, pass(this) :: init_from_components => &
          actuator_line_source_term_init_from_components
     !> Destructor.
     procedure, public, pass(this) :: free => actuator_line_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, public, pass(this) :: compute_ => actuator_line_source_term_compute

     ! ----------------------------------------------------------------------- !
     ! Internal subroutines for the actuator line method
     procedure, pass(this) :: discretize_line
     procedure, pass(this) :: compute_kernel
     procedure, pass(this) :: init_interpolator
     procedure, pass(this) :: F2D
     procedure, pass(this) :: compute_resultant_force
     procedure, nopass :: cross

  end type actuator_line_source_term_t

contains

  ! ========================================================================== !
  ! Public methods

  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine actuator_line_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(actuator_line_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    integer :: N
    real(kind=rp) :: b, eps, x_center, y_center, z_center
    real(kind=rp), allocatable :: gamma_vec(:)
    type(field_t), pointer :: u

    ! ------------------------------------------------------------------------ !
    ! Read the options for the actuator line source term
    call json_get(json, 'N', N)
    call json_get(json, 'b', b)
    call json_get(json, 'eps', eps)
    call json_get(json, 'xcenter', x_center)
    call json_get(json, 'ycenter', y_center)
    call json_get(json, 'zcenter', z_center)
    call json_get(json, 'gamma', gamma_vec)

    ! Initialize interpolator
    allocate(this%interpolator)
    this%owns_interpolator = .true.

    u => neko_registry%get_field_by_name('u') ! Use any field
    call this%interpolator%init(u%dof)

    call actuator_line_source_term_init_from_components(this, fields, coef, N, b, eps, &
        x_center, y_center, z_center, gamma_vec(:), this%interpolator)

  end subroutine actuator_line_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param N Number of actuator line elements.
  !! @param b Wingspan.
  !! @param eps Gaussian kernel core size.
  !! @param x_center Actuator line center X coordinate.
  !! @param y_center Actuator line center Y coordinate.
  !! @param z_center Actuator line center Z coordinate.
  !! @param gamma_vec Circulation vector.
  !! @param interp Global interpolator.
  !! @param id_out Optional output for the actuator line instance id.
  subroutine actuator_line_source_term_init_from_components(this, fields, coef, N, b, eps, &
        x_center, y_center, z_center, gamma_vec, interp, id_out)
    class(actuator_line_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    integer, intent(in) :: N
    real(kind=rp), intent(in) :: b
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: x_center
    real(kind=rp), intent(in) :: y_center
    real(kind=rp), intent(in) :: z_center
    real(kind=rp), intent(in) :: gamma_vec(:)
    type(global_interpolation_t), pointer, intent(in) :: interp
    integer, optional, intent(out) :: id_out

    real(kind=rp) :: start_time, end_time
    character(len=:), allocatable :: resultant_force_name

    ! Assign instance id
    next_actuator_line_id = next_actuator_line_id + 1
    this%instance_id = next_actuator_line_id

    ! Mandatory parameters for the general source term
    start_time = 0.0_rp
    end_time = huge(0.0_rp)

    ! Init base class
    call this%init_base(fields, coef, start_time, end_time)

    ! Check gamma size
    if (size(gamma_vec) .ne. N) then
       call neko_error('Gamma should have N elements.')
    end if

    ! Assign parameters
    this%N = N
    this%b = b
    this%eps = eps
    this%x_center = x_center
    this%y_center = y_center
    this%z_center = z_center
    this%gamma_vec = gamma_vec

    if (.not. this%owns_interpolator) then
       this%interpolator => interp
    end if

    ! Initialize variables
    call this%discretize_line() ! Create actuator line elements
    call this%compute_kernel() ! Pre-compute kernels
    call this%init_interpolator() ! Initialize interpolator

    resultant_force_name = make_registry_name("resultant_force", this%instance_id)
    call neko_registry%add_vector(2, resultant_force_name, .false.)
    this%resultant_force => neko_registry%get_vector(resultant_force_name)

    ! Create aerodynamics coefficients file
    this%coefs_out_num = 10
    open(unit=this%coefs_out_num, file="coefs.csv", status="replace", action="write")
    write(this%coefs_out_num, *) "t, Lift, Drag"
    
    if (present(id_out)) then
      id_out = this%instance_id
    end if

  end subroutine actuator_line_source_term_init_from_components

  !> Destructor.
  subroutine actuator_line_source_term_free(this)
    class(actuator_line_source_term_t), intent(inout) :: this

    if (allocated(this%gamma_vec)) then
      deallocate(this%gamma_vec)
    end if
    if (allocated(this%lift_dist)) then
      deallocate(this%lift_dist)
    end if
    if (allocated(this%drag_dist)) then
      deallocate(this%drag_dist)
    end if

    nullify(this%kernel)
    nullify(this%elem_center)
    nullify(this%nodes_coord)
    nullify(this%resultant_force)
    nullify(this%u_interp)
    nullify(this%v_interp)
    nullify(this%w_interp)

    if (associated(this%interpolator)) then
      if (this%owns_interpolator) then
        call this%interpolator%free()
        deallocate(this%interpolator)
      else
        nullify(this%interpolator)
      end if
    end if

    close(this%coefs_out_num)
    call this%free_base()

  end subroutine actuator_line_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine actuator_line_source_term_compute(this, time)
    class(actuator_line_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    
    integer :: j
    real(kind=rp) :: F_2D(this%N, 3)

    ! Update control points forces
    call this%F2D(time, F_2D)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call actuator_line_source_term_compute_device(this%kernel, F_2D, this%fields)
    else
       call actuator_line_source_term_compute_cpu(this%kernel, F_2D, this%fields)
    end if

  end subroutine actuator_line_source_term_compute

  ! ========================================================================== !
  ! Private methods

  !> Create actuator line elements
  subroutine discretize_line(this)
    class(actuator_line_source_term_t), intent(inout) :: this
    real(kind=rp) :: dy
    integer :: j
    character(len=:), allocatable :: elem_center_name, nodes_coord_name

    ! ------------------------------------------------------------------------ !
    ! Create actuator line elements
    ! 1st dimension: x, y, z coordinates
    ! 2nd dimension: elements
    elem_center_name = make_registry_name("elem_center", this%instance_id)
    call neko_registry%add_matrix(3, this%N, elem_center_name, .false.)
    this%elem_center => neko_registry%get_matrix(elem_center_name)

    nodes_coord_name = make_registry_name("nodes_coord", this%instance_id)
    call neko_registry%add_matrix(3, this%N + 1, nodes_coord_name, .false.)
    this%nodes_coord => neko_registry%get_matrix(nodes_coord_name)

    allocate(this%lift_dist(this%N))
    allocate(this%drag_dist(this%N))

    this%nodes_coord%x = 0.0_rp
    this%elem_center%x = 0.0_rp

    dy = this%b / real(this%N, kind=rp)
    this%nodes_coord%x(2, 1) = - 0.5_rp * this%b

    ! Uniformly spaced nodes
    do j = 1, this%N
      this%nodes_coord%x(2, j + 1) = this%nodes_coord%x(2, 1) + real(j, kind=rp) * dy
      this%elem_center%x(2, j) = 0.5_rp * (this%nodes_coord%x(2, j) + this%nodes_coord%x(2, j + 1))
    end do
    
  end subroutine discretize_line

  !> Computes Gaussian kernel, in the entire domain, for each actuator line element.
  subroutine compute_kernel(this)
    class(actuator_line_source_term_t), intent(inout) :: this

    integer :: i, j
    real(kind=rp) :: x_ij(3), xplus(3), xminus(3), t_hat(3)
    real(kind=rp) :: x, y, z, lminus, lplus, d
    type(field_t), pointer :: u
    character(len=:), allocatable :: kernel_name

    ! Allocate kernels
    u => neko_registry%get_field_by_name('u') ! Use any field
    kernel_name = make_registry_name("kernel", this%instance_id)
    call neko_registry%add_matrix(u%dof%size(), this%N, kernel_name, .false.)
    this%kernel => neko_registry%get_matrix(kernel_name)

    ! Loop through all actuator line elements
    do j = 1, this%N
      xplus = this%nodes_coord%x(:, j + 1) ! Right node coordinate
      xminus = this%nodes_coord%x(:, j) ! Left node coordinate

      t_hat = (xplus - xminus) / norm2(xplus - xminus) ! Actuator line element direction (from left node to right node)

      ! Loop through all mesh nodes
      do i = 1, u%dof%size()
        ! Get coordinates
        x = u%dof%x(i,1,1,1)
        y = u%dof%y(i,1,1,1)
        z = u%dof%z(i,1,1,1)

        x_ij = [x, y, z] ! Domain point to evaluate kernel

        lminus = dot_product(x_ij - xminus, t_hat) ! Length of projection of the distance between the domain point and the element's left node
        lplus = dot_product(x_ij - xplus, t_hat) ! Length of projection of the distance between the domain point and the element's right node
        d = norm2(cross(x_ij - xminus, t_hat)) ! Distance to the element

        this%kernel%x(i,j) = (1.0 / (pi * this%eps**2)) * exp(- d**2 / this%eps**2) * &
                           (erf(lminus / this%eps) - erf(lplus / this%eps)) / 2.0
      end do
    end do

  end subroutine compute_kernel

  !> Initialize interpolator and find points' owners
  subroutine init_interpolator(this)
    class(actuator_line_source_term_t), intent(inout) :: this
    character(len=:), allocatable :: u_name, v_name, w_name

    u_name = make_registry_name("u_interp", this%instance_id)
    v_name = make_registry_name("v_interp", this%instance_id)
    w_name = make_registry_name("w_interp", this%instance_id)

    call this%interpolator%find_points(this%elem_center%x, this%N) ! Find points and their owners

    ! Initialize interpolated velocities vectors
    call neko_registry%add_vector(this%N, u_name, .false.)
    call neko_registry%add_vector(this%N, v_name, .false.)
    call neko_registry%add_vector(this%N, w_name, .false.)

    this%u_interp => neko_registry%get_vector(u_name)
    this%v_interp => neko_registry%get_vector(v_name)
    this%w_interp => neko_registry%get_vector(w_name)

  end subroutine init_interpolator

  !> Updates control points forces.
  subroutine F2D(this, time, F_2D)
    class(actuator_line_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(out) :: F_2D(this%N, 3)

    integer :: j
    real(kind=rp) :: u_vec(3), lift, drag
    type (field_t), pointer :: up, vp, wp

    up => neko_registry%get_field_by_name("u")
    vp => neko_registry%get_field_by_name("v")
    wp => neko_registry%get_field_by_name("w")

    ! Retrieve updated field from GPU
    if (neko_bcknd_device .eq. 1) then
      call device_memcpy(up%x, up%x_d, up%size(), DEVICE_TO_HOST, sync=.true.)
      call device_memcpy(vp%x, vp%x_d, vp%size(), DEVICE_TO_HOST, sync=.true.)
      call device_memcpy(wp%x, wp%x_d, wp%size(), DEVICE_TO_HOST, sync=.true.)
    end if

    ! Interpolate velocities at ALM control points
    call this%interpolator%evaluate(this%u_interp%x, up%x, .false.)
    call this%interpolator%evaluate(this%v_interp%x, vp%x, .false.)
    call this%interpolator%evaluate(this%w_interp%x, wp%x, .false.)

    do j = 1, this%N
      u_vec = [this%u_interp%x(j), this%v_interp%x(j), this%w_interp%x(j)]
      F_2D(j,:) = cross(u_vec, [0.0_rp, this%gamma_vec(j), 0.0_rp])
      this%lift_dist(j) = F_2D(j,3)
      this%drag_dist(j) = F_2D(j,1)
    end do

    ! Print and save to coefs file
    call this%compute_resultant_force()
    lift = this%resultant_force%x(1)
    drag = this%resultant_force%x(2)
    write(this%coefs_out_num, *) time%t, ",", lift, ",", drag

  end subroutine F2D

  !> Compute 3D lift and drag coefficients 
  subroutine compute_resultant_force(this)
    class(actuator_line_source_term_t), intent(in) :: this

    this%resultant_force%x(1) = sum(this%lift_dist) * (this%b / real(this%N, kind=rp)) ! Lift
    this%resultant_force%x(2) = sum(this%drag_dist) * (this%b / real(this%N, kind=rp)) ! Drag

  end subroutine compute_resultant_force

  function make_registry_name(base_name, instance_id) result(registry_name)
    character(len=*), intent(in) :: base_name
    integer, intent(in) :: instance_id

    character(len=64) :: temp_name
    character(len=:), allocatable :: registry_name
 
    write(temp_name, '("alm_", A, "_", I0)') trim(base_name), instance_id
    registry_name = trim(temp_name)

  end function make_registry_name

  !> Computes the cross product of two 3D vectors.
  function cross(a, b) result(c)
    real(kind=rp), intent(in) :: a(3), b(3)
    real(kind=rp) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

end module actuator_line_source_term
