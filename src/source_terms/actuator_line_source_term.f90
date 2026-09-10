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

  !> Actuator line source term.
  type, public, extends(source_term_t) :: actuator_line_source_term_t
     private

     !> Actuator line id
     integer :: instance_id
     !> Number of actuator line elements
     integer :: N
     !> Target lift coefficient
     real(kind=rp) :: CL
     !> Wingspan
     real(kind=rp) :: b
     !> Aspect ratio
     real(kind=rp) :: AR
     !> Gaussian kernel core size
     real(kind=rp) :: eps
     !> Actuator line center X coordinate
     real(kind=rp) :: x_center
     !> Actuator line center Y coordinate
     real(kind=rp) :: y_center
     !> Actuator line center Z coordinate
     real(kind=rp) :: z_center
     !> Actuator line elements' midpoint coordinates
     type(matrix_t), pointer :: x_vec => null()
     !> Actuator line elements' right node coordinates
     real(kind=rp), allocatable :: xplus_vec(:,:)
     !> Actuator line elements' left node coordinates
     real(kind=rp), allocatable :: xminus_vec(:,:)
     !> Gaussian kernel values (entire domain and for each actuator line elements)
     real(kind=rp), allocatable :: kernel(:,:)
     !> Field for holding kernel values and performing operations on the GPU
     type(field_t), pointer :: kernel_j => null()
     !> Velocities interpolated at the actuator line (x-component)
     type(vector_t), pointer :: u_interp => null()
     !> Velocities interpolated at the actuator line (y-component)
     type(vector_t), pointer :: v_interp => null()
     !> Velocities interpolated at the actuator line (z-component)
     type(vector_t), pointer :: w_interp => null()
     !> Circulation vector
     real(kind=rp), allocatable :: gamm(:)
     !> Sectional lift force
     real(kind=rp), allocatable :: lift(:)
     !> Sectional drag force
     real(kind=rp), allocatable :: drag(:)
     !> Force coefficients vector [CL, CD]
     type(vector_t), pointer :: coefs => null()
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
     procedure, public, pass(this) :: init_from_compenents => &
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
     procedure, pass(this) :: compute_coefficients
     procedure, pass(this) :: make_registry_name
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

    integer :: N_in
    real(kind=rp) :: CL_in, b_in, AR_in, eps_in, x_center_in, y_center_in, z_center_in
    real(kind=rp), allocatable :: gamma_in(:)
    type(field_t), pointer :: u

    ! ------------------------------------------------------------------------ !
    ! Read the options for the actuator line source term
    call json_get(json, 'N', N_in)
    call json_get(json, 'CL_target', CL_in)
    call json_get(json, 'b', b_in)
    call json_get(json, 'AR', AR_in)
    call json_get(json, 'eps', eps_in)
    call json_get(json, 'xcenter', x_center_in)
    call json_get(json, 'ycenter', y_center_in)
    call json_get(json, 'zcenter', z_center_in)
    call json_get(json, 'gamma', gamma_in)

    ! Initialize interpolator
    allocate(this%interpolator)
    this%owns_interpolator = .true.

    u => neko_registry%get_field_by_name('u') ! Use any field
    call this%interpolator%init(u%dof)

    call actuator_line_source_term_init_from_components(this, fields, coef, N_in, CL_in, b_in, AR_in, eps_in, &
        x_center_in, y_center_in, z_center_in, gamma_in(:), this%interpolator)

  end subroutine actuator_line_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param N_in Number of actuator line elements.
  !! @param CL_in Target lift coefficient.
  !! @param b_in Wingspan.
  !! @param AR_in Aspect ratio.
  !! @param eps_in Gaussian kernel core size.
  !! @param x_center_in Actuator line center X coordinate.
  !! @param y_center_in Actuator line center Y coordinate.
  !! @param z_center_in Actuator line center Z coordinate.
  !! @param gamma_in Circulation vector.
  !! @param interp Global interpolator.
  !! @param id_out Optional output for the actuator line instance id.
  subroutine actuator_line_source_term_init_from_components(this, fields, coef, N_in, CL_in, b_in, AR_in, eps_in, &
        x_center_in, y_center_in, z_center_in, gamma_in, interp, id_out)
    class(actuator_line_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    integer, intent(in) :: N_in
    real(kind=rp), intent(in) :: CL_in
    real(kind=rp), intent(in) :: b_in
    real(kind=rp), intent(in) :: AR_in
    real(kind=rp), intent(in) :: eps_in
    real(kind=rp), intent(in) :: x_center_in
    real(kind=rp), intent(in) :: y_center_in
    real(kind=rp), intent(in) :: z_center_in
    real(kind=rp), intent(in) :: gamma_in(:)
    type(global_interpolation_t), pointer, intent(in) :: interp
    integer, optional, intent(out) :: id_out

    real(kind=rp) :: start_time, end_time
    character(len=:), allocatable :: coefs_name

    ! Assign instance id
    next_actuator_line_id = next_actuator_line_id + 1
    this%instance_id = next_actuator_line_id

    ! Mandatory parameters for the general source term
    start_time = 0.0_rp
    end_time = huge(0.0_rp)

    ! Init base class
    call this%init_base(fields, coef, start_time, end_time)

    ! Check gamma size
    if (size(gamma_in) .ne. N_in) then
       call neko_error('Gamma should have N elements.')
    end if

    ! Assign parameters
    this%N = N_in
    this%CL = CL_in
    this%b = b_in
    this%AR = AR_in
    this%eps = eps_in
    this%x_center = x_center_in
    this%y_center = y_center_in
    this%z_center = z_center_in
    this%gamm = gamma_in

    if (.not. this%owns_interpolator) then
       this%interpolator => interp
    end if

    ! Initialize variables
    call this%discretize_line() ! Create actuator line elements
    call this%compute_kernel() ! Pre-compute kernels
    call this%init_interpolator() ! Initialize interpolator

    call this%make_registry_name("coefs", coefs_name)
    call neko_registry%add_vector(2, coefs_name, .false.)
    this%coefs => neko_registry%get_vector(coefs_name)

    ! Create aerodynamics coefficients file
    this%coefs_out_num = 10
    open(unit=this%coefs_out_num, file="coefs.csv", status="replace", action="write")
    write(this%coefs_out_num, *) "t, CL, CD"
    
    ! gamma_0 = 2 * CL * V * b / (pi * AR)
    ! Gamma(y) = gamma_0 * sqrt(1 - (2y/b)^2)
    ! gamm(j) = 2.0 * alm_CL * V_inf * alm_b / (pi * alm_AR) * sqrt(1 - (2.0 * alm_x_vec(j, 2) / alm_b)**2)

    if (present(id_out)) then
      id_out = this%instance_id
    end if

  end subroutine actuator_line_source_term_init_from_components

  !> Destructor.
  subroutine actuator_line_source_term_free(this)
    class(actuator_line_source_term_t), intent(inout) :: this

    if (allocated(this%xplus_vec)) then
      deallocate(this%xplus_vec)
    end if
    if (allocated(this%xminus_vec)) then
      deallocate(this%xminus_vec)
    end if
    if (allocated(this%kernel)) then
      deallocate(this%kernel)
    end if
    if (allocated(this%lift)) then
      deallocate(this%lift)
    end if
    if (allocated(this%drag)) then
      deallocate(this%drag)
    end if
    if (allocated(this%gamm)) then
      deallocate(this%gamm)
    end if

    if (associated(this%kernel_j)) nullify(this%kernel_j)
    nullify(this%x_vec)
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
       call actuator_line_source_term_compute_device(this%kernel, this%kernel_j, F_2D, this%fields)
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
    character(len=:), allocatable :: x_vec_name

    ! ------------------------------------------------------------------------ !
    ! Create actuator line elements
    ! 1st dimension: x, y, z coordinates
    ! 2nd dimension: elements
    call this%make_registry_name("x_vec", x_vec_name)
    call neko_registry%add_matrix(3, this%N, x_vec_name, .false.)
    this%x_vec => neko_registry%get_matrix(x_vec_name)

    allocate(this%xplus_vec(3, this%N))
    allocate(this%xminus_vec(3, this%N))
    allocate(this%lift(this%N))
    allocate(this%drag(this%N))

    dy = this%b / this%N

    do j = 1, this%N
      this%x_vec%x(:, j) = [this%x_center, this%y_center - 0.5_rp * this%b + dy * (j - 0.5_rp), this%z_center]
      this%xminus_vec(:, j) = this%x_vec%x(:, j) - [0.0_rp, 0.5_rp * dy, 0.0_rp]
      this%xplus_vec(:, j) = this%x_vec%x(:, j) + [0.0_rp, 0.5_rp * dy, 0.0_rp]
    end do
    
  end subroutine discretize_line

  !> Computes Gaussian kernel, in the entire domain, for each actuator line element.
  subroutine compute_kernel(this)
    class(actuator_line_source_term_t), intent(inout) :: this

    integer :: i, j
    real(kind=rp) :: x_vec_ij(3), xplus(3), xminus(3), xp_xm(3)
    real(kind=rp) :: x, y, z, lminus, lplus, d, xp_xm_norm
    type(field_t), pointer :: u
    character(len=:), allocatable :: kernel_name

    ! Allocate kernels
    u => neko_registry%get_field_by_name('u') ! Use any field
    allocate(this%kernel(u%dof%size(), this%N))
    
    call this%make_registry_name("kernel", kernel_name)
    call neko_registry%add_field(u%dof, kernel_name, .false.)
    this%kernel_j => neko_registry%get_field(kernel_name)

    ! Loop through all actuator line elements
    do j = 1, this%N
      xplus = this%xplus_vec(:, j)
      xminus = this%xminus_vec(:, j)

      xp_xm = xplus - xminus
      xp_xm_norm = norm2(xp_xm)

      ! Loop through all mesh nodes
      do i = 1, u%dof%size()
        ! Get coordinates
        x = u%dof%x(i,1,1,1)
        y = u%dof%y(i,1,1,1)
        z = u%dof%z(i,1,1,1)

        x_vec_ij = [x, y, z]

        lplus = dot_product((x_vec_ij - xplus), xp_xm) / xp_xm_norm
        lminus = lplus + xp_xm_norm
        d = sqrt(norm2(x_vec_ij - xminus)**2 - lminus**2)

        this%kernel(i,j) = (1.0 / (pi * this%eps**2)) * exp(- d**2 / this%eps**2) * &
                           (erf(lminus / this%eps) - erf(lplus / this%eps)) / 2.0
      end do
    end do

  end subroutine compute_kernel

  !> Initialize interpolator and find points' owners
  subroutine init_interpolator(this)
    class(actuator_line_source_term_t), intent(inout) :: this
    character(len=:), allocatable :: u_name, v_name, w_name

    call this%make_registry_name("u_interp", u_name)
    call this%make_registry_name("v_interp", v_name)
    call this%make_registry_name("w_interp", w_name)

    call this%interpolator%find_points(this%x_vec%x, this%N) ! Find points and their owners

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
    real(kind=rp) :: u_vec(3)
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
      call cross(u_vec, [0.0_rp, this%gamm(j), 0.0_rp], F_2D(j,:))
      this%lift(j) = F_2D(j,3)
      this%drag(j) = F_2D(j,1)
    end do

    ! Print and save to coefs file
    call this%compute_coefficients()
    write(this%coefs_out_num, *) time%t, ",", this%coefs%x(1), ",", this%coefs%x(2)
  end subroutine F2D

  !> Compute 3D lift and drag coefficients 
  subroutine compute_coefficients(this)
    class(actuator_line_source_term_t), intent(in) :: this

    this%coefs%x(1) = this%AR * sum(this%lift / this%N) / (0.5_rp * this%b) ! CL
    this%coefs%x(2) = this%AR * sum(this%drag / this%N) / (0.5_rp * this%b) ! CD

  end subroutine compute_coefficients

  subroutine make_registry_name(this, base_name, registry_name)
    class(actuator_line_source_term_t), intent(in) :: this
    character(len=*), intent(in) :: base_name
    character(len=:), allocatable, intent(out) :: registry_name
    character(len=64) :: temp_name

    write(temp_name, '("alm_", A, "_", I0)') trim(base_name), this%instance_id
    registry_name = trim(temp_name)

  end subroutine make_registry_name

  !> Computes the cross product of two 3D vectors.
  subroutine cross(a, b, c)
    real(kind=rp), intent(in) :: a(3), b(3)
    real(kind=rp), intent(out) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
  end subroutine cross

end module actuator_line_source_term
