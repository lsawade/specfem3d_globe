!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!----
!---- Green function database: expand tagged elements to include neighbor shells
!----

  subroutine gf_expand_neighbors()

! Expands the set of tagged elements by GF_NEIGHBOR_SHELLS layers of neighbors.
!
! Local expansion uses the xadj/adjncy CSR adjacency graph (within-rank).
! Cross-MPI expansion uses ibool_interfaces_crust_mantle to find elements
! on neighboring ranks that share nodes with tagged boundary elements.
!
! Must be called BEFORE xadj/adjncy are deallocated.

  use constants_solver, only: IMAIN,NGLLX,NGLLY,NGLLZ,NPROCTOT_VAL,itag

  use shared_parameters, only: GF_NEIGHBOR_SHELLS

  use specfem_par, only: &
    myrank, &
    nspec => NSPEC_CRUST_MANTLE, &
    nglob => NGLOB_CRUST_MANTLE, &
    xadj, adjncy

  use specfem_par, only: &
    num_interfaces_crust_mantle, &
    max_nibool_interfaces_cm, &
    my_neighbors_crust_mantle, &
    nibool_interfaces_crust_mantle, &
    ibool_interfaces_crust_mantle

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle

  use green_function_par, only: &
    gf_nelem_local, gf_local_elements

  implicit none

  ! local parameters
  integer :: ishell, ispec, j, ier
  integer :: nelem_after

  ! boolean mask for tagged elements (O(n) deduplication)
  logical, allocatable, dimension(:) :: ispec_tagged

  ! nothing to do if no expansion requested or no local elements
  if (GF_NEIGHBOR_SHELLS <= 0) then
    if (myrank == 0) then
      write(IMAIN,*) 'Green function database: GF_NEIGHBOR_SHELLS = 0, no expansion'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    return
  endif

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' Green function database: expanding neighbor shells'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,*) 'GF_NEIGHBOR_SHELLS = ',GF_NEIGHBOR_SHELLS
    call flush_IMAIN()
  endif

  ! build boolean mask from current tagged elements
  allocate(ispec_tagged(nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating ispec_tagged in gf_expand_neighbors')
  ispec_tagged(:) = .false.

  do j = 1, gf_nelem_local
    ispec_tagged(gf_local_elements(j)) = .true.
  enddo

  ! expand by GF_NEIGHBOR_SHELLS layers
  do ishell = 1, GF_NEIGHBOR_SHELLS

    ! --- local expansion using xadj/adjncy ---
    call gf_expand_local(nspec,ispec_tagged,xadj,adjncy)

    ! --- cross-MPI expansion ---
    if (NPROCTOT_VAL > 1) then
      call gf_expand_cross_mpi(nspec,nglob,ispec_tagged,ibool, &
                               num_interfaces_crust_mantle, &
                               max_nibool_interfaces_cm, &
                               my_neighbors_crust_mantle, &
                               nibool_interfaces_crust_mantle, &
                               ibool_interfaces_crust_mantle)
    endif

  enddo

  ! rebuild gf_local_elements from updated mask
  nelem_after = count(ispec_tagged)

  ! reallocate if size changed
  if (nelem_after /= gf_nelem_local) then
    if (allocated(gf_local_elements)) deallocate(gf_local_elements)
    gf_nelem_local = nelem_after
    if (gf_nelem_local > 0) then
      allocate(gf_local_elements(gf_nelem_local),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating gf_local_elements after expansion')
      j = 0
      do ispec = 1, nspec
        if (ispec_tagged(ispec)) then
          j = j + 1
          gf_local_elements(j) = ispec
        endif
      enddo
    endif
  endif

  deallocate(ispec_tagged)

  ! gather and report per-rank element counts and IDs on rank 0
  call gf_report_tagged_elements('after neighbor expansion')

  end subroutine gf_expand_neighbors

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_expand_local(nspec,ispec_tagged,xadj,adjncy)

! Expands tagged elements by one layer using the local xadj/adjncy CSR adjacency.
! Uses a snapshot of currently tagged elements to avoid chain expansion within
! a single layer.

  implicit none

  integer, intent(in) :: nspec
  logical, dimension(nspec), intent(inout) :: ispec_tagged
  integer, dimension(nspec+1), intent(in) :: xadj
  integer, dimension(*), intent(in) :: adjncy

  ! local parameters
  integer :: ispec, iadj, ineighbor

  ! snapshot: which elements are tagged at the START of this layer
  logical, allocatable :: tagged_snapshot(:)
  integer :: ier

  allocate(tagged_snapshot(nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating tagged_snapshot'
  tagged_snapshot(:) = ispec_tagged(:)

  do ispec = 1, nspec
    if (.not. tagged_snapshot(ispec)) cycle

    ! expand to all neighbors of this element
    do iadj = xadj(ispec) + 1, xadj(ispec + 1)
      ineighbor = adjncy(iadj)
      if (ineighbor >= 1 .and. ineighbor <= nspec) then
        ispec_tagged(ineighbor) = .true.
      endif
    enddo
  enddo

  deallocate(tagged_snapshot)

  end subroutine gf_expand_local

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_expand_cross_mpi(nspec,nglob,ispec_tagged,ibool, &
                                 num_interfaces,max_nibool_interfaces, &
                                 my_neighbors,nibool_interfaces, &
                                 ibool_interfaces)

! Cross-MPI neighbor expansion.
!
! For tagged boundary elements, finds which of their GLL nodes are shared
! with neighboring MPI ranks via ibool_interfaces. Sends those shared node
! indices to the neighbor rank, which then tags any of its elements that
! contain those nodes.
!
! Uses send_i/recv_i with careful ordering (lower rank sends first) to
! avoid deadlock.

  use constants_solver, only: NGLLX,NGLLY,NGLLZ,itag

  use specfem_par, only: myrank

  implicit none

  integer, intent(in) :: nspec, nglob
  logical, dimension(nspec), intent(inout) :: ispec_tagged
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  integer, intent(in) :: num_interfaces, max_nibool_interfaces
  integer, dimension(num_interfaces), intent(in) :: my_neighbors, nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces), intent(in) :: ibool_interfaces

  ! local parameters
  integer :: iinterface, ipoin, iglob, ispec, i, j, k, ier
  integer :: nsend, nrecv, neighbor_rank
  integer :: ntagged_on_interface

  ! node-is-tagged lookup (built from tagged elements)
  logical, allocatable :: node_tagged(:)

  ! send/recv buffers for shared node flags
  integer, allocatable :: send_flags(:), recv_flags(:)

  ! reverse lookup: node -> element mapping for boundary nodes only
  ! We use ibool to find which elements contain received nodes
  logical, allocatable :: node_is_shared(:)

  ! nothing to do with single process
  if (num_interfaces == 0) return

  ! build node_tagged: mark all GLL nodes of tagged elements
  allocate(node_tagged(nglob),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating node_tagged')
  node_tagged(:) = .false.

  do ispec = 1, nspec
    if (.not. ispec_tagged(ispec)) cycle
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          node_tagged(ibool(i,j,k,ispec)) = .true.
        enddo
      enddo
    enddo
  enddo

  ! for each MPI interface, exchange flags indicating which shared nodes
  ! belong to tagged elements
  do iinterface = 1, num_interfaces
    neighbor_rank = my_neighbors(iinterface)
    nsend = nibool_interfaces(iinterface)
    nrecv = nsend  ! symmetric interface

    allocate(send_flags(nsend), recv_flags(nrecv), stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating MPI exchange buffers')

    ! pack: 1 if this shared node belongs to a tagged element, 0 otherwise
    do ipoin = 1, nsend
      iglob = ibool_interfaces(ipoin, iinterface)
      if (node_tagged(iglob)) then
        send_flags(ipoin) = 1
      else
        send_flags(ipoin) = 0
      endif
    enddo

    ! exchange using send/recv with ordering to avoid deadlock
    if (myrank < neighbor_rank) then
      call send_i(send_flags, nsend, neighbor_rank, itag)
      call recv_i(recv_flags, nrecv, neighbor_rank, itag)
    else
      call recv_i(recv_flags, nrecv, neighbor_rank, itag)
      call send_i(send_flags, nsend, neighbor_rank, itag)
    endif

    ! process received flags: for shared nodes flagged by the neighbor,
    ! find which local elements contain those nodes and tag them
    ntagged_on_interface = 0
    do ipoin = 1, nrecv
      if (recv_flags(ipoin) == 1) then
        iglob = ibool_interfaces(ipoin, iinterface)
        ntagged_on_interface = ntagged_on_interface + 1
        ! mark this node so we can find its elements
        node_tagged(iglob) = .true.
      endif
    enddo

    ! if any nodes were flagged, find elements containing them
    if (ntagged_on_interface > 0) then
      ! build a quick set of received flagged node IDs
      allocate(node_is_shared(nglob), stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating node_is_shared')
      node_is_shared(:) = .false.

      do ipoin = 1, nrecv
        if (recv_flags(ipoin) == 1) then
          iglob = ibool_interfaces(ipoin, iinterface)
          node_is_shared(iglob) = .true.
        endif
      enddo

      ! scan all elements: if any corner node is in the received set, tag it
      do ispec = 1, nspec
        if (ispec_tagged(ispec)) cycle  ! already tagged

        ! check 8 corner nodes (sufficient for adjacency — same as xadj/adjncy)
        if (node_is_shared(ibool(1,1,1,ispec)) .or. &
            node_is_shared(ibool(NGLLX,1,1,ispec)) .or. &
            node_is_shared(ibool(NGLLX,NGLLY,1,ispec)) .or. &
            node_is_shared(ibool(1,NGLLY,1,ispec)) .or. &
            node_is_shared(ibool(1,1,NGLLZ,ispec)) .or. &
            node_is_shared(ibool(NGLLX,1,NGLLZ,ispec)) .or. &
            node_is_shared(ibool(NGLLX,NGLLY,NGLLZ,ispec)) .or. &
            node_is_shared(ibool(1,NGLLY,NGLLZ,ispec))) then
          ispec_tagged(ispec) = .true.
        endif
      enddo

      deallocate(node_is_shared)
    endif

    deallocate(send_flags, recv_flags)
  enddo

  deallocate(node_tagged)

  end subroutine gf_expand_cross_mpi
