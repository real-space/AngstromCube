// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // uint16_t
#include <vector> // std::vector<T>

#include "green_parallel.hxx" // rank_int_t

#include "simple_stats.hxx" // ::Stats<double>
#include "inline_math.hxx" // set
#include "mpi_parallel.hxx" // ::init, ::size, ::rank, ::finalize, ::max, ::allreduce
#include "sho_tools.hxx" // ::nSHO
#include "global_coordinates.hxx" // ::get
#include "print_tools.hxx" // printf_vector
#include "recorded_warnings.hxx" // warn
#include "data_view.hxx" // view3D

namespace green_parallel {

  int init(int argc, char **argv) {
      return mpi_parallel::init(argc, argv);
  } // init

  MPI_Comm comm() { return MPI_COMM_WORLD; }

  unsigned size(void) {
      return mpi_parallel::size(); }

  int rank(void) {
      return mpi_parallel::rank(); }

  int finalize(void) {
      return mpi_parallel::finalize(); }

  int max(uint16_t data[], size_t const n) {
      return mpi_parallel::max(data, n); }

  int allreduce(simple_stats::Stats<double> & stats) {
      return mpi_parallel::allreduce(stats, comm()); }

  inline char const * spin_name(int const Noco, int const spin) {
      if (1 == Noco && 0 == spin) return "";
      if (2 == Noco) {
          switch (spin) {
            case 0: return " V_down";
            case 1: return " V_up";
            case 2: return " V_x";
            case 3: return " V_y";
          } // spin
      } // 2 == Noco
      return " ???";
  } // spin_name

  // For MPI parallel calculations the potential values need to be exchanged

  status_t potential_exchange(
        double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco^2][nrows][64]
      , std::vector<int64_t> const & requests  // indices requested by this MPI process,         [nrows]
      , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco^2 ][64]
      , std::vector<int64_t> const & offerings //  indices offered  by this MPI process, [ncols]
      , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
      , uint32_t const nb[3] // global bounding box, maybe group this with owner_rank into a translator object global_index --> owner_rank
      , int const Noco // =1, 1:no spin, 2: (non-collinear) spin
      , bool const debug // =true
      , int const echo // =0, log-level
  ) {
      if (echo > 0) std::printf("# old MPI exchange of potential, MPI one-sided communication, Noco=%d\n", Noco);
      assert(1 == Noco || 2 == Noco);

      assert(Veff && "may not be called with a nullptr for output");
      for (int spin = 0; spin < Noco*Noco; ++spin) {
          assert(Veff[spin] && "may not be called with a nullptr for spin output");
      } // spin
      assert(Vinp && "may not be called with a nullptr for input");

      status_t status(0);

      auto const comm = MPI_COMM_WORLD;
      auto const me = mpi_parallel::rank(comm);
      auto const ncols = offerings.size(); // number of offerings
      auto const nrows =  requests.size(); // number of requests
      size_t stats[] = {0, 0};

#ifndef   HAS_NO_MPI

      int constexpr X=0, Y=1, Z=2;
      auto const nall = nb[Z]*size_t(nb[Y])*size_t(nb[X]);

      std::vector<uint16_t> local_index(nall, 0);
      std::vector<uint16_t> local_check(nall*int(debug), 0);

      assert(ncols <= (1ul << 16) && "each process can hold max 2^16 Green function columns");
      for (size_t col = 0; col < ncols; ++col) {
          auto const global_id = offerings[col];

          // translate global_id into an box index iall
          uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
          for (int d = 0; d < 3; ++d) {
              assert(xyz[d] < nb[d]);
          } // d
          auto const iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
          assert(iall < nall);

          assert(me == owner_rank[iall] && "offerings must be owned");

          local_index[iall] = col;
          if (debug) ++local_check[iall];
      } // col

      if (debug) {
          if (echo > 7) { std::printf("# local_check before "); printf_vector(" %i", local_check); }
          for (size_t iall = 0; iall < nall; ++iall) {
              assert(local_check[iall] <= 1 && "duplicates found");
          } // iall
          status += MPI_Allreduce(MPI_IN_PLACE, local_check.data(), nall, MPI_UINT16, MPI_SUM, comm);

          if (echo > 7) { std::printf("# local_check after  "); printf_vector(" %i", local_check); }
          for (size_t iall = 0; iall < nall; ++iall) {
              assert(1 == local_check[iall] && "not all covered");
          } // iall
          local_check.resize(0); // DEBUG
      } // debug

      // get a global list of which local index is where
      status += MPI_Allreduce(MPI_IN_PLACE, local_index.data(), nall, MPI_UINT16, MPI_MAX, comm);

      // set up a memory window to read from
      MPI_Win window;
      size_t const disp_unit = 64*sizeof(double); // in Bytes
      size_t const win_size  = ncols*Noco*Noco*disp_unit; // in Bytes
      status += MPI_Win_create((void*)Vinp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);

      // synchronize processes
      int const assertions = MPI_MODE_NOPUT; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
      status += MPI_Win_fence(assertions, window);

#endif // HAS_NO_MPI

      // look up who owns the requested cube
      for (size_t row = 0; row < nrows; ++row) {
          auto const global_id = requests[row];

#ifndef   HAS_NO_MPI

          // translate global_id into an owner rank
          uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
          for (int d = 0; d < 3; ++d) {
              assert(xyz[d] < nb[d]);
          } // d
          auto const iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
          assert(iall < nall);

          auto const owner = owner_rank[iall];
          auto const iloc = local_index[iall];

#else  // HAS_NO_MPI

          // version without MPI
          auto const owner = me;
          // search global_id in offerings
          int64_t iloc{-1};
          for (size_t col = 0; col < ncols && iloc < 0; ++col) {
              if (global_id == offerings[col]) iloc = col;
          } // col
          assert(iloc > -1 && "index not found in local offerings");

#endif // HAS_NO_MPI

          for (int spin = 0; spin < Noco*Noco; ++spin) {
              auto const target_disp = iloc*Noco*Noco + spin;
              if (me == owner) {
                  if (echo > 7 + 5*spin) { std::printf("# exchange: rank#%i get potential at cube id o%llo  copy local element %i\n", me, global_id, iloc); std::fflush(stdout); }
                  set(Veff[spin][row], 64, Vinp[target_disp]);
                  ++stats[0];
              } else { // me == owner
                  ++stats[1];
#ifndef   HAS_NO_MPI
                  if (echo > 7 + 5*spin) { std::printf("# exchange: rank#%i get potential at cube id o%llo from rank#%i element %i\n", me, global_id, owner, iloc); std::fflush(stdout); }
                  status += MPI_Get(Veff[spin][row], 64, MPI_DOUBLE, owner, target_disp, 64, MPI_DOUBLE, window);
                  /*
                   *  Criticism for MPI_Get-solution:
                   *    each time called it pushes a request into the remote process and receives an MPI_Put with the data from that process.
                   *    After the construction of the Green function, the structure of requested potential elements does not change any more.
                   *    Latencies would be lower if we separated it into two phases:
                   *    a) communicate where to push, e.g. with MPI_Alltoall and MPI_Alltoallv before the loop over energy parameters.
                   *    b) MPI_Put the data
                   */
#else  // HAS_NO_MPI
                  error("Without MPI all potential elements must reside in the same process, me=%i, owner=%i", me, owner);
#endif // HAS_NO_MPI
              } // me == owner
          } // spin
      } // row

#ifndef   HAS_NO_MPI
      // synchronize processes
      status += MPI_Barrier(comm); // { std::printf("# rank#%d hit MPI_Barrier at line %d\n", me, __LINE__); std::fflush(stdout); }
      status += MPI_Win_fence(assertions, window);
      status += MPI_Win_free(&window);
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == Noco*Noco*nrows);

      return status;
  } // potential_exchange (deprecated)


  // For MPI parallel calculations the atom matrix values need to be exchanged
  status_t dyadic_exchange(                                    // here, nSHO is the global maximum of nsho[ia]
        double       *const mat_out // output effective atom matrices, data layout mat_out[nrows*Noco^2*2*nSHO^2]
      , std::vector<int64_t> const & requests  // indices requested by this MPI process,  [nrows]
      , double const *const mat_inp //  input effective atom matrices, data layout mat_inp[ncols*Noco^2*2*nSHO^2]
      , std::vector<int64_t> const & offerings // indices offered   by this MPI process,  [ncols]
      , rank_int_t const owner_rank[] // where to find it, [nall]
      , uint32_t const nall // number of all atoms
      , int const count // number of doubles per package
      , bool const debug // =true
      , int const echo // =0, log-level
  ) {
      // The number of local atoms is limited to 2^16 == 65536
      if (echo > 0) std::printf("# old MPI exchange of atom matrices, MPI one-sided communication, packages of %.3f k doubles\n", count*.001);

      assert(mat_out && "may not be called with a nullptr for output");
      assert(mat_inp && "may not be called with a nullptr for input");

      status_t status(0);

      auto const comm = MPI_COMM_WORLD;
      auto const me = mpi_parallel::rank(comm);
      auto const ncols = offerings.size(); // number of offerings
      auto const nrows =  requests.size(); // number of requests
      size_t stats[] = {0, 0};

#ifndef   HAS_NO_MPI

      std::vector<uint16_t> local_index(nall, 0);
      std::vector<uint16_t> local_check(nall*int(debug), 0);

      assert(ncols <= (1ul << 16) && "each process can hold max 2^16 locally owner atoms");
      for (size_t col = 0; col < ncols; ++col) {
          auto const global_id = offerings[col];

          // translate global_id into an index iall
          assert(global_id > -1);
          auto const iall = global_id;
          assert(iall < nall);

          assert(me == owner_rank[iall] && "offerings must be owned");

          local_index[iall] = col;
          if (debug) ++local_check[iall];
      } // col

      if (debug) {
          if (echo > 7) { std::printf("# local_check before "); printf_vector(" %i", local_check); }
          for (size_t iall = 0; iall < nall; ++iall) {
              assert(local_check[iall] <= 1 && "duplicates found");
          } // iall
          status += MPI_Allreduce(MPI_IN_PLACE, local_check.data(), nall, MPI_UINT16, MPI_SUM, comm);
//        status += mpi_parallel::sum(local_check.data(), nall, comm);

          if (echo > 7) { std::printf("# local_check after  "); printf_vector(" %i", local_check); }
          for (size_t iall = 0; iall < nall; ++iall) {
              assert(1 == local_check[iall] && "not all covered");
          } // iall
          local_check.resize(0); // DEBUG
      } // debug

      // get a global list of which local index is where
      status += MPI_Allreduce(MPI_IN_PLACE, local_index.data(), nall, MPI_UINT16, MPI_MAX, comm);
//    status += mpi_parallel::max(local_index.data(), nall, comm);

      // set up a memory window to read from
      MPI_Win window;
      uint const disp_unit = count*sizeof(double); // in Bytes
      size_t const win_size  = ncols*disp_unit; // in Bytes
      status += MPI_Win_create((void*)mat_inp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);

      // synchronize processes
      int const assertions = 0; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
      status += MPI_Win_fence(assertions, window);

#endif // HAS_NO_MPI

      // look up who owns the requested cube
      for (size_t row = 0; row < nrows; ++row) {
          auto const global_id = requests[row];

#ifndef   HAS_NO_MPI

          // translate global_id into an index iall
          assert(global_id > -1);
          auto const iall = global_id;
          assert(iall < nall);

          auto const owner = owner_rank[iall];
          auto const iloc = local_index[iall];

#else  // HAS_NO_MPI

          // version without MPI
          auto const owner = me;
          // search global_id in offerings
          int64_t iloc{-1};
          for (size_t col = 0; col < ncols && iloc < 0; ++col) {
              if (global_id == offerings[col]) iloc = col;
          } // col
          assert(iloc > -1 && "index not found in local offerings");

#endif // HAS_NO_MPI

          if (me == owner) {
              if (echo > 7) std::printf("# exchange: rank#%i get matrices of atom#%i  copy local element %i\n", me, global_id, iloc);
              set(&mat_out[row*count], count, &mat_inp[iloc*count]);
              ++stats[0];
          } else { // me == owner
              ++stats[1];
#ifndef   HAS_NO_MPI
              if (echo > 7) std::printf("# exchange: rank#%i get matrices of atom#%i from rank#%i element %i\n", me, global_id, owner, iloc);
              status += MPI_Get(&mat_out[row*count], count, MPI_DOUBLE, owner, iloc, count, MPI_DOUBLE, window);
              /*
                *  Criticism for MPI_Get-solution:
                *    each time called it pushes a request into the remote process and receives an MPI_Put with the data from that process.
                *    After the construction of the Green function, the structure of requested potential elements does not change any more.
                *    Latencies would be lower if we separated it into two phases:
                *    a) communicate where to push, e.g. with MPI_Alltoall and MPI_Alltoallv before the loop over energy parameters.
                *    b) MPI_Put the data
                */
#else  // HAS_NO_MPI
              error("Without MPI all atom matrices must reside in the same process, me=%i, owner=%i", me, owner);
#endif // HAS_NO_MPI
          } // me == owner
      } // row

#ifndef   HAS_NO_MPI
      // synchronize processes
      status += MPI_Barrier(comm); // { std::printf("# rank#%d hit MPI_Barrier at line %d\n", me, __LINE__); std::fflush(stdout); }
      status += MPI_Win_fence(assertions, window);
      status += MPI_Win_free(&window);
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == nrows);

      return status;
  } // dyadic_exchange (deprecated)





    // New Grouping: requests are first sorted into a RequestList_t and then, data exchange happens

    RequestList_t::RequestList_t( // constructor implementation
        std::vector<int64_t> const & requests
      , std::vector<int64_t> const & offerings // do we need this at all? the offerings should match owner_rank[]==me anyway
      , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
      , uint32_t const nb[3] // global bounding box or {natoms,0,0}
      , int const echo // =0 // log-level
    ) {
        int constexpr X=0, Y=1, Z=2;
        auto const grid = size_t(nb[Z])*size_t(nb[Y])*size_t(nb[X]);
        auto const nall = grid ? grid : nb[X] + nb[Y] + nb[Z];
        if (echo > 7) std::printf("# RequestList_t %d %d %d, nall= %ld\n", nb[X], nb[Y], nb[Z], nall);
        auto const nown = offerings.size(); // number of offerings
        auto const nreq = requests.size(); // number of requests

#ifndef   HAS_NO_MPI
        int constexpr debug = 1;
        auto const comm = MPI_COMM_WORLD;
        auto const me = mpi_parallel::rank(comm);

        std::vector<uint16_t> local_index(nall, 0);
        std::vector<uint16_t> local_check(nall*int(debug), 0);

        assert(nown <= (1ul << 16) && "each process can hold max 2^16 locally owned items");
        for (size_t iown = 0; iown < nown; ++iown) {
            auto const global_id = offerings[iown];

            // translate global_id into an index iall
            assert(global_id > -1);
            size_t iall = global_id;
            if (grid) {
                uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
                for (int d = 0; d < 3; ++d) assert(xyz[d] < nb[d] && "requested coordinates exceed box");
                iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
            } // grid
            assert(iall < nall);

            assert(me == owner_rank[iall] && "all offerings must be owned");

            local_index[iall] = iown;
            if (debug) ++local_check[iall];
        } // iown

        if (debug) {
            if (echo > 7) { std::printf("# local_check before "); printf_vector(" %i", local_check); }
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(local_check[iall] <= 1 && "duplicates found");
            } // iall

            auto const stat = MPI_Allreduce(MPI_IN_PLACE, local_check.data(), nall, MPI_UINT16, MPI_SUM, comm);
            if (stat) warn("MPI_Allreduce failed with status= %d", int(stat));

            if (echo > 7) { std::printf("# local_check after  "); printf_vector(" %i", local_check); }
            for (size_t iall = 0; iall < nall; ++iall) {
                assert(1 == local_check[iall] && "not all covered");
            } // iall
            local_check.resize(0); // DEBUG
        } // debug

        // get a global list of which local index is where
        auto const stat = MPI_Allreduce(MPI_IN_PLACE, local_index.data(), nall, MPI_UINT16, MPI_MAX, comm);
        if (stat) warn("MPI_Allreduce failed with status= %d", int(stat));
        // if this is too expensive see ALTERNATIVE
        //
        // ALTERNATIVE:
        // use MPI_Send and MPI_Recv to distribute the necessary info about local_index 
        //

#endif // HAS_NO_MPI

        // initialize member fields
        owner = std::vector<int32_t>(nreq, 0); // initialize with master rank for the serial version
        index = std::vector<uint32_t>(nreq, 0);
        requested_id = std::vector<int64_t>(nreq, 0);
        window_size = nown; // number of owned data items

        for (size_t ireq = 0; ireq < nreq; ++ireq) {
            auto const global_id = requests[ireq];

#ifndef   HAS_NO_MPI

            // translate global_id into an index iall
            assert(global_id > -1 && "invalid request id");
            size_t iall = global_id;
            if (grid) {
                uint32_t xyz[3]; global_coordinates::get(xyz, global_id);
                for (int d = 0; d < 3; ++d) assert(xyz[d] < nb[d] && "requested coordinates exceed box");
                iall = (xyz[Z]*size_t(nb[Y]) + xyz[Y])*nb[X] + xyz[X];
            } // grid
            assert(iall < nall && "internal index exceeded");

            owner[ireq]     =  owner_rank[iall];
            auto const iloc = local_index[iall];

#else  // HAS_NO_MPI

            // without MPI search global_id in offerings
            int64_t iloc{-1};
            for (size_t iown = 0; iown < nown && iloc < 0; ++iown) {
                if (global_id == offerings[iown]) iloc = iown;
            } // col
            assert(iloc > -1 && "index not found in local offerings");

#endif // HAS_NO_MPI

            index[ireq] = iloc;
            requested_id[ireq] = global_id; // copy the array of requests
        } // ireq

    } // constructor implementation


  status_t exchange(
        double       *const data_out // output data, data layout data_out[nrequests*count]
      , double const *const data_inp //  input data, data layout data_inp[nowned   *count]
      , RequestList_t const & requests
      , int const count // number of doubles per package
      , int const echo // =0, log-level
  ) {
      // The number of local atoms is limited to 2^16 == 65536
      if (echo > 0) std::printf("# exchange using MPI one-sided communication, packages of %.3f k doubles\n", count*.001);

      assert(data_out && "may not be called with a nullptr for output");
      assert(data_inp && "may not be called with a nullptr for input");

      status_t status(0);

      auto const comm = MPI_COMM_WORLD;
      auto const me = mpi_parallel::rank(comm);

#ifndef   HAS_NO_MPI

      // set up a memory window to read from
      MPI_Win window;
      uint const disp_unit = count*sizeof(double); // in Bytes
      size_t const win_size = requests.window()*disp_unit; // in Bytes
      status += MPI_Win_create((void*)data_inp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);

      // synchronize processes
      int const assertions = 0; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
      status += MPI_Win_fence(assertions, window);

#endif // HAS_NO_MPI

      size_t stats[] = {0, 0}; // get element from {0:local, 1:remote}
      auto const nreq = requests.size(); // number of requests
      for (size_t ireq = 0; ireq < nreq; ++ireq) {
          auto const global_id = requests.requested_id[ireq];
          auto const owner     = requests.owner[ireq];
          auto const iloc      = requests.index[ireq];

          if (me == owner) {
              if (echo > 7) std::printf("# exchange: rank#%i get matrices of atom#%i  copy local element %i\n", me, global_id, iloc);
              set(&data_out[ireq*count], count, &data_inp[iloc*count]); // copy
              ++stats[0]; // local
          } else { // me == owner
              ++stats[1]; // remote
#ifndef   HAS_NO_MPI
              if (echo > 7) std::printf("# exchange: rank#%i get matrices of atom#%i from rank#%i element %i\n", me, global_id, owner, iloc);
              status += MPI_Get(&data_out[ireq*count], count, MPI_DOUBLE, owner, iloc, count, MPI_DOUBLE, window);
#else  // HAS_NO_MPI
              error("Without MPI all atom matrices must reside in the same process, me=%i, owner=%i", me, owner);
#endif // HAS_NO_MPI
          } // me == owner
      } // ireq

#ifndef   HAS_NO_MPI
      // synchronize processes
      status += MPI_Barrier(comm); // { std::printf("# rank#%d hit MPI_Barrier at line %d\n", me, __LINE__); std::fflush(stdout); }
      status += MPI_Win_fence(assertions, window);
      status += MPI_Win_free(&window);
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == nreq);

      return status;
  } // exchange


  status_t potential_exchange(
        double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco^2][nreq][64]
      , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco^2][64]
      , RequestList_t const & requests
      , int const Noco // =1, 1:no spin, 2: (non-collinear) spin
      , int const echo // =0, log-level
  ) {
      if (echo > 0) std::printf("# new MPI exchange of potential, MPI one-sided communication, Noco=%d\n", Noco);
      assert(1 == Noco || 2 == Noco);

      assert(Veff && "may not be called with a nullptr for output");
      for (int spin = 0; spin < Noco*Noco; ++spin) {
          assert(Veff[spin] && "may not be called with a nullptr for spin output");
      } // spin
      assert(Vinp && "may not be called with a nullptr for input");

      auto const nreq = requests.size();
      view3D<double> Vout(nreq,Noco*Noco,64, 0.0); // get temporary CPU memory in [Noco^2][64] layout

      auto const status = exchange(Vout(0,0), Vinp[0], requests, Noco*Noco*64, echo);

      // convert Vout into special data layout of Veff (GPU memory) 
      for (size_t ireq = 0; ireq < nreq; ++ireq) {
          for (int spin = 0; spin < Noco*Noco; ++spin) {
              set(Veff[spin][ireq], 64, Vout(ireq,spin));
          } // spin
      } // ireq

      return status;
  } // potential_exchange with RequestList_t



} // namespace green_parallel

#include "status.hxx" // status_t

namespace green_parallel {

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  status_t test_exchanges(int echo=0, int const nSHO=20) {
      status_t stat(0);
      uint32_t const nb[] = {2, 2, 2}; // grid box is 2x2x2 or 8 atoms
      std::vector<int64_t> requests = {0,7,6,1,5,2,4,3}; // every process requests all of these ids

      auto const nrows = requests.size();
      int const me = mpi_parallel::rank(),
                np = mpi_parallel::size(); assert(np > 0);

      std::vector<int64_t> offerings(0);
      uint32_t const nall[] = {nb[2]*nb[1]*nb[0], 0, 0};
      std::vector<uint16_t> owner_rank(nall[0], 0);
      for(int id{0}; id < nall[0]; ++id) {
          int const rank = id % np; // block-cyclic distribution
          owner_rank[id] = rank;
          if (me == rank) offerings.push_back(id);
      } // id
      int const ncols = offerings.size();
      RequestList_t rlV(requests, offerings, owner_rank.data(), nb, echo);
      RequestList_t rlD(requests, offerings, owner_rank.data(), nall, echo);

      double (*pot_out[2*2])[64]; view3D<double> ptr(4,nrows,64);
      for (int spin = 0; spin < 2*2; ++spin) { pot_out[spin] = (double(*)[64]) ptr(spin,0); }

      int nerr{0};
      for (int Noco = 1; Noco <= 2; ++Noco) {
          if (1) {
            auto const pot_inp = new double[nrows*Noco*Noco][64];
            for (int col = 0; col < ncols; ++col) pot_inp[col*Noco*Noco][0] = 0.5 + me; // ear-mark with owner rank
            stat += green_parallel::potential_exchange(pot_out, requests, pot_inp, offerings, owner_rank.data(), nb, Noco, true, echo);
            for (int row = 0; row < nrows; ++row) nerr += (pot_out[0][row][0] != (0.5 + owner_rank[requests[row]]));
            if (echo > 0) std::printf("# potential_exchange Noco= %d status= %d errors= %d\n\n", Noco, int(stat), nerr);
            stat += green_parallel::potential_exchange(pot_out, pot_inp, rlV, Noco, echo);
            for (int row = 0; row < nrows; ++row) nerr += (pot_out[0][row][0] != (0.5 + owner_rank[requests[row]]));
            if (echo > 0) std::printf("# new potential_exchange Noco= %d status= %d errors= %d\n\n", Noco, int(stat), nerr);
            delete[] pot_inp;
          }

          int const count = Noco*Noco*2*nSHO*nSHO; // number of doubles per package
          std::vector<double> mat_out(nrows*count), mat_inp(ncols*count);
          for (int col = 0; col < ncols; ++col) mat_inp[col*count] = 0.5 + me; // ear-mark with owner rank
          stat += green_parallel::dyadic_exchange(mat_out.data(), requests, mat_inp.data(), offerings, owner_rank.data(), nall[0], count, true, echo);
          for (int row = 0; row < nrows; ++row) nerr += (mat_out[row*count] != (0.5 + owner_rank[requests[row]]));
          if (echo > 0) std::printf("# dyadic_exchange Noco= %d status= %d errors= %d\n\n", Noco, int(stat), nerr);
          stat += green_parallel::exchange(mat_out.data(), mat_inp.data(), rlD, count, echo);
          for (int row = 0; row < nrows; ++row) { nerr += (mat_out[row*count] != (0.5 + owner_rank[requests[row]])); }
          if (echo > 0) std::printf("# new exchange Noco= %d status= %d errors= %d\n\n", Noco, int(stat), nerr);
      } // Noco
      return stat;
  } // test_exchanges

  status_t all_tests(int const echo) {
      status_t stat(0);
      bool const already_initialized = mpi_parallel::init();
      stat += test_exchanges(echo);
      if (!already_initialized) mpi_parallel::finalize();
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_parallel
