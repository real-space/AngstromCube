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
#include "green_memory.hxx" // get_memory
#include "recorded_warnings.hxx" // warn

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
      if (echo > 0) std::printf("# MPI exchange of potential, MPI one-sided communication, Noco=%d\n", Noco);
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
                  if (echo > 9 + 5*spin) { std::printf("# exchange: rank#%i get potential at cube id o%llo from rank#%i element %i\n", me, global_id, owner, iloc); std::fflush(stdout); }
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
      status += MPI_Barrier(comm); std::printf("# rank#%d hit MPI_Barrier at line %d\n", me, __LINE__);
      status += MPI_Win_fence(assertions, window);
      // MPI_Win_free is only needed if we use MPI_Win_alloc above
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == Noco*Noco*nrows);

      return status;
  } // potential_exchange


  // For MPI parallel calculations the atom matrix values need to be exchanged
  status_t dyadic_exchange(                                    // here, nSHO is the global maximum of nsho[ia]
        double       *const mat_out // output effective atom matrices, data layout mat_out[nrows*Noco^2*2*nSHO^2]
      , std::vector<int64_t> const & requests  // indices requested by this MPI process,  [nrows]
      , double const *const mat_inp //  input effective atom matrices, data layout mat_inp[ncols*Noco^2*2*nSHO^2]
      , std::vector<int64_t> const & offerings // indices offered   by this MPI process,  [ncols]
      , rank_int_t const owner_rank[] // where to find it, [nall]
      , uint32_t const nall // number of all atoms
      , int const count
      , bool const debug // =true
      , int const echo // =0, log-level
  ) {
      // The number of local atoms is limited to 2^16 == 65536
      if (echo > 0) std::printf("# MPI exchange of atom matrices, MPI one-sided communication, packages of %.3f k doubles\n", count*.001);

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
              if (echo > 5) std::printf("# exchange: rank #%i local copy of matrices of atom#%i\n", me, global_id);
              set(&mat_out[row*count], count, &mat_inp[iloc*count]);
              ++stats[0];
          } else { // me == owner
              ++stats[1];
#ifndef   HAS_NO_MPI
              if (echo > 7) std::printf("# exchange: rank #%i get matrices of atom#%i from rank #%i element %i\n", me, global_id, owner, iloc);
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
      status += MPI_Win_fence(assertions, window);
      // MPI_Win_free is only needed if we use MPI_Win_alloc above
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == nrows);

      return status;
  } // dyadic_exchange





    // New Grouping: requests are first sorted into a RequestList_t and then, data exchange happens

    RequestList_t::RequestList_t( // constructor implementation
        std::vector<int64_t> const & requests
      , std::vector<int64_t> const & offerings
      , rank_int_t const owner_rank[] // where to find it, [nb[Z]*nb[Y]*nb[X]]
      , uint32_t const nb[3] // global bounding box or {natoms,0,0}
      , int const echo // =0 // log-level
    ) {
        int constexpr X=0, Y=1, Z=2;
        auto const grid = size_t(nb[Z])*size_t(nb[Y])*size_t(nb[X]);
        auto const nall = grid ? grid : nb[X] + nb[Y] + nb[Z];
        if (echo > 7) std::printf("# RequestList_t %d %d %d, nall= %ld\n", nb[X], nb[Y], nb[Z], nall);
        auto const ncols = offerings.size(); // number of offerings
        auto const nrows =  requests.size(); // number of requests

#ifndef   HAS_NO_MPI
        int constexpr debug = 1;
        auto const comm = MPI_COMM_WORLD;
        auto const me = mpi_parallel::rank(comm);

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
        owner = std::vector<int32_t>(nrows, 0);
        index = std::vector<uint32_t>(nrows, 0);
        requested_id = std::vector<int64_t>(nrows, 0);
        window_size = ncols;

        for (size_t row = 0; row < nrows; ++row) {
            auto const global_id = requests[row];

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

            owner[row]      =  owner_rank[iall];
            auto const iloc = local_index[iall];

#else  // HAS_NO_MPI

            // without MPI search global_id in offerings
            int64_t iloc{-1};
            for (size_t col = 0; col < ncols && iloc < 0; ++col) {
                if (global_id == offerings[col]) iloc = col;
            } // col
            assert(iloc > -1 && "index not found in local offerings");

#endif // HAS_NO_MPI

            index[row] = iloc;
            requested_id[row] = global_id; // copy the array of requests
        } // row

    } // constructor implementation





  status_t potential_exchange(
        double    (*const Veff[4])[64]  // output effective potentials,  data layout Veff[Noco^2][nrows][64]
      , double const (*const Vinp)[64]  //  input effective potentials,  data layout Vinp[ncols*Noco^2 ][64]
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

      status_t status(0);

      auto const comm = MPI_COMM_WORLD;
      auto const me = mpi_parallel::rank(comm);
      auto const nreq = requests.size(); // number of requests
      size_t stats[] = {0, 0};

#ifndef   HAS_NO_MPI

      // set up a memory window to read from
      MPI_Win window;
      uint const disp_unit = Noco*Noco*64*sizeof(double); // in Bytes
      size_t const win_size = requests.window()*disp_unit; // in Bytes
      status += MPI_Win_create((void*)Vinp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);

      // synchronize processes
      int const assertions = 0; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
      status += MPI_Win_fence(assertions, window);

#endif // HAS_NO_MPI

      for (size_t ireq = 0; ireq < nreq; ++ireq) {
          auto const global_id = requests.requested_id[ireq];
          auto const owner     = requests.owner[ireq];
          auto const iloc      = requests.index[ireq];

          auto const target_disp = iloc*Noco*Noco;
          if (me == owner) {
              if (echo > 7) std::printf("# exchange: rank #%i local copy of potential at cube id o%llo\n", me, global_id);
              for (int spin = 0; spin < Noco*Noco; ++spin) {
                  set(Veff[spin][ireq], 64, Vinp[target_disp + spin]);
              } // spin
              ++stats[0];
          } else { // me == owner
              ++stats[1];
#ifndef   HAS_NO_MPI
              if (echo > 9) std::printf("# exchange: rank #%i get potential at cube id o%llo from rank #%i element %i\n", me, global_id, owner, iloc);
              std::vector<double> Vget(Noco*Noco*64);
              status += MPI_Get(Vget.data(), Noco*Noco*64, MPI_DOUBLE, owner, iloc, Noco*Noco*64, MPI_DOUBLE, window);
              for (int spin = 0; spin < Noco*Noco; ++spin) {
                  set(Veff[spin][ireq], 64, &Vget[spin*64]);
              } // spin
#else  // HAS_NO_MPI
              error("Without MPI all potential elements must reside in the same process, me=%i, owner=%i", me, owner);
#endif // HAS_NO_MPI
          } // me == owner
      } // ireq

#ifndef   HAS_NO_MPI
      // synchronize processes
      status += MPI_Win_fence(assertions, window);
      // MPI_Win_free is only needed if we use MPI_Win_alloc above
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == nreq);

      return status;
  } // potential_exchange with RequestList_t



  status_t dyadic_exchange(                                    // here, nSHO is the global maximum of nsho[ia]
        double       *const mat_out // output effective atom matrices, data layout mat_out[nrows*Noco^2*2*nSHO^2]
      , double const *const mat_inp //  input effective atom matrices, data layout mat_inp[ncols*Noco^2*2*nSHO^2]
      , RequestList_t const & requests
      , int const count
      , int const echo // =0, log-level
  ) {
      // The number of local atoms is limited to 2^16 == 65536
      if (echo > 0) std::printf("# new MPI exchange of atom matrices, MPI one-sided communication, packages of %.3f k doubles\n", count*.001);

      assert(mat_out && "may not be called with a nullptr for output");
      assert(mat_inp && "may not be called with a nullptr for input");

      status_t status(0);

      auto const comm = MPI_COMM_WORLD;
      auto const me = mpi_parallel::rank(comm);
      auto const nreq = requests.size(); // number of requests
      size_t stats[] = {0, 0};

#ifndef   HAS_NO_MPI

      // set up a memory window to read from
      MPI_Win window;
      uint const disp_unit = count*sizeof(double); // in Bytes
      size_t const win_size  = requests.window()*disp_unit; // in Bytes
      status += MPI_Win_create((void*)mat_inp, win_size, disp_unit, MPI_INFO_NULL, comm, &window);

      // synchronize processes
      int const assertions = 0; // use bitwise or, e.g. MPI_MODE_NOSTORE | MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE | MPI_MODE_NOSUCCEED;
      status += MPI_Win_fence(assertions, window);

#endif // HAS_NO_MPI

      for (size_t ireq = 0; ireq < nreq; ++ireq) {
          auto const global_id = requests.requested_id[ireq];
          auto const owner     = requests.owner[ireq];
          auto const iloc      = requests.index[ireq];

          if (me == owner) {
              if (echo > 5) std::printf("# exchange: rank #%i local copy of matrices of atom#%i\n", me, global_id);
              set(&mat_out[ireq*count], count, &mat_inp[iloc*count]);
              ++stats[0];
          } else { // me == owner
              ++stats[1];
#ifndef   HAS_NO_MPI
              if (echo > 7) std::printf("# exchange: rank #%i get matrices of atom#%i from rank #%i element %i\n", me, global_id, owner, iloc);
              status += MPI_Get(&mat_out[ireq*count], count, MPI_DOUBLE, owner, iloc, count, MPI_DOUBLE, window);
#else  // HAS_NO_MPI
              error("Without MPI all atom matrices must reside in the same process, me=%i, owner=%i", me, owner);
#endif // HAS_NO_MPI
          } // me == owner
      } // ireq

#ifndef   HAS_NO_MPI
      // synchronize processes
      status += MPI_Win_fence(assertions, window);
      // MPI_Win_free is only needed if we use MPI_Win_alloc above
#endif // HAS_NO_MPI

      if (echo > 7) std::printf("# rank #%i copied %.3f k and pulled %.3f k elements\n", me, stats[0]*.001, stats[1]*.001);
      assert(stats[0] + stats[1] == nreq);

      return status;
  } // dyadic_exchange with RequestList_t





} // namespace green_parallel

#include "status.hxx" // status_t

namespace green_parallel {

#ifdef    NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else  // NO_UNIT_TESTS

  status_t test_exchanges(int echo=0, int const nSHO=220) {
      status_t stat(0);
      uint32_t const nb[] = {2, 2, 1};
      std::vector<int64_t> requests = {0, 3, 2, 1};
      auto const nrows = requests.size();
      int const me = mpi_parallel::rank(),
                np = mpi_parallel::size(); assert(np > 0);

      std::vector<int64_t> offerings(0);
      std::vector<uint16_t> owner_rank(nrows, 0);
      for(int row{0}; row < nrows; ++row) {
          int const rank = row % np; // block-cyclic distribution
          owner_rank[row] = rank;
          if (me == rank) offerings.push_back(row);
      } // row
      RequestList_t rlV(requests, offerings, owner_rank.data(), nb, echo); // test_new_structure
      uint32_t const nall[] = {nb[2]*nb[1]*nb[0], 0, 0};
      RequestList_t rlD(requests, offerings, owner_rank.data(), nall, echo); // test_new_structure

      double (*pot_out[2*2])[64];
      for (int spin = 0; spin < 2*2; ++spin) pot_out[spin] = get_memory<double[64]>(nrows);
      for (int Noco = 1; Noco <= 2; ++Noco) {
          auto const pot_inp = new double[nrows*Noco*Noco][64];
          for (int row = 0; row < nrows; ++row) pot_inp[row*Noco*Noco][0] = 0.5 + me;
          stat += green_parallel::potential_exchange(pot_out, requests, pot_inp, offerings, owner_rank.data(), nb, Noco, true, echo);
          stat += green_parallel::potential_exchange(pot_out, pot_inp, rlV, Noco, echo);
          delete[] pot_inp;
          for (int row = 0; row < nrows; ++row) stat += (pot_out[0][row][0] != 0.5 + row % np);

          int const count = Noco*Noco*2*nSHO*nSHO;
          std::vector<double> mat_out(nrows*count), mat_inp(nrows*count);
          stat += green_parallel::dyadic_exchange(mat_out.data(), requests, mat_inp.data(), offerings, owner_rank.data(), nall[0], count, true, echo);
          stat += green_parallel::dyadic_exchange(mat_out.data(), mat_inp.data(), rlD, count, echo);
      } // Noco

      for (int spin = 0; spin < 2*2; ++spin) free_memory(pot_out[spin]);
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
