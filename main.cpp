// -*- C++ -*-

#include "diagnoser.hpp"
#include "expic3d.hpp"

constexpr int order = 1;

class MainChunk;
class MainApplication;
using MainDiagnoser = Diagnoser;

class MainChunk : public ExChunk3D<order>
{
public:
  using ExChunk3D<order>::ExChunk3D; // inherit constructors

  virtual void setup(json& config) override
  {
    // parameter for load balancing
    field_load = config.value("field_load", 1.0);

    // check validity of assumptions
    {
      constexpr int Ns_mustbe = 3;

      Ns = config["Ns"].get<int>();

      if (Ns != Ns_mustbe) {
        ERROR << "Assumption of Ns = 2 is violated";
        exit(-1);
      }
    }

    // speed of light
    cc = config["cc"].get<float64>();

    int     nppc  = config["nppc"].get<int>();
    float64 wp    = config["wp"].get<float64>();
    float64 delt  = config["delt"].get<float64>();
    float64 delh  = config["delh"].get<float64>();
    float64 mime  = config["mime"].get<float64>();
    float64 mach  = config["mach"].get<float64>();
    float64 theta = config["theta"].get<float64>();
    float64 phi   = config["phi"].get<float64>();
    float64 sigma = config["sigma"].get<float64>();
    float64 alpha = config["alpha"].get<float64>();
    float64 betae = config["betae"].get<float64>();
    float64 betai = config["betai"].get<float64>();
    float64 betar = config["betar"].get<float64>();
    float64 mele  = 1.0 / (sigma * nppc);
    float64 qele  = -wp * sqrt(sigma) * mele;
    float64 mion  = mele * mime;
    float64 qion  = -qele;
    float64 b0    = cc * sqrt(sigma) / std::abs(qele / mele);
    float64 vae   = cc * sqrt(sigma);
    float64 vai   = cc * sqrt(sigma / mime);
    float64 vte   = vae * sqrt(0.5 * betae);
    float64 vti   = vai * sqrt(0.5 * betai);
    float64 vtr   = vai * sqrt(0.5 * betar);
    float64 vdi   = -vai * mach * alpha;
    float64 vdr   = +vai * mach * (1 - alpha);

    // set grid size and coordinate
    set_coordinate(delh, delh, delh);

    //
    // initialize field
    //
    {
      float64 Bx = b0 * cos(theta / 180 * nix::math::pi);
      float64 By = b0 * sin(theta / 180 * nix::math::pi) * cos(phi / 180 * nix::math::pi);
      float64 Bz = b0 * sin(theta / 180 * nix::math::pi) * sin(phi / 180 * nix::math::pi);

      // memory allocation
      allocate();

      for (int iz = Lbz; iz <= Ubz; iz++) {
        for (int iy = Lby; iy <= Uby; iy++) {
          for (int ix = Lbx; ix <= Ubx; ix++) {
            uf(iz, iy, ix, 0) = 0;
            uf(iz, iy, ix, 1) = 0;
            uf(iz, iy, ix, 2) = 0;
            uf(iz, iy, ix, 3) = Bx;
            uf(iz, iy, ix, 4) = By;
            uf(iz, iy, ix, 5) = Bz;
          }
        }
      }

      // allocate MPI buffer for field
      this->set_mpi_buffer(mpibufvec[BoundaryEmf], 0, 0, sizeof(float64) * 6);
      this->set_mpi_buffer(mpibufvec[BoundaryCur], 0, 0, sizeof(float64) * 4);
      this->set_mpi_buffer(mpibufvec[BoundaryMom], 0, 0, sizeof(float64) * Ns * 11);
    }

    //
    // initialize particles
    //
    {
      // random number generators
      int                                     random_seed = 0;
      std::mt19937                            mtp(0);
      std::mt19937                            mtv(0);
      std::uniform_real_distribution<float64> uniform(0.0, 1.0);
      std::normal_distribution<float64>       normal(0.0, 1.0);

      // random seed
      {
        std::string seed_type = config.value("seed_type", "random"); // random by default

        if (seed_type == "random") {
          random_seed = std::random_device()();
        } else if (seed_type == "chunkid") {
          random_seed = this->myid; // chunk ID
        } else {
          ERROR << tfm::format("Ignoring invalid seed_type: %s", seed_type);
        }

        mtp.seed(random_seed);
        mtv.seed(random_seed);
      }

      {
        int   nz  = dims[0] + 2 * Nb;
        int   ny  = dims[1] + 2 * Nb;
        int   nx  = dims[2] + 2 * Nb;
        int   mp  = nppc * dims[0] * dims[1] * dims[2];
        int   mp1 = mp * (1 - alpha);
        int   mp2 = mp - mp1;
        int64 id  = static_cast<int64>(mp) * static_cast<int64>(this->myid);

        up.resize(Ns);

        // electron
        up[0]     = std::make_shared<Particle>(2 * mp, nz * ny * nx);
        up[0]->m  = mele;
        up[0]->q  = qele;
        up[0]->Np = mp;

        // incoming ion
        up[1]     = std::make_shared<Particle>(2 * mp1, nz * ny * nx);
        up[1]->m  = mion;
        up[1]->q  = qion;
        up[1]->Np = mp1;

        // reflected ion
        up[2]     = std::make_shared<Particle>(2 * mp2, nz * ny * nx);
        up[2]->m  = mion;
        up[2]->q  = qion;
        up[2]->Np = mp2;

        // initialize particle distribution
        std::vector<int>     mp_ele{0, mp1};
        std::vector<int>     mp_ion{mp1, mp2};
        std::vector<float64> vt_ion{vti, vtr};
        std::vector<float64> vd_ion{vdi, vdr};

        for (int is = 0; is < 2; is++) {
          const int is_ele      = 0;
          const int is_ion      = is + 1;
          const int ip_ele_zero = mp_ele[is];
          const int ip_ion_zero = 0;

          for (int ip = 0; ip < mp_ion[is]; ip++) {
            const int ip_ele = ip + ip_ele_zero;
            const int ip_ion = ip + ip_ion_zero;

            // position: using these guarantees charge neutrality
            float64 x = uniform(mtp) * xlim[2] + xlim[0];
            float64 y = uniform(mtp) * ylim[2] + ylim[0];
            float64 z = uniform(mtp) * zlim[2] + zlim[0];

            // electrons
            up[is_ele]->xu(ip_ele, 0) = x;
            up[is_ele]->xu(ip_ele, 1) = y;
            up[is_ele]->xu(ip_ele, 2) = z;
            up[is_ele]->xu(ip_ele, 3) = normal(mtv) * vte;
            up[is_ele]->xu(ip_ele, 4) = normal(mtv) * vte;
            up[is_ele]->xu(ip_ele, 5) = normal(mtv) * vte;

            // ions
            up[is_ion]->xu(ip_ion, 0) = x;
            up[is_ion]->xu(ip_ion, 1) = y;
            up[is_ion]->xu(ip_ion, 2) = z;
            up[is_ion]->xu(ip_ion, 3) = normal(mtv) * vt_ion[is] + vd_ion[is];
            up[is_ion]->xu(ip_ion, 4) = normal(mtv) * vt_ion[is];
            up[is_ion]->xu(ip_ion, 5) = normal(mtv) * vt_ion[is];

            // ID
            int64* ele_id64 = reinterpret_cast<int64*>(&up[is_ele]->xu(ip_ele, 0));
            int64* ion_id64 = reinterpret_cast<int64*>(&up[is_ion]->xu(ip_ion, 0));
            ele_id64[6]     = id + ip_ele;
            ion_id64[6]     = id + ip_ele;
          }
        }
      }

      // initial sort
      this->sort_particle(up);

      // use default MPI buffer allocator for particle
      float64 fraction = config.value("mpi_buffer_fraction", cc * delt / delh);
      setup_particle_mpi_buffer(fraction);
    }
  }
};

class MainApplication : public ExPIC3D<order, MainDiagnoser>
{
public:
  using ExPIC3D<order, MainDiagnoser>::ExPIC3D; // inherit constructors

  std::unique_ptr<ExChunk3D<order>> create_chunk(const int dims[], int id) override
  {
    return std::make_unique<MainChunk>(dims, id);
  }
};

//
// main
//
int main(int argc, char** argv)
{
  MainApplication app(argc, argv);
  return app.main(std::cout);
}

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
