// -*- C++ -*-

#include "expic3d.hpp"

constexpr int order = 1;

class MainChunk;
class MainApplication;

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
      constexpr int Ns_mustbe = 2;

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
        int   nz = dims[0] + 2 * Nb;
        int   ny = dims[1] + 2 * Nb;
        int   nx = dims[2] + 2 * Nb;
        int   mp = nppc * dims[0] * dims[1] * dims[2];
        int64 id = static_cast<int64>(mp) * static_cast<int64>(this->myid);

        up.resize(Ns);

        // electron
        up[0]     = std::make_shared<Particle>(2 * mp, nz * ny * nx);
        up[0]->m  = mele;
        up[0]->q  = qele;
        up[0]->Np = mp;

        // ion
        up[1]     = std::make_shared<Particle>(2 * mp, nz * ny * nx);
        up[1]->m  = mion;
        up[1]->q  = qion;
        up[1]->Np = mp;

        // position
        for (int ip = 0; ip < mp; ip++) {
          float64* ele = &up[0]->xu(ip, 0);
          float64* ion = &up[1]->xu(ip, 0);
          float64  x   = uniform(mtp) * xlim[2] + xlim[0];
          float64  y   = uniform(mtp) * ylim[2] + ylim[0];
          float64  z   = uniform(mtp) * zlim[2] + zlim[0];

          // this guarantees charge neutrality
          ele[0] = x;
          ele[1] = y;
          ele[2] = z;
          ion[0] = x;
          ion[1] = y;
          ion[2] = z;

          // ID
          int64* ele_id64 = reinterpret_cast<int64*>(ele);
          int64* ion_id64 = reinterpret_cast<int64*>(ion);
          ele_id64[6]     = id + ip;
          ion_id64[6]     = id + ip;
        }

        // velocity
        for (int ip = 0; ip < mp; ip++) {
          float64* ele = &up[0]->xu(ip, 0);
          float64* ion = &up[1]->xu(ip, 0);

          // electrons
          ele[3] = normal(mtv) * vte;
          ele[4] = normal(mtv) * vte;
          ele[5] = normal(mtv) * vte;

          // incoming or reflected ions
          float64 r  = uniform(mtv);
          float64 vt = r > alpha ? vti : vtr;
          float64 vd = r > alpha ? vdi : vdr;

          ion[3] = normal(mtv) * vt + vd;
          ion[4] = normal(mtv) * vt;
          ion[5] = normal(mtv) * vt;
        }
      }

      // initial sort
      this->sort_particle(up);

      // allocate MPI buffer for particle
      int npmax = static_cast<int>(nppc * cc * delt / delh);
      this->set_mpi_buffer(mpibufvec[BoundaryParticle], 0, sizeof(int) * Ns,
                           Ns * npmax * sizeof(float64) * Particle::Nc);
    }
  }
};

class MainApplication : public ExPIC3D<order>
{
public:
  using ExPIC3D<order>::ExPIC3D; // inherit constructors

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
