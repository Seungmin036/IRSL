//for ROS
#include <Controller/ControllerNode.hpp>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include <mujoco/mujoco.h>
#include "glfw_dispatch.h"
#include "simulate.h"
#include "array_safety.h"

#include <ros/ros.h>

#define MUJOCO_PLUGIN_DIR "mujoco_plugin"

#define CONTROLLER
// #define VISUALIZATION_PARTICLE

extern "C" {
#if defined(_WIN32) || defined(__CYGWIN__)
  #include <windows.h>
#else
  #if defined(__APPLE__)
    #include <mach-o/dyld.h>
  #endif
  #include <sys/errno.h>
  #include <unistd.h>
#endif
}

namespace {
namespace mj = ::mujoco;
namespace mju = ::mujoco::sample_util;

using ::mujoco::Glfw;

// constants
const double syncMisalign = 0.1;        // maximum mis-alignment before re-sync (simulation seconds)
const double simRefreshFraction = 0.7;  // fraction of refresh available for simulation
const int kErrorLength = 1024;          // load error string length

// model and data
mjModel* m = nullptr;
mjData* d = nullptr;

// control noise variables
mjtNum* ctrlnoise = nullptr;



//---------------------------------------- plugin handling -----------------------------------------

// return the path to the directory containing the current executable
// used to determine the location of auto-loaded plugin libraries
std::string getExecutableDir() {
#if defined(_WIN32) || defined(__CYGWIN__)
  constexpr char kPathSep = '\\';
  std::string realpath = [&]() -> std::string {
    std::unique_ptr<char[]> realpath(nullptr);
    DWORD buf_size = 128;
    bool success = false;
    while (!success) {
      realpath.reset(new(std::nothrow) char[buf_size]);
      if (!realpath) {
        std::cerr << "cannot allocate memory to store executable path\n";
        return "";
      }

      DWORD written = GetModuleFileNameA(nullptr, realpath.get(), buf_size);
      if (written < buf_size) {
        success = true;
      } else if (written == buf_size) {
        // realpath is too small, grow and retry
        buf_size *=2;
      } else {
        std::cerr << "failed to retrieve executable path: " << GetLastError() << "\n";
        return "";
      }
    }
    return realpath.get();
  }();
#else
  constexpr char kPathSep = '/';
#if defined(__APPLE__)
  std::unique_ptr<char[]> buf(nullptr);
  {
    std::uint32_t buf_size = 0;
    _NSGetExecutablePath(nullptr, &buf_size);
    buf.reset(new char[buf_size]);
    if (!buf) {
      std::cerr << "cannot allocate memory to store executable path\n";
      return "";
    }
    if (_NSGetExecutablePath(buf.get(), &buf_size)) {
      std::cerr << "unexpected error from _NSGetExecutablePath\n";
    }
  }
  const char* path = buf.get();
#else
  const char* path = "/proc/self/exe";
#endif
  std::string realpath = [&]() -> std::string {
    std::unique_ptr<char[]> realpath(nullptr);
    std::uint32_t buf_size = 128;
    bool success = false;
    while (!success) {
      realpath.reset(new(std::nothrow) char[buf_size]);
      if (!realpath) {
        std::cerr << "cannot allocate memory to store executable path\n";
        return "";
      }

      std::size_t written = readlink(path, realpath.get(), buf_size);
      if (written < buf_size) {
        realpath.get()[written] = '\0';
        success = true;
      } else if (written == -1) {
        if (errno == EINVAL) {
          // path is already not a symlink, just use it
          return path;
        }

        std::cerr << "error while resolving executable path: " << strerror(errno) << '\n';
        return "";
      } else {
        // realpath is too small, grow and retry
        buf_size *= 2;
      }
    }
    return realpath.get();
  }();
#endif

  if (realpath.empty()) {
    return "";
  }

  for (std::size_t i = realpath.size() - 1; i > 0; --i) {
    if (realpath.c_str()[i] == kPathSep) {
      return realpath.substr(0, i);
    }
  }

  // don't scan through the entire file system's root
  return "";
}



// scan for libraries in the plugin directory to load additional plugins
void scanPluginLibraries() {
  // check and print plugins that are linked directly into the executable
  int nplugin = mjp_pluginCount();
  if (nplugin) {
    std::printf("Built-in plugins:\n");
    for (int i = 0; i < nplugin; ++i) {
      std::printf("    %s\n", mjp_getPluginAtSlot(i)->name);
    }
  }

  // define platform-specific strings
#if defined(_WIN32) || defined(__CYGWIN__)
  const std::string sep = "\\";
#else
  const std::string sep = "/";
#endif


  // try to open the ${EXECDIR}/plugin directory
  // ${EXECDIR} is the directory containing the simulate binary itself
  const std::string executable_dir = getExecutableDir();
  if (executable_dir.empty()) {
    return;
  }

  const std::string plugin_dir = getExecutableDir() + sep + MUJOCO_PLUGIN_DIR;
  mj_loadAllPluginLibraries(
      plugin_dir.c_str(), +[](const char* filename, int first, int count) {
        std::printf("Plugins registered by library '%s':\n", filename);
        for (int i = first; i < first + count; ++i) {
          std::printf("    %s\n", mjp_getPluginAtSlot(i)->name);
        }
      });
}


//------------------------------------------- simulation -------------------------------------------


mjModel* LoadModel(const char* file, mj::Simulate& sim) {
  // this copy is needed so that the mju::strlen call below compiles
  char filename[mj::Simulate::kMaxFilenameLength];
  mju::strcpy_arr(filename, file);

  // make sure filename is not empty
  if (!filename[0]) {
    return nullptr;
  }

  // load and compile
  char loadError[kErrorLength] = "";
  mjModel* mnew = 0;
  if (mju::strlen_arr(filename)>4 &&
      !std::strncmp(filename + mju::strlen_arr(filename) - 4, ".mjb",
                    mju::sizeof_arr(filename) - mju::strlen_arr(filename)+4)) {
    mnew = mj_loadModel(filename, nullptr);
    if (!mnew) {
      mju::strcpy_arr(loadError, "could not load binary model");
    }
  } else {
    mnew = mj_loadXML(filename, nullptr, loadError, mj::Simulate::kMaxFilenameLength);
    // remove trailing newline character from loadError
    if (loadError[0]) {
      int error_length = mju::strlen_arr(loadError);
      if (loadError[error_length-1] == '\n') {
        loadError[error_length-1] = '\0';
      }
    }
  }

  mju::strcpy_arr(sim.loadError, loadError);

  if (!mnew) {
    std::printf("%s\n", loadError);
    return nullptr;
  }

  // compiler warning: print and pause
  if (loadError[0]) {
    // mj_forward() below will print the warning message
    std::printf("Model compiled, but simulation warning (paused):\n  %s\n", loadError);
    sim.run = 0;
  }
  return mnew;
}

// simulate in background thread (while rendering in main thread)
void PhysicsLoop(mj::Simulate& sim, const std::string & pinocchio_urdf_path) {
  // cpu-sim syncronization point
  double syncCPU = 0;
  mjtNum syncSim = 0;

  //for import pinocchio model
  #ifdef CONTROLLER
  ros::NodeHandle nh;

  pinocchio::Model pinocchio_model;

  pinocchio::urdf::buildModel(pinocchio_urdf_path,pinocchio_model);

  pinocchio_model.gravity.linear(pinocchio::Model::gravity981);


  ControllerNode controller_node(nh, pinocchio_model);
  controller_node.InitControllerNode(m,d);
  #endif


  // run until asked to exit
  while (!sim.exitrequest.load() && ros::ok()) {
    if (sim.droploadrequest.load()) {
      mjModel* mnew = LoadModel(sim.dropfilename, sim);
      sim.droploadrequest.store(false);

      mjData* dnew = nullptr;
      if (mnew) dnew = mj_makeData(mnew);
      if (dnew) {
        sim.load(sim.dropfilename, mnew, dnew);

        mj_deleteData(d);
        mj_deleteModel(m);

        m = mnew;
        d = dnew;
        mj_forward(m, d);

        // allocate ctrlnoise
        free(ctrlnoise);
        ctrlnoise = (mjtNum*) malloc(sizeof(mjtNum)*m->nu);
        mju_zero(ctrlnoise, m->nu);
      }
    }

    if (sim.uiloadrequest.load()) {
      sim.uiloadrequest.fetch_sub(1);
      mjModel* mnew = LoadModel(sim.filename, sim);
      mjData* dnew = nullptr;
      if (mnew) dnew = mj_makeData(mnew);
      if (dnew) {
        sim.load(sim.filename, mnew, dnew);

        mj_deleteData(d);
        mj_deleteModel(m);

        m = mnew;
        d = dnew;
        mj_forward(m, d);

        // allocate ctrlnoise
        free(ctrlnoise);
        ctrlnoise = static_cast<mjtNum*>(malloc(sizeof(mjtNum)*m->nu));
        mju_zero(ctrlnoise, m->nu);
      }
    }
    
    // sleep for 1 ms or yield, to let main thread run
    //  yield results in busy wait - which has better timing but kills battery life
    if (sim.run && sim.busywait) {
      std::this_thread::yield();
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    {
      // lock the sim mutex
      const std::lock_guard<std::mutex> lock(sim.mtx);

      // run only if model is present
      if (m) {
        // running
        if (sim.run) {
          // record cpu time at start of iteration
          double startCPU = Glfw().glfwGetTime();

          // elapsed CPU and simulation time since last sync
          double elapsedCPU = startCPU - syncCPU;
          double elapsedSim = d->time - syncSim;

          // inject noise
          if (sim.ctrlnoisestd) {
            // convert rate and scale to discrete time (Ornstein–Uhlenbeck)
            mjtNum rate = mju_exp(-m->opt.timestep / sim.ctrlnoiserate);
            mjtNum scale = sim.ctrlnoisestd * mju_sqrt(1-rate*rate);

            for (int i=0; i<m->nu; i++) {
              // update noise
              ctrlnoise[i] = rate * ctrlnoise[i] + scale * mju_standardNormal(nullptr);

              // apply noise
              d->ctrl[i] = ctrlnoise[i];
            }
          }

          // requested slow-down factor
          double slowdown = 100 / sim.percentRealTime[sim.realTimeIndex];

          // misalignment condition: distance from target sim time is bigger than syncmisalign
          bool misaligned = mju_abs(elapsedCPU/slowdown - elapsedSim) > syncMisalign;

          // Run controller
        //   robot.init(d, m);
        //   robot.updateData(d, m);

          // out-of-sync (for any reason): reset sync times, step
          if (elapsedSim < 0 || elapsedCPU < 0 || syncCPU == 0 || misaligned || sim.speedChanged) {
            // re-sync
            syncCPU = startCPU;
            syncSim = d->time;
            sim.speedChanged = false;

            // clear old perturbations, apply new
            mju_zero(d->xfrc_applied, 6*m->nbody);
            sim.applyposepertubations(0);  // move mocap bodies only
            sim.applyforceperturbations();


            #ifdef CONTROLLER
            controller_node.Control_Loop(m,d);
            ros::spinOnce();
            #endif

            // run single step, let next iteration deal with timing
            mj_step(m, d);
          }

          // in-sync: step until ahead of cpu
          else {
            bool measured = false;
            mjtNum prevSim = d->time;
            double refreshTime = simRefreshFraction/sim.refreshRate;

            // step while sim lags behind cpu and within refreshTime
            while ((d->time - syncSim)*slowdown < (Glfw().glfwGetTime()-syncCPU) &&
                   (Glfw().glfwGetTime()-startCPU) < refreshTime) {
              // measure slowdown before first step
              if (!measured && elapsedSim) {
                sim.measuredSlowdown = elapsedCPU / elapsedSim;
                measured = true;
              }

              // clear old perturbations, apply new
              mju_zero(d->xfrc_applied, 6*m->nbody);
              sim.applyposepertubations(0);  // move mocap bodies only
              sim.applyforceperturbations();

              #ifdef CONTROLLER
              controller_node.Control_Loop(m,d);
              ros::spinOnce();
              #endif

              // call mj_step
              mj_step(m, d);

              // break if reset
              if (d->time < prevSim) {
                break;
              }
            }
          }
        }

        // paused
        else {
          // apply pose perturbation
          sim.applyposepertubations(1);  // move mocap and dynamic bodies

          // run mj_forward, to update rendering and joint sliders
          mj_forward(m, d);
          ros::spinOnce();
        }
      }
    }  // release std::lock_guard<std::mutex>
  }
}
}  // namespace

//-------------------------------------- physics_thread --------------------------------------------

void PhysicsThread(mj::Simulate* sim, const char* mujoco_xml_path, const std::string &pinocchio_urdf_path) {
  // request loadmodel if file given (otherwise drag-and-drop)
  if (mujoco_xml_path != nullptr) {
    m = LoadModel(mujoco_xml_path, *sim);
    if (m) d = mj_makeData(m);
    if (d) {
      sim->load(mujoco_xml_path, m, d);
      mj_forward(m, d);
      // allocate ctrlnoise
      free(ctrlnoise);
      ctrlnoise = static_cast<mjtNum*>(malloc(sizeof(mjtNum)*m->nu));
      mju_zero(ctrlnoise, m->nu);
    }
  }

  PhysicsLoop(*sim, pinocchio_urdf_path);

  // delete everything we allocated
  free(ctrlnoise);
  mj_deleteData(d);
  mj_deleteModel(m);
}

//------------------------------------------ main --------------------------------------------------

// machinery for replacing command line error by a macOS dialog box when running under Rosetta
#if defined(__APPLE__) && defined(__AVX__)
extern void DisplayErrorDialogBox(const char* title, const char* msg);
static const char* rosetta_error_msg = nullptr;
__attribute__((used, visibility("default"))) extern "C" void _mj_rosettaError(const char* msg) {
  rosetta_error_msg = msg;
}
#endif

// run event loop
int main(int argc, char** argv) 
{
    ros::init(argc, argv, "mujoco_kinova");

    // print version, check compatibility
    std::printf("MuJoCo version %s\n", mj_versionString());
    if (mjVERSION_HEADER!=mj_version()) {
        mju_error("Headers and library have different versions");
    }
    // scan for libraries in the plugin directory to load additional plugins
    scanPluginLibraries();

    // simulate object encapsulates the UI
    auto sim = std::make_unique<mj::Simulate>();

    // init GLFW
    if (!Glfw().glfwInit()) {
        mju_error("could not initialize GLFW");
    }

    std::string pinocchio_urdf_path;
    char* mujoco_xml_path = nullptr;

    pinocchio_urdf_path = "/home/irsl/catkin_ws/src/Mujoco-1dof-main/model/kinova/pinocchio_model/pinocchio_urdf/kinova.urdf";
    mujoco_xml_path = "/home/irsl/catkin_ws/src/Mujoco-1dof-main/model/kinova/mujoco_model/kinova_gen3.xml";


    // start physics thread
    std::thread physicsthreadhandle = std::thread(&PhysicsThread, sim.get(), mujoco_xml_path, pinocchio_urdf_path);

    // start simulation UI loop (blocking call)
    sim->renderloop();

    physicsthreadhandle.join();

    // terminate GLFW (crashes with Linux NVidia drivers)
    #if defined(__APPLE__) || defined(_WIN32)
    Glfw().glfwTerminate();
    #endif

    return 0;
}