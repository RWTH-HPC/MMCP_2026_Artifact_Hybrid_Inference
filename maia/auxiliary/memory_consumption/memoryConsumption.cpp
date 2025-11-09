#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <mpi.h>

#include <chrono>
#include <thread>

using namespace std;

/// \brief Write memory statistics
void writeMemoryStatistics(const MPI_Comm comm, const int& globalDomainId,
                                                                     const int& globalNoDomains) {

  // OBTAIN MEMORY USAGE
  int physMem = 0;
  int virtMem = 0;
  {
    // Open status file of the process
    ifstream fin;
    fin.open("/proc/self/status");
    // Terminate if file is not found
    if (!fin) {
      cout << "Error while opening file!" << endl;
      exit(1);
    }

    string name;
    int buffer;
    bool foundVmrss = false;
    bool foundVmdata = false;
    bool foundVmstk = false;

    // Read all lines and get memory statistics
    while (!fin.eof()) {
      buffer = 0;
      fin >> name;
      // Physical memory currently used by current process
      if (name == "VmRSS:") {
        fin >> buffer;
        physMem = buffer;
        foundVmrss = true;
      } else if (name == "VmData:") {
        fin >> buffer;
        virtMem += buffer;
        foundVmdata = true;
      } else if (name == "VmStk:") {
        fin >> buffer;
        virtMem += buffer;
        foundVmstk = true;
      }
    }
    fin.close();

    // Terminate if memory usage couldn't be gather from file
    if (!foundVmrss | !foundVmdata | !foundVmstk) {
      cout << "Could not gather memory usage!" << endl;
      exit(1);
    }
  }

  vector<int> physMemPerProcess(globalNoDomains, 0);
  vector<int> virtMemPerProcess(globalNoDomains, 0);

  // Gather memory from each process
  MPI_Gather(&physMem, 1, MPI_INT, &physMemPerProcess[0], 1, MPI_INT, 0, comm);
  MPI_Gather(&virtMem, 1, MPI_INT, &virtMemPerProcess[0], 1, MPI_INT, 0, comm);

  // Get global memory statistics
  if (globalDomainId == 0) {
    int totalPhysMem =  0;
    int totalVirtMem =  0;
    for (int i = 0; i < globalNoDomains; i++) {
      totalPhysMem += physMemPerProcess[i];
      totalVirtMem += virtMemPerProcess[i];
    }

    // Write memory statistics
    // Memory per process
    cout << "\n";
    for (int i = 0; i < globalNoDomains; i++) {
      cout << " Process " << i << " - Current physical memory usage = "
           << (float)physMemPerProcess[i] / 1024.0
           << " MB (allocation size = " << (float)virtMemPerProcess[i] / 1024.0 << " MB)"
           << std::endl;
    }

    // Average memory
    cout << " Average physical memory usage: "
         << (float)totalPhysMem/(globalNoDomains*1024) << " MB\n";

    // Min/Max memory
    cout << " Minimun physical memory usage: "
         << (float)*min_element(physMemPerProcess.begin(), physMemPerProcess.end())/1024
         << " MB\n Maximum physical memory usage: "
         << (float)*max_element(physMemPerProcess.begin(), physMemPerProcess.end())/1024
         << " MB\n\n";

    // Total memory
    cout << " Total physical memory usage (RAM): "
         << (float)totalPhysMem/1024 << " MB\n";
    cout << " Total allocation size (Virtual memory): "
         << (float)totalVirtMem/1024 << " MB\n" << endl;

  }
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Pass base memory to be allocated on the heap (in MB) as argument." << std::endl;
    return 0;
  }

  MPI_Init(&argc, &argv);

  int domainId, noDomains;
  MPI_Comm_rank(MPI_COMM_WORLD, &domainId);
  MPI_Comm_size(MPI_COMM_WORLD, &noDomains);

  // Allocate memory, increase allocation size with domain id
  const long mem = atoi(argv[1]);
  const long memSize = (mem + domainId*50) *1024*1024;
  if (domainId == 0) {
    std::cout << "Base allocation size: " << mem << "; memSize=" << memSize << std::endl;
  }

  char *bufferHeap = new char[memSize];

  // Memory reported only in "allocation size" (not in RAM)
  writeMemoryStatistics(MPI_COMM_WORLD, domainId, noDomains);

  for (int i=0; i < memSize; i++) {
    bufferHeap[i] = 'a';
  }

  // Memory now reported in "Current physical memory usage"
  writeMemoryStatistics(MPI_COMM_WORLD, domainId, noDomains);

  // Sleep for 30 seconds (to validate memory check with top/...)
  std::this_thread::sleep_for(std::chrono::seconds(30));

  delete [] bufferHeap;

  MPI_Finalize();
}
