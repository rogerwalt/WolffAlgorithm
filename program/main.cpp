#include <omp.h>        // openmp
#include <vector>
#include <cassert>
#include <random>
#include <cmath>
#include <iostream>
#include <stack>
#include <set>
#include <fstream>      // file operations
#include <map>
#include <chrono>

#include "myMeasurement.hpp"

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

class IsingLattice {
    private:
        const double J_ = 1;
        const double H_ = 0;
        const double T_;
        //const double k_ = 1.38065e-23; // joules/kelin
        const double k_ = 1;
        const double beta_;                                     // beta_ = 1/(k_*T_)
        const unsigned int latticeLength_;                      // length in all directions
        int spinSum_;                                           // used to track and calculate magnetization
        double energySum_;                                      // same for energy
        
        unsigned int steps_ = 0;                                // these two are used to calculate the acceptance rate for single spin flips
        unsigned int failed_ = 0;
        
        std::vector<int> latticeNode_;                          // 3-dim lattice
        std::mt19937 mt_;                                       // random generator
        std::uniform_real_distribution<double> real_d_;         // random number distributions
        std::uniform_int_distribution<unsigned int> int_i_x_;
        std::uniform_int_distribution<unsigned int> int_i_y_;
        std::uniform_int_distribution<unsigned int> int_i_z_;

        // flips the spin of the given node
        void flipSpin(unsigned int x, unsigned int y, unsigned int z) {
            assert(x<latticeLength_); assert(y<latticeLength_); assert(z<latticeLength_);
            
            // flip
            latticeNode_[x+y*latticeLength_+z*latticeLength_*latticeLength_] *= -1;
        }

    public:
        struct coordinate {
            unsigned int x,y,z;
            coordinate(unsigned int x, unsigned int y, unsigned int z) {
                this->x = x; this->y = y; this->z = z;
            }
            // operator used by set.count()
            bool operator<(const coordinate &other) const {
                if (x<other.x) {
                    return true;
                } else if (x==other.x && y<other.y) {
                    return true;
                } else if (y==other.y && x==other.x && z<other.z) {
                    return true;
                } else {
                    return false;
                }
            }
        };
        
        // constructor, initialize lattice_ with spin up
        IsingLattice(const double T, const unsigned int length)
        : T_(T)
        , latticeLength_(length)
        , beta_(1/(k_*T_))
        , latticeNode_(pow(latticeLength_, 3), 1)
        , real_d_(0., 1.)
        , int_i_x_(0, length-1)
        , int_i_y_(0, length-1)
        , int_i_z_(0, length-1)
        {
        }

        /*################*
         * System metrics *
         *################*/
        
        // returns the spin of the given node
        int getSpin(unsigned int x, unsigned int y, unsigned int z) {
            assert(x<latticeLength_); assert(y<latticeLength_); assert(z<latticeLength_);
            return latticeNode_[x+y*latticeLength_+z*latticeLength_*latticeLength_];
        }
        
        unsigned int getLatticeLength() {
            return latticeLength_;
        }

        // only nearest neighbours interact with each other
        // energy -J for each parallel pair, J for each antiparallel pair
        // energy -H resulting resulting from magnetic field
        double computeEnergyOfNode(unsigned int x, unsigned int y, unsigned int z) {
            double energy = 0;
            int localSpin = getSpin(x,y,z);

            // down
            energy += -J_*localSpin*getSpin(x,y,(z-1+latticeLength_)%latticeLength_);
            // up
            energy += -J_*localSpin*getSpin(x,y,(z+1)%latticeLength_);
            // left
            energy += -J_*localSpin*getSpin(x,(y-1+latticeLength_)%latticeLength_,z);
            // right
            energy += -J_*localSpin*getSpin(x,(y+1)%latticeLength_,z);
            // front
            energy += -J_*localSpin*getSpin((x-1+latticeLength_)%latticeLength_,y,z);
            // behind
            energy += -J_*localSpin*getSpin((x+1)%latticeLength_,y,z);

            // field energy
            energy += -H_*localSpin;

            return energy;
        }
        
        /*######################################*
         * Single spinflip metropolis functions *
         *######################################*/

        // do timestep (single flip metropolis)
        void timeStep() {
            double energyX, energyY, deltaE, prob;
            unsigned int x,y,z;

            // choose random node
            x=int_i_x_(mt_); y=int_i_y_(mt_); z=int_i_z_(mt_);

            // compute energy in configuration X (before flipping spin)
            energyX = computeEnergyOfNode(x,y,z);

            // new configuration Y has flipped spin of this node
            flipSpin(x,y,z);

            // compute energy in configuration Y (after flipping spin)
            energyY = computeEnergyOfNode(x,y,z);

            // compute energy difference of configurations Y and X
            deltaE = energyY - energyX;
            
            // if deltaE < 0 => spin stays flipped
            // if deltaE > 0 => probability that spin stays flipped is exp(-deltaE/k*T)
            if (deltaE > 0) {
                // calculate probability that spin stays flipped
                // probability is already normalized (or at least between 0 and 1)
                prob = exp(-deltaE / (k_*T_));

                // if the random double is lower than the probability, flip it back
                if (!(real_d_(mt_) < prob)) {
                    flipSpin(x,y,z);
                    ++failed_;
                }
            }
            ++steps_;
        }

        // do a time sweep (one sweep is N single flips, N = length³)
        void timeSweep() {
            for (unsigned int i=0; i<latticeNode_.size(); ++i) {
                timeStep();
            }
        }

        // do 3 sweeps
        void time3Sweep() {
            for (unsigned int i=0; i<3; ++i) {
                timeSweep();
            }
        }
        
        // does 5k*(number of nodes) single spin flip sweeps
        void do5kSweeps() {
            unsigned int numberOfSweeps = 5000 * latticeNode_.size();
            for (unsigned int i=0; i<numberOfSweeps; ++i) {
                timeSweep();
            }
        }
        
        /*#################*
         * Wolff functions *
         *#################*/
        
        void addNodeToClusterAndFlipSpinIfProbable(
                std::stack<coordinate> &stack,
                std::set<coordinate> &cluster,
                int startingNodeSpin,
                unsigned int x,
                unsigned int y,
                unsigned int z)
        {
            // if node already in cluster, dont check again
            if (cluster.count(coordinate(x,y,z))) return;
            
            // probability to connect nodes
            double prob = 1 - exp(2 * beta_ * J_ * startingNodeSpin * getSpin(x,y,z));
            
            if (real_d_(mt_) < prob) {
                stack.push(coordinate(x,y,z));
                cluster.insert(coordinate(x,y,z));
                flipSpin(x,y,z);
            }
        }
        
        unsigned int doWolffStep() {
            unsigned int x,y,z;
            std::stack<coordinate> stack;
            std::set<coordinate> cluster;
            
            // choose random node
            x=int_i_x_(mt_); y=int_i_y_(mt_); z=int_i_z_(mt_);
            
            // add current node coordinates to stack at position 0
            stack.push(coordinate(x,y,z));
            
            cluster.insert(coordinate(x,y,z));
            
            flipSpin(x,y,z);
            
            int startingNodeSpin = getSpin(x,y,z);
            
            // while stack not empty
            while (!stack.empty()) {
                // take last node coordinates from stack
                coordinate currentCoordinates = stack.top();
                stack.pop();
                
                // check all neighbours of node stack[i] with same spin state,
                  // add them to stack with probability p
                    
                // x+1
                addNodeToClusterAndFlipSpinIfProbable(stack, cluster, startingNodeSpin,
                        (currentCoordinates.x+1)%latticeLength_, currentCoordinates.y, currentCoordinates.z);

                // x-1
                addNodeToClusterAndFlipSpinIfProbable(stack, cluster, startingNodeSpin,
                        (currentCoordinates.x+latticeLength_-1)%latticeLength_, currentCoordinates.y, currentCoordinates.z);

                
                // y+1
                addNodeToClusterAndFlipSpinIfProbable(stack, cluster, startingNodeSpin,
                        currentCoordinates.x, (currentCoordinates.y+1)%latticeLength_, currentCoordinates.z);

                // y-1
                addNodeToClusterAndFlipSpinIfProbable(stack, cluster, startingNodeSpin,
                        currentCoordinates.x, (currentCoordinates.y+latticeLength_-1)%latticeLength_, currentCoordinates.z);

                
                // z+1
                addNodeToClusterAndFlipSpinIfProbable(stack, cluster, startingNodeSpin,
                        currentCoordinates.x, currentCoordinates.y, (currentCoordinates.z+1)%latticeLength_);

                // z-1
                addNodeToClusterAndFlipSpinIfProbable(stack, cluster, startingNodeSpin,
                        currentCoordinates.x, currentCoordinates.y, (currentCoordinates.z+latticeLength_-1)%latticeLength_);
            }
            
            return cluster.size();
        }
};

void measure() {
    unsigned int numWolffMeasurementValues = 1e4;
    unsigned int numSingleMeasurementValues = 5e5;
    
    const double startingTemperature = 0;
    const double temperatureStep = 0.5;
    const double endTemperature = 6;
    
    const double startingSystemSize = 5;
    const double systemSizeStep = 5;
    const double endSystemSize = 15;
    
    const unsigned int numberOfTemperatureLoops = (endTemperature - startingTemperature) / temperatureStep + 1;
    const unsigned int numberOfSystemSizeLoops = (endSystemSize - startingSystemSize) / systemSizeStep + 1;
    
    // open file
    std::ofstream myfile;
    myfile.open("tdata.py");
    
    // prepare file
    myfile << "data = [" << std::endl;

    #pragma omp parallel
    { // omp parallel
        // initialize measurement containers
        #pragma omp for schedule(dynamic)
        for (int i=0;i<numberOfTemperatureLoops;++i) { // temperature loop
            for (int j=0;j<numberOfSystemSizeLoops;++j) { // system size loop
                #pragma omp critical
                { std::cerr << "[" << omp_get_thread_num() << "] measuring for temperature: " << startingTemperature+i*temperatureStep << ", systemSize: " << startingSystemSize+j*systemSizeStep << std::endl; }

                // initialize the systems
                    IsingLattice wolffLattice(1, startingSystemSize+j*systemSizeStep);
                    IsingLattice singleLattice(1, startingSystemSize+j*systemSizeStep);

                // evolve in time and measure!
                auto start = std::chrono::steady_clock::now();
                unsigned int clusterSum = 0;
                for (unsigned k=0; k<numWolffMeasurementValues; ++k) {
                    clusterSum += wolffLattice.doWolffStep();
                    
                }
                auto end = std::chrono::steady_clock::now();

                int us = std::chrono::duration_cast<std::chrono::microseconds> (end-start).count();
                double timePerWolff = us/double(clusterSum)/double(numWolffMeasurementValues);

                #pragma omp critical
                { std::cerr << "[" << omp_get_thread_num() << "] wolff done" << std::endl; }

                
                start = std::chrono::steady_clock::now();
                for (unsigned k=0; k<numSingleMeasurementValues; ++k) {
                    singleLattice.timeStep();
                }
                end = std::chrono::steady_clock::now();

                us = std::chrono::duration_cast<std::chrono::microseconds> (end-start).count();         
                double timePerSingle = us/numSingleMeasurementValues;

                // output measurement data
                // write data
                #pragma omp critical
                {
                    std::cerr << "[" << omp_get_thread_num() << "] writing output" << std::endl;

                    myfile  << "\t{" << std::endl
                            << "\t\t" << "'temperature': " << startingTemperature+i*temperatureStep << "," << std::endl
                            << "\t\t" << "'systemSize': " << startingSystemSize+j*systemSizeStep << "," << std::endl
                            << "\t\t" << "'results': {" << std::endl
                    // wolff
                            << "\t\t\t" << "'wolff': {" << std::endl
                        // time
                                << "\t\t\t\t" << "'usPerflip': {" << std::endl
                                << timePerWolff
                                << "\t\t\t\t}," << std::endl
                            << "\t\t\t" << "}," << std::endl
                    // single spin flip
                            << "\t\t\t" << "'single': {" << std::endl
                        // time
                                << "\t\t\t\t" << "'usPerflip': {" << std::endl
                                << timePerSingle
                                << "\t\t\t\t}," << std::endl
                            << "\t\t\t" << "}" << std::endl
                    // end
                            << "\t\t}" << std::endl
                            << "\t}," << std::endl;
                }
            } // system size loop
        } // temperature loop
    } // omp parallel
    
    // finish file
    myfile << "]" << std::endl;
    
    // close file
    myfile.close();
}

// to plot data fast:
// ./ising_model | gnuplot -p -e 'plot "< cat /proc/$$/fd/0"'

// add 'gomp' to linker settings in codeblocks
// add '-fopenmp' to additional compiler options in netbeans
int main()
{
    measure();

    return 0;
}
