#include <omp.h>        // openmp
#include <vector>
#include <cassert>
#include <random>
#include <cmath>
#include <iostream>
#include <algorithm>    // std::transform


class IsingLattice {
    private:
        const double J_ = 1;
        const double H_ = 0;
        const double T_;
        //const double k_ = 1.38065e-23; // joules/kelin
        const double k_ = 1;
        const double beta_;                                     // beta_ = 1/(k_*T_)
        const unsigned int latticeLength_;                      // length in all directions
        std::vector<int> latticeNode_;                          // 3-dim lattice
        std::mt19937 mt_;                                       // random generator
        std::uniform_real_distribution<double> real_d_;         // random number distributions
        std::uniform_int_distribution<unsigned int> int_i_x_;
        std::uniform_int_distribution<unsigned int> int_i_y_;
        std::uniform_int_distribution<unsigned int> int_i_z_;

        // flips the spin of the given node
        void flipSpin(unsigned int x, unsigned int y, unsigned int z) {
            assert(x<latticeLength_); assert(y<latticeLength_); assert(z<latticeLength_);
            latticeNode_[x+y*latticeLength_+z*latticeLength_*latticeLength_] *= -1;
        }

    public:
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
        { }

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

        // sum up the eneriges of all nodes
        double computeEnergyOfSystem() {
            double energy = 0;

            // loop over all nodes
            for (unsigned int x=0; x<latticeLength_; ++x) {
                for (unsigned int y=0; y<latticeLength_; ++y) {
                    for (unsigned int z=0; z<latticeLength_; ++z) {
                        energy += computeEnergyOfNode(x,y,z);
                    }
                }
            }

            return energy;
        }

        double computeMagnetization() {
            double sum = 0;
            // loop over all nodes
            for (unsigned int x=0; x<latticeLength_; ++x) {
                for (unsigned int y=0; y<latticeLength_; ++y) {
                    for (unsigned int z=0; z<latticeLength_; ++z) {
                        sum += getSpin(x,y,z);
                    }
                }
            }
            return 1/pow(latticeLength_, 3) * sum;
        }

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
                }
            }
        }

        // do a time sweep (one sweep is N single flips, N = length³)
        void timeSweep() {
            for (unsigned int i=0; i<pow(latticeLength_, 3); ++i) {
                timeStep();
            }
        }

        // do 3 sweeps
        void time3Sweep() {
            for (unsigned int i=0; i<3; ++i) {
                timeSweep();
            }
        }

        // chi: magnetic susceptibility, cap: heat capacity
        // does sweeps and calculate mean values when its done
        void doSweeps(const unsigned int sweeps, double& meanEnergy, double& meanMagnetization, double& chi, double& cap) {
            std::vector<double> magnetization(sweeps, 0);
            std::vector<double> energy(sweeps, 0);

            for (unsigned int i=0; i<sweeps; ++i) {
                magnetization[i] = computeMagnetization();
                energy[i] = computeEnergyOfSystem();

                // let the system evolve for 3 sweeps
                time3Sweep();
            }
            // for chi we have to compute the mean of the (magnetization²) and the (mean of the magnetization)²
            chi = beta_ * pow(latticeLength_, 3) * (computeSquaredMean(magnetization) - computeMeanSquared(magnetization));
            // for cap we have to compute the mean of the (energy)² and the (mean of the energy)²
            cap = pow(beta_, 2) / pow(latticeLength_, 3) * (computeSquaredMean(energy) - computeMeanSquared(energy));

            meanEnergy = computeMean(energy);
            meanMagnetization = computeMean(magnetization);
            return;
        }

        double computeMean(std::vector<double>& nums) {
            double sum = 0;
            for (unsigned i=0; i<nums.size(); ++i) {
                sum += nums[i];
            }
            return sum/nums.size();
        }

        double computeSquaredMean(std::vector<double>& nums) {
            double sum = 0;
            for (unsigned i=0; i<nums.size(); ++i) {
                sum += pow(nums[i], 2);
            }
            return sum/nums.size();
        }

        double computeMeanSquared(std::vector<double>& nums) {
            double sum = 0;
            for (unsigned i=0; i<nums.size(); ++i) {
                sum += nums[i];
            }
            return pow(sum/nums.size(), 2);
        }
        
        struct coordinate {
            unsigned int x,y,z;
            coordinate(unsigned int x, unsigned int y, unsigned int z) {
                this->x = x; this->y = y; this->z = z;
            }
        };
        
        void addNodeToClusterAndFlipSpinIfSameStartingSpin(
                std::vector<coordinate> &stack,
                int startingNodeSpin,
                double prob,
                unsigned int x,
                unsigned int y,
                unsigned int z)
        {
            if (getSpin(x,y,z) == startingNodeSpin) {
                if (real_d_(mt_) < prob) {
                    stack.push_back(coordinate(x,y,z));
                    flipSpin(x, y, z);
                }
            }
        }
        
        void doWolffStep() {
            int startingNodeSpin;
            double prob;
            unsigned int x,y,z;
            std::vector<coordinate> stack;

            // probability to connect nodes
            prob = 1-exp(-2*beta_*J_);
            
            // choose random node
            x=int_i_x_(mt_); y=int_i_y_(mt_); z=int_i_z_(mt_);
            
            // save its spin
            startingNodeSpin = getSpin(x,y,z);
            
            // add current node coordinates to stack at position 0
            stack.push_back(coordinate(x,y,z));
            
            // while stack not empty
            while (!stack.empty()) {
                // take last node coordinates from stack
                coordinate currentCoordinates = stack.back();
                stack.pop_back();

                // check all neighbours of node stack[i] with same spin state,
                  // add them to stack with probability p
                    
                // x+1
                addNodeToClusterAndFlipSpinIfSameStartingSpin(stack, startingNodeSpin, prob,
                        (currentCoordinates.x+1)%latticeLength_, currentCoordinates.y, currentCoordinates.z);

                // x-1
                addNodeToClusterAndFlipSpinIfSameStartingSpin(stack, startingNodeSpin, prob,
                        (currentCoordinates.x+latticeLength_-1)%latticeLength_, currentCoordinates.y, currentCoordinates.z);

                
                // y+1
                addNodeToClusterAndFlipSpinIfSameStartingSpin(stack, startingNodeSpin, prob,
                        currentCoordinates.x, (currentCoordinates.y+1)%latticeLength_, currentCoordinates.z);

                // y-1
                addNodeToClusterAndFlipSpinIfSameStartingSpin(stack, startingNodeSpin, prob,
                        currentCoordinates.x, (currentCoordinates.y+latticeLength_-1)%latticeLength_, currentCoordinates.z);

                
                // z+1
                addNodeToClusterAndFlipSpinIfSameStartingSpin(stack, startingNodeSpin, prob,
                        currentCoordinates.x, currentCoordinates.y, (currentCoordinates.z+1)%latticeLength_);

                // z-1
                addNodeToClusterAndFlipSpinIfSameStartingSpin(stack, startingNodeSpin, prob,
                        currentCoordinates.x, currentCoordinates.y, (currentCoordinates.z+latticeLength_-1)%latticeLength_);
            }
        }
};

void plot1() {
    IsingLattice myLattice(6, 10);
    for (unsigned int i=0; i<127; ++i) {
        std::cout << i << " " << myLattice.computeEnergyOfSystem() << " " << myLattice.computeMagnetization() << std::endl;
        myLattice.time3Sweep();
    }
}

void plot2() {
    double meanEnergy, meanMagnetization, chi, cap;
    double temperature = 3;

    for (int i=0;i<201;++i) {
        meanEnergy = meanMagnetization = chi = cap = 0;
        IsingLattice myLattice(temperature, 10);

        // let the system relax into equilibirum for 3 sweeps
        myLattice.time3Sweep();

        // do 2k*3 sweeps and calculate the values
        myLattice.doSweeps(2000*3, meanEnergy, meanMagnetization, chi, cap);
        std::cout << temperature << " " << meanEnergy << " " << meanMagnetization << " " << chi << " " << cap << std::endl;
        temperature += 0.01;
    }
}

// finite size scaling: measure ChiMax = chi(Tc) for different system sizes to get an estimate for gamma/nu (lecture slide 72)
// system sizes 3 - 20, delta 1
void measureChiMax() {
    const unsigned int startSize = 3;
    const unsigned int endSize = 20;

    double meanEnergy, meanMagnetization, chi, cap;

    const double Tc = 4.51;
    for (unsigned int i=startSize; i<=endSize; ++i) {
        std::cerr << "measuring for size " << i << std::endl;
        IsingLattice myLattice(Tc, i);
        // let the system relax into equilibirum for 3 sweeps
        myLattice.time3Sweep();

        // do 300*3 sweeps and calculate the values
        myLattice.doSweeps(300*3, meanEnergy, meanMagnetization, chi, cap);

        std::cout << i << " " << chi << std::endl;
    }
}

// measure magnetic susceptibility for different system sizes and temperatures
// temperature range 4.0 - 4.8 K, delta 0.01
// system sizes 3 - 20, delta 1
void finiteSizeScaling() {
    const unsigned int startSize = 3;
    const unsigned int endSize = 20;

    double meanEnergy, meanMagnetization, chi, cap;

    for (unsigned int i=startSize; i<=endSize; ++i) {
        std::cerr << std::endl << "measuring for size: " << i << std::endl << std::endl;

        #pragma omp parallel for
        for (int j=0;j<80;++j) {
            std::cerr << "   measuring for temperature: " << 4.0+j*0.01 << std::endl << "   ";

            IsingLattice myLattice(4.0+j*0.01, i);
            // let the system relax into equilibirum for 3 sweeps
            myLattice.time3Sweep();

            // do 5000*3 sweeps and calculate the values
            myLattice.doSweeps(5000*3, meanEnergy, meanMagnetization, chi, cap);

            std::cout << i << " " << 4.0+j*0.01 << " " << chi << std::endl;
        }
    }
}

// generates data to a E(T) curve using the wolff algorithm
void wolff1() {
    unsigned int systemSize = 10;
    unsigned int numWolffSteps = 10;
    unsigned int numMeasurementValues = 1e6;
    
//    unsigned int numWolffSteps = 1000;
//    unsigned int numMeasurementValues = 2000;
//    stderr: 200

//    unsigned int numWolffSteps = 10;
//    unsigned int numMeasurementValues = 1e6;
//    stderr: 280 --> shit

    std::vector<double> energy(numMeasurementValues, 0);
    
//    #pragma omp parallel for
//    for (int j=0;j<80;++j) {
    for (int j=0;j<8;++j) {
        std::cerr << "    measuring for temperature: " << 4.0+j*0.01 << std::endl << "   ";

        IsingLattice myLattice(4.0+j*0.01, systemSize);
        // numMeasurementValues measurement values will be averaged
        for (unsigned k=0; k<numMeasurementValues; ++k) {
            // do numWolffSteps wolff steps
            for (long l=0; l<numWolffSteps; ++l) {
                myLattice.doWolffStep();
            }
            // measure energy
            energy[k] = myLattice.computeEnergyOfSystem();
        }
        
        // calculate mean and standard deviation (http://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos)
        double sum = std::accumulate(energy.begin(), energy.end(), 0.0);
        double mean = sum / energy.size();
        std::vector<double> diff(energy.size());
        std::transform(energy.begin(), energy.end(), diff.begin(),
                       std::bind2nd(std::minus<double>(), mean));
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / energy.size());
        
        // average energy and print result
        std::cout << 4.0+j*0.01 << " " << mean << " " << stdev << std::endl;
    }
}

// to plot data fast:
// ./ising_model | gnuplot -p -e 'plot "< cat /proc/$$/fd/0"'

// add 'gomp' to linker settings in codeblocks
// add '-fopenmp' to additional compiler options in netbeans
int main()
{

//    plot1();
//    plot2();
//    measureChiMax();
//    finiteSizeScaling();
    wolff1();

    return 0;
}
